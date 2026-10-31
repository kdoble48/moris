/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_WRK_perform_refinement.cpp
 *
 */

#include "cl_WRK_Performer.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_WRK_perform_refinement.hpp"
#include "HMR_Globals.hpp"
#include "cl_HMR_Background_Element_Base.hpp"
#include <unordered_map>
#include <utility>

// Logging package
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

#include "cl_MTK_Field.hpp"

#include "cl_Parameter_List.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_HMR_Element.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_GEN_Cell_Feature_Detector.hpp"

#include <algorithm>
#include <limits>

namespace moris::wrk
{
    namespace
    {
        /**
         * Corner coordinates (exodus order) of the axis-aligned box spanned by an HMR Lagrange element, and the
         * element's Lagrange order (4 / 8 nodes: 1, more: 2).
         */
        Matrix< DDRMat >
        element_corners( const hmr::Element* aElement, uint& aLagrangeOrder )
        {
            Vector< mtk::Vertex* > tVertices = aElement->get_vertex_pointers();
            Matrix< DDRMat >       tFirst    = tVertices( 0 )->get_coords();
            uint                   tNumDims  = tFirst.numel();

            Matrix< DDRMat > tLower( 1, tNumDims, std::numeric_limits< real >::max() );
            Matrix< DDRMat > tUpper( 1, tNumDims, -std::numeric_limits< real >::max() );
            for ( const mtk::Vertex* iVertex : tVertices )
            {
                Matrix< DDRMat > tCoords = iVertex->get_coords();
                for ( uint iDim = 0; iDim < tNumDims; iDim++ )
                {
                    tLower( iDim ) = std::min( tLower( iDim ), tCoords( iDim ) );
                    tUpper( iDim ) = std::max( tUpper( iDim ), tCoords( iDim ) );
                }
            }
            aLagrangeOrder = tVertices.size() == ( 1u << tNumDims ) ? 1 : 2;

            uint             tNumCorners = 1u << tNumDims;
            Matrix< DDRMat > tCorners( tNumCorners, tNumDims );
            for ( uint iCorner = 0; iCorner < tNumCorners; iCorner++ )
            {
                // exodus order: (0,0), (1,0), (1,1), (0,1) in x-y, bottom face before top face
                uint tInFace = iCorner % 4;
                bool tHigh[ 3 ] = { tInFace == 1 or tInFace == 2, tInFace == 2 or tInFace == 3, iCorner >= 4 };
                for ( uint iDim = 0; iDim < tNumDims; iDim++ )
                {
                    tCorners( iCorner, iDim ) = tHigh[ iDim ] ? tUpper( iDim ) : tLower( iDim );
                }
            }
            return tCorners;
        }

        uint
        max_entry( const Vector< uint >& aVector )
        {
            uint tMax = 0;
            for ( uint iEntry : aVector )
            {
                tMax = std::max( tMax, iEntry );
            }
            return tMax;
        }
    }    // namespace


    Refinement_Mini_Performer::Refinement_Mini_Performer(
            Parameter_List&               aParameterlist,
            std::shared_ptr< Library_IO > aLibrary )
            : mLibrary( std::move( aLibrary ) )
    {
        // set field names
        Vector< std::string > tFieldNames;
        string_to_vector(
                aParameterlist.get< std::string >( "field_names" ),
                mParameters.mFieldNames );

        if ( !aParameterlist.get< std::string >( "levels_of_refinement" ).empty() )
        {
            // set refinement level
            string_to_cell_mat(
                    aParameterlist.get< std::string >( "levels_of_refinement" ),
                    mParameters.mRefinementLevel );
        }

        // set refinementpattern
        string_to_cell_mat(
                aParameterlist.get< std::string >( "refinement_pattern" ),
                mParameters.mRefinementPattern );

        if ( !aParameterlist.get< std::string >( "remeshing_copy_old_pattern_to_pattern" ).empty() )
        {
            // set refinement pattern
            string_to_cell_mat(
                    aParameterlist.get< std::string >( "remeshing_copy_old_pattern_to_pattern" ),
                    mParameters.mRefinemenCopytPatternToPattern );
        }

        string_to_vector(
                aParameterlist.get< std::string >( "refinement_function_name" ),
                mParameters.mRefinementFunctionName );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Refinement_Mini_Performer::perform_refinement(
            Vector< std::shared_ptr< mtk::Field > >& aFields,
            const std::shared_ptr< hmr::HMR >&       aHMR )
    {
        Tracer tTracer( "WRK", "Refinement Mini Performer", "Perform refinement" );
        // create field name to index map
        moris::map< std::string, moris_index > tFieldNameToIndexMap;

        for ( uint Ik = 0; Ik < aFields.size(); Ik++ )
        {
            tFieldNameToIndexMap[ aFields( Ik )->get_label() ] = Ik;
        }

        Vector< moris_index >                          tPattern;
        Vector< Vector< std::string > >                tFieldNames;
        Vector< Vector< uint > >                       tRefinements;
        Vector< sint >                                 tMaxRefinementPerLevel;
        Vector< Vector< hmr::Refinement_Function_2 > > tRefinementFunctions;

        this->prepare_input_for_refinement(
                tPattern,
                tFieldNames,
                tRefinements,
                tMaxRefinementPerLevel,
                tRefinementFunctions );

        if ( ( mParameters.mRefinementFunctionName.size() > 0 ) && mLibrary != nullptr )
        {
            mParameters.mRefinementFunction = mLibrary->load_function< moris::hmr::Refinement_Function >( mParameters.mRefinementFunctionName( 0 ) );
        }

        for ( uint Ik = 0; Ik < tPattern.size(); Ik++ )
        {
            // get pattern
            uint tActivationPattern = tPattern( Ik );

            for ( uint Ii = 0; Ii < tMaxRefinementPerLevel.size(); Ii++ )
            {
                // Mesh has changed after first refinement and therefore has to be rebuilt
                // This will only work with analytic fields. When remeshing the L2 id done in the remesh mini performer
                if ( Ii > 0 )
                {
                    // Interpolation_Mesh_HMR * tInterpolationMesh = aHMR->create_interpolation_mesh( tOrder, tPattern, tPattern );
                    ////   mtk::Mesh_Pair

                    for ( uint Ia = 0; Ia < tFieldNames( Ik ).size(); Ia++ )
                    {
                        //  mtk::Field * tField = aFields( aFieldNameToIndexMap.find( tFieldNames( Ia ) ) );
                        //  tField->
                    }
                }

                for ( uint Ia = 0; Ia < tFieldNames( Ik ).size(); Ia++ )
                {
                    const std::shared_ptr< mtk::Field >& tField = aFields( tFieldNameToIndexMap.find( tFieldNames( Ik )( Ia ) ) );

                    uint tOrder = tField->get_lagrange_order();

                    const Matrix< DDRMat >& tFieldValues = tField->get_values();

                    // Put elements on queue and set flag for refinement
                    aHMR->based_on_field_put_elements_on_queue(
                            tFieldValues,
                            tActivationPattern,
                            tOrder,
                            mParameters.mRefinementFunction );
                }
            }

            aHMR->perform_refinement( tActivationPattern );
            aHMR->update_refinement_pattern( tActivationPattern );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Refinement_Mini_Performer::perform_refinement_2(
            Vector< std::shared_ptr< mtk::Field > >& aFields,
            const std::shared_ptr< hmr::HMR >&       aHMR )
    {
        Tracer tTracer( "WRK", "Refinement Mini Performer", "Perform advanced refinement" );
        // create field name to index map
        moris::map< std::string, moris_index > tFieldNameToIndexMap;

        for ( uint Ik = 0; Ik < aFields.size(); Ik++ )
        {
            tFieldNameToIndexMap[ aFields( Ik )->get_label() ] = Ik;
        }

        Vector< moris_index >                          tPattern;
        Vector< Vector< std::string > >                tFieldNames;
        Vector< Vector< uint > >                       tRefinements;
        Vector< sint >                                 tMaxRefinementPerLevel;
        Vector< Vector< hmr::Refinement_Function_2 > > tRefinementFunctions;

        this->prepare_input_for_refinement(
                tPattern,
                tFieldNames,
                tRefinements,
                tMaxRefinementPerLevel,
                tRefinementFunctions );

        // get pattern to save
        uint tNumPatternToSave = mParameters.mRefinemenCopytPatternToPattern.size();

        moris::map< sint, sint > tReferencePattern;
        for ( uint Ik = 0; Ik < tNumPatternToSave; Ik++ )
        {
            tReferencePattern[ mParameters.mRefinemenCopytPatternToPattern( Ik )( 0 ) ] =
                    mParameters.mRefinemenCopytPatternToPattern( Ik )( 1 );
        }

        MORIS_ERROR( mParameters.mRefinementFunctionName.size() > 0, "Refinement function names not set" );
        MORIS_ERROR( mLibrary != nullptr, "mLibrary not set" );

        uint tWorkingPattern = aHMR->get_parameters()->get_working_pattern();    // mParameters->get_working_pattern();

        // create a map with ids
        std::unordered_map< moris_index, luint > tMap;

        Vector< hmr::Background_Element_Base* > tBGElements;
        aHMR->get_database()->get_background_mesh()->collect_all_elements( tBGElements );

        for ( uint Ib = 0; Ib < tBGElements.size(); Ib++ )
        {
            auto tHMRId_2 = tBGElements( Ib )->get_hmr_id();

            tMap[ tHMRId_2 ] = Ib;
        }

        for ( uint Ik = 0; Ik < tPattern.size(); Ik++ )
        {

            aHMR->get_database()->get_background_mesh()->reset_pattern( tWorkingPattern );

            // get pattern
            uint tActivationPattern = tPattern( Ik );

            sint tThisRefernecePattern = tReferencePattern.find( tActivationPattern );

            for ( uint Ia = 0; Ia < tFieldNames( Ik ).size(); Ia++ )
            {
                std::string tFieldName = tFieldNames( Ik )( Ia );

                // a field may be withheld because it lives on an outdated mesh (MB-033): refine on the others
                if ( not tFieldNameToIndexMap.key_exists( tFieldName ) )
                {
                    MORIS_LOG_WARNING( "perform_refinement_2(): refinement field '%s' not available; pattern %d is refined without it",
                            tFieldName.c_str(),
                            (int)tActivationPattern );
                    continue;
                }

                const std::shared_ptr< mtk::Field >& tField = aFields( tFieldNameToIndexMap.find( tFieldName ) );

                // get interpolation mesh from mesh pair
                moris::mtk::Mesh* tSourceMesh = tField->get_mesh_pair().get_interpolation_mesh();

                MORIS_ERROR( tSourceMesh != nullptr, "Source mesh of field %s is nullptr", tFieldName.c_str() );

                uint tNumElements = tSourceMesh->get_num_elems();

                for ( uint Ii = 0; Ii < tNumElements; Ii++ )
                {
                    mtk::Cell& tCell = tSourceMesh->get_mtk_cell( Ii );

                    uint tMaxLevel = MORIS_UINT_MAX;

                    enum hmr::ElementalRefienmentIndicator tRefIndicator =
                            tRefinementFunctions( Ik )( Ia )( &tCell, tField, tActivationPattern, tMaxLevel );

                    hmr::Background_Element_Base* tBGElementOld =
                            reinterpret_cast< hmr::Element* >( tCell.get_base_cell() )->get_background_element();

                    luint tHMRId = tBGElementOld->get_hmr_id();

                    auto tIter = tMap.find( tHMRId );

                    MORIS_ERROR( tIter != tMap.end(),
                            "perform_refinement_2(): cell of field '%s' has no background element in the current mesh "
                            "(field older than the last remesh?)",
                            tFieldName.c_str() );

                    moris_index tIndex = tIter->second;

                    hmr::Background_Element_Base* tBGElementNew = tBGElements( tIndex );

                    while ( !tBGElementNew->is_active( tThisRefernecePattern ) )
                    {
                        // jump to parent
                        tBGElementNew = tBGElementNew->get_parent();

                        MORIS_ERROR( tBGElementNew != nullptr,
                                "perform_refinement_2(): field '%s' lives on a mesh not covered by reference pattern %d "
                                "(field older than the last remesh?)",
                                tFieldName.c_str(),
                                (int)tThisRefernecePattern );
                    }

                    // in the case that the BG level is equal the max level, hold this element
                    if ( tRefIndicator == hmr::ElementalRefienmentIndicator::REFINE && tBGElementNew->get_level() >= tMaxLevel )
                    {
                        tRefIndicator = hmr::ElementalRefienmentIndicator::HOLD;
                    }

                    switch ( tRefIndicator )
                    {
                        case hmr::ElementalRefienmentIndicator::REFINE:
                        {
                            tBGElementNew->set_refined_flag( tWorkingPattern );
                            break;
                        }
                        case hmr::ElementalRefienmentIndicator::HOLD:
                        {
                            if ( tBGElementNew->get_level() > 0 )
                            {
                                tBGElementNew->get_parent()    //
                                        ->set_refined_flag( tWorkingPattern );
                            }

                            break;
                        }
                        case hmr::ElementalRefienmentIndicator::COARSEN:
                        {
                            if ( tBGElementNew->get_level() > 1 )
                            {
                                tBGElementNew->get_parent()    //
                                        ->get_parent()         //
                                        ->set_refined_flag( tWorkingPattern );
                            }

                            break;
                        }
                        case hmr::ElementalRefienmentIndicator::DROP:
                        {
                            // do nothing
                            break;
                        }
                        default:
                        {
                            MORIS_ERROR( false,
                                    "Refinement_Mini_Performer::perform_refinement_2 " );
                        }
                    }
                }
            }

            aHMR->perform_refinement_based_on_working_pattern( tActivationPattern );
            aHMR->update_refinement_pattern( tActivationPattern );

            aHMR->get_database()->get_background_mesh()->update_database();
            aHMR->get_database()->update_bspline_meshes();
            aHMR->get_database()->update_lagrange_meshes();
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Refinement_Mini_Performer::perform_refinement_based_on_working_pattern(
            Vector< std::shared_ptr< mtk::Field > >& aFields,
            const std::shared_ptr< hmr::HMR >&       aHMR )
    {
        // create field name to index map
        moris::map< std::string, moris_index > tFieldNameToIndexMap;

        for ( uint Ik = 0; Ik < aFields.size(); Ik++ )
        {
            tFieldNameToIndexMap[ aFields( Ik )->get_label() ] = Ik;
        }

        if ( ( mParameters.mRefinementFunctionName.size() > 0 ) && mLibrary != nullptr )
        {
            mParameters.mRefinementFunction = mLibrary->load_function< moris::hmr::Refinement_Function >( mParameters.mRefinementFunctionName( 0 ) );
        }

        MORIS_ERROR( mParameters.mRefinementFunction != nullptr, "Refinement function not set!" );

        for ( uint Ik = 0; Ik < mParameters.mRefinementPattern( 0 ).numel(); Ik++ )
        {
            // get pattern
            uint tActivationPattern = mParameters.mRefinementPattern( 0 )( Ik );

            for ( uint Ia = 0; Ia < mParameters.mFieldNames.size(); Ia++ )
            {
                const std::shared_ptr< mtk::Field >& tField = aFields( tFieldNameToIndexMap.find( mParameters.mFieldNames( Ia ) ) );

                uint tOrder = tField->get_lagrange_order();

                const Matrix< DDRMat >& tFieldValues = tField->get_values();

                // Put elements on queue and set flag for refinement
                aHMR->based_on_field_flag_elements_on_working_pattern(
                        tFieldValues,
                        4,    // tActivationPattern
                        tOrder,
                        mParameters.mRefinementFunction );
            }

            aHMR->perform_refinement_based_on_working_pattern( tActivationPattern );
            aHMR->update_refinement_pattern( tActivationPattern );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    uint
    Refinement_Mini_Performer::perform_refinement_low_level_elements(
            Vector< std::shared_ptr< mtk::Field > >& aFields,
            const std::shared_ptr< hmr::HMR >&       aHMR )
    {
        Tracer tTracer( "WRK", "Refinement Mini Performer", "Perform refinement of low level elements" );
        uint   tNumQueuedElements = 0;

        // create field name to index map
        moris::map< std::string, moris_index > tFieldNameToIndexMap;

        for ( uint Ik = 0; Ik < aFields.size(); Ik++ )
        {
            tFieldNameToIndexMap[ aFields( Ik )->get_label() ] = Ik;
        }

        Vector< moris_index >                          tPattern;
        Vector< Vector< std::string > >                tFieldNames;
        Vector< Vector< uint > >                       tRefinements;
        Vector< sint >                                 tMaxRefinementPerLevel;
        Vector< Vector< hmr::Refinement_Function_2 > > tRefinementFunctions;

        this->prepare_input_for_refinement(
                tPattern,
                tFieldNames,
                tRefinements,
                tMaxRefinementPerLevel,
                tRefinementFunctions );

        for ( uint Ik = 0; Ik < tPattern.size(); Ik++ )
        {
            // get pattern
            uint tActivationPattern = tPattern( Ik );

            for ( uint Ia = 0; Ia < tFieldNames( Ik ).size(); Ia++ )
            {
                const std::shared_ptr< mtk::Field >& tField = aFields( tFieldNameToIndexMap.find( tFieldNames( Ik )( Ia ) ) );

                uint tOrder = tField->get_lagrange_order();

                const Matrix< DDRMat >& tFieldValues = tField->get_values();

                // Put elements on queue and set flag for refinement
                tNumQueuedElements += aHMR->based_on_field_put_low_level_elements_on_queue(
                        tFieldValues,
                        tActivationPattern,
                        tOrder,
                        mParameters.mRefinementFunction );
            }

            aHMR->perform_refinement( tActivationPattern );
            aHMR->update_refinement_pattern( tActivationPattern );
        }

        return tNumQueuedElements;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Refinement_Mini_Performer::prepare_input_for_refinement(
            Vector< moris_index >&                          aPatternForRefinement,
            Vector< Vector< std::string > >&                aFieldsForRefinement,
            Vector< Vector< uint > >&                       aRefinements,
            Vector< sint >&                                 aMaxRefinementPerPattern,
            Vector< Vector< hmr::Refinement_Function_2 > >& aRefinementFunctions )
    {

        // produce unique list of pattern which will be refined
        for ( uint Ik = 0; Ik < mParameters.mRefinementPattern.size(); Ik++ )
        {
            for ( uint Ii = 0; Ii < mParameters.mRefinementPattern( Ik ).numel(); Ii++ )
            {
                aPatternForRefinement.push_back( mParameters.mRefinementPattern( Ik )( Ii ) );
            }
        }

        // Sort this created list
        std::sort( ( aPatternForRefinement.data() ).data(), ( aPatternForRefinement.data() ).data() + aPatternForRefinement.size() );

        // use std::unique and std::distance to create list containing all used dof types. This list is unique
        auto last = std::unique( ( aPatternForRefinement.data() ).data(), ( aPatternForRefinement.data() ).data() + aPatternForRefinement.size() );
        auto pos  = std::distance( ( aPatternForRefinement.data() ).data(), last );

        aPatternForRefinement.resize( pos );

        uint tNumberOfRefinementPattern = aPatternForRefinement.size();

        // resize
        aFieldsForRefinement.resize( tNumberOfRefinementPattern );
        aRefinements.resize( tNumberOfRefinementPattern );
        aRefinementFunctions.resize( tNumberOfRefinementPattern );
        aMaxRefinementPerPattern.resize( tNumberOfRefinementPattern, 0 );

        // create list with field pointers and refinements per pattern
        for ( uint Ik = 0; Ik < tNumberOfRefinementPattern; Ik++ )
        {
            moris_index tPattern = aPatternForRefinement( Ik );

            // loop over all fields and corresponding patterns. Find the pattern which corresponds to tPattern and put it in list.
            // This is kind of a brute force algorithm. however there will be only a few fields
            for ( uint Ii = 0; Ii < mParameters.mRefinementPattern.size(); Ii++ )
            {
                for ( uint Ia = 0; Ia < mParameters.mRefinementPattern( Ii ).numel(); Ia++ )
                {
                    if ( tPattern == mParameters.mRefinementPattern( Ii )( Ia ) )
                    {
                        aFieldsForRefinement( Ik ).push_back( mParameters.mFieldNames( Ii ) );
                        if ( mParameters.mRefinementLevel.size() > 0 )
                        {
                            aRefinements( Ik ).push_back( mParameters.mRefinementLevel( Ii )( Ia ) );
                            aMaxRefinementPerPattern( Ik ) = std::max( aMaxRefinementPerPattern( Ik ), mParameters.mRefinementLevel( Ii )( Ia ) );
                        }
                        if ( mParameters.mRefinementFunctionName.size() > 0 )
                        {
                            aRefinementFunctions( Ik ).push_back(
                                    mLibrary->load_function< moris::hmr::Refinement_Function_2 >( mParameters.mRefinementFunctionName( Ii ) ) );
                        }
                    }
                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Refinement_Mini_Performer::perform_refinement_old(
            const std::shared_ptr< hmr::HMR >&     aHMR,
            Vector< std::shared_ptr< Performer > > aPerformers,
            bool                                   aSimultaneous )
    {
        Tracer tTracer( "WRK", "Refinement Mini Performer", "Perform GEN refinement" );

        // sub-cell feature refinement (MB-084) starts unresolved for every field, on every cell
        mResolvedFeatureFields.clear();
        mFeatureMinLevel.clear();

        // get the maximum number of refinements that need to be performed and therefore refinement iterations
        moris::sint tMaxNumRefinements = this->get_max_refinement_level( aPerformers );

        // get number of Lagrange meshes
        moris::Matrix< DDSMat > tMeshIndices;
        uint                    tNumMeshes = 0;
        this->get_all_refinement_mesh_indices( //  list of all the mesh indices listed for refinement in GEN
                aPerformers,
                tMeshIndices,
                tNumMeshes );

        // get number of performers that request refinement
        moris_index tNumPerformers = aPerformers.size();

        // get access to the database and background mesh
        std::shared_ptr< hmr::Database > tDataBase       = aHMR->get_database();
        hmr::Background_Mesh_Base*       tBackgroundMesh = tDataBase->get_background_mesh();

        // loop over the Lagrange meshes which are listed for refinement in GEN
        for ( uint iMesh = 0; iMesh < tNumMeshes; iMesh++ )
        {
            // get the Lagrange mesh index
            moris_index tMeshIndex = tMeshIndices( iMesh );

            // iteratively refine meshes
            for ( moris_index iRefinement = 0; iRefinement < tMaxNumRefinements; iRefinement++ )
            {
                // Create an hmr::Mesh (which is a wrapper for the Lagrange mesh to interface via the mtk::Mesh)
                std::shared_ptr< hmr::Mesh > tMesh = aHMR->create_mesh( tMeshIndex );

                uint tLagrangeMeshPattern = tMesh->get_lagrange_mesh_pattern(); // get the pattern associated with the Lagrange mesh which was requested to be refined in GEN

                for ( moris_index iPerformer = 0; iPerformer < tNumPerformers; iPerformer++ )
                {
                    // find all elements active on the Lagrange mesh for which there is a sign-change in the LS-field
                    this->queue_single_refinement(
                            aHMR,
                            tMesh,
                            aPerformers( iPerformer ),
                            iRefinement,
                            tMeshIndex );
                }

                //  perform refinement on the geometric meshes (replacing elements with their children, where necessary)
                aHMR->perform_refinement( tLagrangeMeshPattern );
                aHMR->update_refinement_pattern( tLagrangeMeshPattern ); // update the B-spline and Lagrange meshes

                // TODO: ?
                tBackgroundMesh->update_database();
                tDataBase->update_bspline_meshes();
                tDataBase->update_lagrange_meshes();

            }    // end for: each refinement step

            // TODO: ?
            tBackgroundMesh->clear_refinement_queue();

        }    // end for: each pattern which is listed for local (geometric) refinement

        // refinement loop to ensure that all intersected elements are on same refinement level
        for ( uint iMesh = 0; iMesh < tNumMeshes; iMesh++ )
        {
            moris_index tMeshIndex = tMeshIndices( iMesh );

            bool tRefinedAllLowLevelElements = false;

            // FIXME: ...
            const uint tMaxLowLevelRefinementSteps = 20;

            // set to true for using low level element refinement
            for ( uint iLowLevelRefinement = 0; iLowLevelRefinement < tMaxLowLevelRefinementSteps; ++iLowLevelRefinement )
            {
                // Create mesh //FIXME
                std::shared_ptr< hmr::Mesh > tMesh = aHMR->create_mesh( tMeshIndex );

                uint tLagrangeMeshPattern = tMesh->get_lagrange_mesh_pattern();

                uint tNumQueuedElements = 0;

                for ( moris_index iPerformer = 0; iPerformer < tNumPerformers; iPerformer++ )
                {
                    // Queue refinement
                    tNumQueuedElements += this->queue_low_level_elements_for_refinement(
                            aHMR,
                            tMesh,
                            aPerformers( iPerformer ),
                            tMeshIndex );
                }

                // compute number of queued elements across all processors
                tNumQueuedElements = sum_all( tNumQueuedElements );

                MORIS_LOG_INFO( "Found %d elements for low level refinement.", tNumQueuedElements );

                if ( tNumQueuedElements == 0 )
                {
                    tRefinedAllLowLevelElements = true;
                    break;
                }

                // Perform refinement and update indices of elements
                aHMR->perform_refinement( tLagrangeMeshPattern );
                aHMR->update_refinement_pattern( tLagrangeMeshPattern );

                // FIXME should be removed such that loop is continued until all elements are refined

            }    // end for: each lower level refinement

            // check that all low level elements were refined
            MORIS_ERROR(
                    tRefinedAllLowLevelElements,
                    "Refinement_Mini_Performer::perform_refinement_old() - "
                    "could not refine all low level elements in %d steps.",
                    tMaxLowLevelRefinementSteps );

        }    // end for: each mesh // TODO: what kind of meshes?

        // put every cell that XTK sees as cut on one level (sub-cell feature refinement; D-205), then report whether
        // the sub-cell features are represented on the final mesh (MB-084)
        for ( uint iMesh = 0; iMesh < tNumMeshes; iMesh++ )
        {
            moris_index tMeshIndex = tMeshIndices( iMesh );

            this->equalize_cut_cell_levels( aHMR, aPerformers, tMeshIndex );

            for ( moris_index iPerformer = 0; iPerformer < tNumPerformers; iPerformer++ )
            {
                for ( uint iField = 0; iField < aPerformers( iPerformer )->get_num_refinement_fields(); iField++ )
                {
                    if ( aPerformers( iPerformer )->get_feature_refinement_samples( iField ) == 0 )
                    {
                        continue;
                    }
                    const Vector< uint >& tLagrangeMeshIndices = aPerformers( iPerformer )->get_refinement_mesh_indices( iField );
                    if ( std::find( tLagrangeMeshIndices.begin(), tLagrangeMeshIndices.end(), uint( tMeshIndex ) ) == tLagrangeMeshIndices.end() )
                    {
                        continue;
                    }
                    this->queue_feature_refinement( aHMR, aPerformers( iPerformer ), iField, tMaxNumRefinements, tMeshIndex, true );
                }
            }
        }

        // TODO: ?
        tDataBase->update_bspline_meshes();
        tDataBase->update_lagrange_meshes();

    }    // end function: Refinement_Mini_Performer::perform_refinement_old()

    //--------------------------------------------------------------------------------------------------------------

    void
    Refinement_Mini_Performer::queue_single_refinement(
            const std::shared_ptr< hmr::HMR >&  aHMR,
            const std::shared_ptr< hmr::Mesh >& aMesh,
            const std::shared_ptr< Performer >& aPerformer,
            uint                                aRefinementNumber,
            uint                                aMeshIndex )
    {
        // uint tLagrangeMeshPattern = aMesh->get_lagrange_mesh_pattern();

        // Loop over fields
        for ( uint iField = 0; iField < aPerformer->get_num_refinement_fields(); iField++ )
        {
            const Vector< uint >& tNumRefinements      = aPerformer->get_num_refinements( iField );
            const Vector< uint >& tLagrangeMeshIndices = aPerformer->get_refinement_mesh_indices( iField );

            // loop over tLagrangeMeshIndices // if aMeshIndex put in queue
            for ( uint iLagMesh = 0; iLagMesh < tLagrangeMeshIndices.size(); iLagMesh++ )
            {
                if ( tLagrangeMeshIndices( iLagMesh ) == aMeshIndex )
                {
                    // sub-cell feature refinement (MB-084): number_of_refinements is the most passes
                    if ( aPerformer->get_feature_refinement_samples( iField ) > 0 )
                    {
                        if ( tNumRefinements( iLagMesh ) > aRefinementNumber
                                and mResolvedFeatureFields.count( { aPerformer.get(), iField } ) == 0 )
                        {
                            this->queue_feature_refinement( aHMR, aPerformer, iField, aRefinementNumber, aMeshIndex, false );
                        }
                        continue;
                    }

                    if ( tNumRefinements( iLagMesh ) > aRefinementNumber )
                    {
                        // Loop over all nodes in the Lagrange mesh and get field values
                        Matrix< DDRMat > tFieldValues( aMesh->get_num_nodes(), 1 );
                        for ( uint tNodeIndex = 0; tNodeIndex < aMesh->get_num_nodes(); tNodeIndex++ )
                        {
                            Matrix< DDRMat > tNodeCoordinates = aMesh->get_node_coordinate( tNodeIndex );

                            // FIXME right now the geometric region is converted into a real value, but this should be handled by the refinement interface
                            tFieldValues( tNodeIndex ) = (real)aPerformer->get_geometric_region( iField, tNodeIndex, tNodeCoordinates );
                        }

                        // Put elements on queue and set flag for refinement
                        sint tRefinementFunctionIndex = aPerformer->get_refinement_function_index( iField, aRefinementNumber );
                        aHMR->based_on_field_put_elements_on_queue( tFieldValues, aMeshIndex, tRefinementFunctionIndex );
                    }
                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    uint
    Refinement_Mini_Performer::queue_low_level_elements_for_refinement(
            const std::shared_ptr< hmr::HMR >&  aHMR,
            const std::shared_ptr< hmr::Mesh >& aMesh,
            const std::shared_ptr< Performer >& aPerformer,
            uint                                aMeshIndex )
    {
        uint tNumElements = 0;

        // Loop over fields
        for ( uint iField = 0; iField < aPerformer->get_num_refinement_fields(); iField++ )
        {
            const Vector< uint >& tLagrangeMeshIndices = aPerformer->get_refinement_mesh_indices( iField );

            // loop over tLagrangeMesh indices // if aMeshIndex put in queue
            for ( uint iLagMesh = 0; iLagMesh < tLagrangeMeshIndices.size(); iLagMesh++ )
            {
                if ( tLagrangeMeshIndices( iLagMesh ) == aMeshIndex )
                {
                    // get the number of nodes on the mesh
                    uint tNumNodes = aMesh->get_num_nodes();

                    // Loop over nodes and get field values
                    Matrix< DDRMat > tFieldValues( tNumNodes, 1 );

                    for ( uint tNodeIndex = 0; tNodeIndex < tNumNodes; tNodeIndex++ )
                    {
                        Matrix< DDRMat > tNodeCoordinates = aMesh->get_node_coordinate( tNodeIndex );

                        // FIXME right now the geometric region is converted into a real value, but this should be handled by the refinement interface
                        tFieldValues( tNodeIndex ) = (real)aPerformer->get_geometric_region(
                                iField,
                                tNodeIndex,
                                tNodeCoordinates );
                    }

                    // Put elements on queue and set flag for refinement // FIXME: this is untested for a refinement function,
                    sint tRefineFuncIndex = aPerformer->get_refinement_function_index( iField, 0 );
                    tNumElements += aHMR->based_on_field_put_low_level_elements_on_queue(
                            tFieldValues,
                            aMeshIndex,
                            tRefineFuncIndex );
                }
            }
        }

        return tNumElements;
    }

    //--------------------------------------------------------------------------------------------------------------

    uint
    Refinement_Mini_Performer::queue_feature_refinement(
            const std::shared_ptr< hmr::HMR >&  aHMR,
            const std::shared_ptr< Performer >& aPerformer,
            uint                                aFieldIndex,
            uint                                aPass,
            uint                                aMeshIndex,
            bool                                aFinalCheck )
    {
        uint tSamples = aPerformer->get_feature_refinement_samples( aFieldIndex );

        // the active (leaf) cells of the Lagrange mesh
        Vector< hmr::Element* > tCandidates;
        aHMR->get_active_candidates_for_refinement( tCandidates, aMeshIndex );

        Vector< hmr::Element* > tQueue;
        tQueue.reserve( tCandidates.size() );
        uint tNumCrossing      = 0;
        uint tNumHidden        = 0;
        uint tNumInsignificant = 0;
        uint tMinLevel         = std::numeric_limits< uint >::max();
        uint tMaxLevel         = 0;

        // cells below this level were evaluated in an earlier pass without a crossing and were not refined since
        uint tMinLevel0 = mFeatureMinLevel.count( { aPerformer.get(), aFieldIndex } ) ? mFeatureMinLevel[ { aPerformer.get(), aFieldIndex } ] : 0;
        uint tNumSkipped = 0;

        for ( hmr::Element* iElement : tCandidates )
        {
            if ( uint( iElement->get_level() ) < tMinLevel0 )
            {
                tNumSkipped++;
                continue;
            }
            uint             tOrder   = 1;
            Matrix< DDRMat > tCorners = element_corners( iElement, tOrder );
            uint             tFlags   = aPerformer->detect_cell_features( aFieldIndex, tCorners, tOrder, tSamples );

            // every cell with a crossing is refined, also those XTK already cuts: XTK needs all cut cells on one level
            if ( tFlags & ( gen::CELL_FEATURE_CROSSING | gen::CELL_FEATURE_XTK_CUT ) )
            {
                tQueue.push_back( iElement );
                tNumCrossing++;
                tMinLevel = std::min( tMinLevel, uint( iElement->get_level() ) );
                tMaxLevel = std::max( tMaxLevel, uint( iElement->get_level() ) );
            }
            tNumHidden += ( tFlags & gen::CELL_FEATURE_HIDDEN ) ? 1 : 0;
            tNumInsignificant += ( tFlags & gen::CELL_FEATURE_INSIGNIFICANT ) ? 1 : 0;
        }

        tNumCrossing      = sum_all( tNumCrossing );
        tNumHidden        = sum_all( tNumHidden );
        tNumInsignificant = sum_all( tNumInsignificant );
        tMinLevel         = min_all( tMinLevel );
        tMaxLevel         = max_all( tMaxLevel );

        if ( aFinalCheck )
        {
            if ( mResolvedFeatureFields.count( { aPerformer.get(), aFieldIndex } ) > 0 or tNumHidden == 0 )
            {
                MORIS_LOG_INFO( "GEN feature refinement - geometry %u: all sub-cell features represented on the final mesh "
                                "(%u cut cells on levels %u to %u).",
                        aFieldIndex, tNumCrossing, tNumCrossing ? tMinLevel : 0, tMaxLevel );
                if ( tNumInsignificant > 0 )
                {
                    MORIS_LOG_WARNING( "GEN feature refinement - geometry %u: %u cells hold features shallower than the depth "
                                       "tolerance (feature_refinement_depth_tolerance x cell size). They are UNRESOLVED: not "
                                       "refined and not represented by the cut (below the stated resolution, not certified).",
                            aFieldIndex, tNumInsignificant );
                }
            }
            else
            {
                MORIS_LOG_WARNING( "GEN feature refinement - geometry %u: UNRESOLVED. %u cells still hide sub-cell features after "
                                   "%u refinement passes (number_of_refinements reached; MB-084). The cut does not represent "
                                   "the geometry in these cells.",
                        aFieldIndex, tNumHidden, aPass );
            }
            return tNumHidden;
        }

        MORIS_LOG_INFO( "GEN feature refinement - geometry %u, pass %u: %u cells with a crossing (levels %u to %u), %u hide a "
                        "feature the cut misses, %u hold only features below the depth tolerance.",
                aFieldIndex, aPass, tNumCrossing, tNumCrossing ? tMinLevel : 0, tMaxLevel, tNumHidden, tNumInsignificant );

        if ( tNumHidden == 0 )
        {
            mResolvedFeatureFields.insert( { aPerformer.get(), aFieldIndex } );
            MORIS_LOG_INFO( "GEN feature refinement - geometry %u: resolved after %u passes.", aFieldIndex, aPass );
            return 0;
        }

        aHMR->put_elements_on_refinement_queue( tQueue );

        // every cell with a crossing is refined: crossings can only remain on the new level
        mFeatureMinLevel[ { aPerformer.get(), aFieldIndex } ] = tMaxLevel + 1;
        MORIS_LOG_INFO( "GEN feature refinement - geometry %u, pass %u: %u cells below level %u skipped (evaluated before).",
                aFieldIndex, aPass, sum_all( tNumSkipped ), tMinLevel0 );
        return tNumHidden;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Refinement_Mini_Performer::equalize_cut_cell_levels(
            const std::shared_ptr< hmr::HMR >&            aHMR,
            const Vector< std::shared_ptr< Performer > >& aPerformers,
            uint                                          aMeshIndex )
    {
        // fields refined on this mesh that XTK can see cut inside cells whose nodes share one region
        Vector< std::pair< uint, uint > > tFields;    // (performer, field)
        bool                              tAnyFeature = false;
        for ( uint iPerformer = 0; iPerformer < aPerformers.size(); iPerformer++ )
        {
            for ( uint iField = 0; iField < aPerformers( iPerformer )->get_num_refinement_fields(); iField++ )
            {
                const Vector< uint >& tMeshIndices    = aPerformers( iPerformer )->get_refinement_mesh_indices( iField );
                const Vector< uint >& tNumRefinements = aPerformers( iPerformer )->get_num_refinements( iField );
                bool                  tOnThisMesh     = std::find( tMeshIndices.begin(), tMeshIndices.end(), aMeshIndex ) != tMeshIndices.end();
                bool                  tFeature        = aPerformers( iPerformer )->get_feature_refinement_samples( iField ) > 0;
                // default criterion: only analytic fields on the decomposition mesh (index 0), which can be cut inside a
                // cell whose nodes HMR saw uniform
                bool tDefault = aMeshIndex == 0 and max_entry( tNumRefinements ) > 0
                            and aPerformers( iPerformer )->can_change_region_inside_cell( iField );
                if ( tOnThisMesh and ( tFeature or tDefault ) )
                {
                    tFields.push_back( { iPerformer, iField } );
                    tAnyFeature = tAnyFeature or tFeature;
                }
            }
        }
        if ( tFields.size() == 0 )
        {
            return;
        }

        // XTK checks the levels of the cells cut by any geometry: once this mesh is refined, also test the other
        // geometries that XTK can see cut inside a cell (exact fields; e.g. fixed pads next to a refined design)
        for ( uint iPerformer = 0; iPerformer < aPerformers.size(); iPerformer++ )
        {
            for ( uint iField = 0; iField < aPerformers( iPerformer )->get_num_refinement_fields(); iField++ )
            {
                bool tListed = false;
                for ( const auto& iListed : tFields )
                {
                    tListed = tListed or ( iListed.first == iPerformer and iListed.second == iField );
                }
                if ( not tListed and aMeshIndex == 0 and aPerformers( iPerformer )->can_change_region_inside_cell( iField ) )
                {
                    tFields.push_back( { iPerformer, iField } );
                }
            }
        }

        const uint tMaxSteps = 20;
        for ( uint iStep = 0; iStep < tMaxSteps; iStep++ )
        {
            std::shared_ptr< hmr::Mesh > tMesh                = aHMR->create_mesh( aMeshIndex );
            uint                         tLagrangeMeshPattern = tMesh->get_lagrange_mesh_pattern();

            Vector< hmr::Element* > tCandidates;
            aHMR->get_active_candidates_for_refinement( tCandidates, aMeshIndex );

            // cells XTK will see as cut (its own samples: lattice with 2 sub-intervals)
            Vector< hmr::Element* > tCutCells;
            sint                    tMaxLevel = -1;
            for ( hmr::Element* iElement : tCandidates )
            {
                uint             tOrder   = 1;
                Matrix< DDRMat > tCorners = element_corners( iElement, tOrder );
                for ( const auto& iField : tFields )
                {
                    // a feature-refined field has no crossing below its minimum level (see queue_feature_refinement)
                    auto tKey = std::make_pair( static_cast< const Performer* >( aPerformers( iField.first ).get() ), iField.second );
                    if ( mFeatureMinLevel.count( tKey ) and uint( iElement->get_level() ) < mFeatureMinLevel[ tKey ] )
                    {
                        continue;
                    }
                    if ( aPerformers( iField.first )->detect_cell_features( iField.second, tCorners, tOrder, 2 ) & gen::CELL_FEATURE_XTK_CUT )
                    {
                        tCutCells.push_back( iElement );
                        tMaxLevel = std::max( tMaxLevel, sint( iElement->get_level() ) );
                        break;
                    }
                }
            }
            tMaxLevel = max_all( tMaxLevel );

            Vector< hmr::Element* > tBelow;
            for ( hmr::Element* iElement : tCutCells )
            {
                if ( sint( iElement->get_level() ) < tMaxLevel )
                {
                    tBelow.push_back( iElement );
                }
            }
            uint tNumBelow = sum_all( uint( tBelow.size() ) );
            if ( tNumBelow == 0 )
            {
                MORIS_LOG_INFO( "GEN refinement: %u cells that XTK will see as cut, all on level %d (%u geometries tested).",
                        sum_all( uint( tCutCells.size() ) ), tMaxLevel, uint( tFields.size() ) );
                return;
            }

            MORIS_LOG_INFO( "GEN refinement: %u cells that XTK sees as cut lie below level %d; refining them so that all cut "
                            "cells share one level (%s).",
                    tNumBelow, tMaxLevel, tAnyFeature ? "sub-cell feature refinement" : "D-205, analytic field cut inside cells" );

            aHMR->put_elements_on_refinement_queue( tBelow );
            aHMR->perform_refinement( tLagrangeMeshPattern );
            aHMR->update_refinement_pattern( tLagrangeMeshPattern );
            aHMR->get_database()->get_background_mesh()->update_database();
            aHMR->get_database()->update_bspline_meshes();
            aHMR->get_database()->update_lagrange_meshes();
        }

        MORIS_LOG_WARNING( "GEN refinement: cut cells still lie on different levels after %u equalisation steps; XTK will "
                           "refuse the decomposition.",
                tMaxSteps );
    }

    //--------------------------------------------------------------------------------------------------------------

    uint
    Refinement_Mini_Performer::get_max_refinement_level( const Vector< std::shared_ptr< Performer > >& aPerformers )
    {
        uint tMaxNumRefinements = 0;

        for ( uint iPerformer = 0; iPerformer < aPerformers.size(); iPerformer++ )
        {
            uint tNumFields = aPerformers( iPerformer )->get_num_refinement_fields();

            // Loop over fields
            for ( uint iField = 0; iField < tNumFields; iField++ )
            {
                const Vector< uint >& tNumRefinements = aPerformers( iPerformer )->get_num_refinements( iField );

                for ( uint iRefinement = 0; iRefinement < tNumRefinements.size(); iRefinement++ )
                {
                    tMaxNumRefinements = std::max( tMaxNumRefinements, tNumRefinements( iRefinement ) );
                }
            }
        }

        return tMaxNumRefinements;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Refinement_Mini_Performer::get_all_refinement_mesh_indices(
            const Vector< std::shared_ptr< Performer > >& aPerformers,
            moris::Matrix< DDSMat >&                      aAllPatternMap,
            moris::uint&                                  aNumPattern )
    {
        moris::Matrix< DDSMat > tCombinedPattern( hmr::gNumberOfPatterns, 1, -1 );

        aNumPattern = 0;

        for ( uint iPerformer = 0; iPerformer < aPerformers.size(); iPerformer++ )
        {
            uint tNumFields = aPerformers( iPerformer )->get_num_refinement_fields();

            // Loop over fields
            for ( uint iRefineField = 0; iRefineField < tNumFields; iRefineField++ )
            {
                const Vector< uint >& tRefinementMeshIndices = aPerformers( iPerformer )->get_refinement_mesh_indices( iRefineField );

                for ( uint iMesh = 0; iMesh < tRefinementMeshIndices.size(); iMesh++ )
                {
                    sint tRefinementMeshIndex = tRefinementMeshIndices( iMesh );

                    if ( tCombinedPattern( tRefinementMeshIndex ) == -1 )
                    {
                        tCombinedPattern( tRefinementMeshIndex ) = 1;

                        aNumPattern++;
                    }
                }
            }
        }

        aAllPatternMap.set_size( aNumPattern, 1, -1 );

        uint tCounter = 0;

        for ( uint iPattern = 0; iPattern < tCombinedPattern.numel(); iPattern++ )
        {
            if ( tCombinedPattern( iPattern ) == 1 )
            {
                aAllPatternMap( tCounter++ ) = iPattern;
            }
        }

    }    // end function: Refinement_Mini_Performer::get_all_refinement_mesh_indices()

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::wrk
