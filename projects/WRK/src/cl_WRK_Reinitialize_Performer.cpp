/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_WRK_Reinitialize_Performer.cpp
 *
 */

#include "cl_WRK_Reinitialize_Performer.hpp"

#include <memory>
#include <sstream>
#include <set>
#include "cl_Matrix.hpp"
#include "cl_Tracer.hpp"
#include "cl_Logger.hpp"
#include "cl_MTK_Mesh_Pair.hpp"

#include "cl_MTK_Field.hpp"
#include "cl_MTK_Field_Discrete.hpp"
#include "cl_MTK_Mapper.hpp"
#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_MTK_Reader_Exodus.hpp"

#include "cl_MSI_Dof_Type_Enums.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_HMR_File.hpp"
#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR_Mesh_Integration.hpp"
#include "HMR_Globals.hpp"

#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_MDL_Model.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_MSI_Solver_Interface.hpp"

#include "cl_SOL_Dist_Vector.hpp"
#include "cl_SOL_Dist_Map.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Warehouse.hpp"
#include "cl_MSI_Equation_Model.hpp"
#include "fn_sort.hpp"

// XTK includes for interface-based SDF
#include "cl_XTK_Cut_Integration_Mesh.hpp"

// SDF includes for distance computation
#include "cl_SDF_From_Interface.hpp"

namespace moris::wrk
{
    //------------------------------------------------------------------------------

    Reinitialize_Performer::Reinitialize_Performer( const std::shared_ptr< Library_IO >& aLibrary )
            : mLibrary( aLibrary )
    {
        // get the parameter lists
        Module_Parameter_Lists tMORISParameterList = aLibrary->get_parameters_for_module( Module_Type::MORISGENERAL );
        Module_Parameter_Lists tMSIParameterList   = aLibrary->get_parameters_for_module( Module_Type::MSI );
        Module_Parameter_Lists tOPTParameterList   = aLibrary->get_parameters_for_module( Module_Type::OPT );

        // First determine the SDF reinit mode from OPT params
        std::string tReinitModeStr = tOPTParameterList( 0 )( 0 ).get< std::string >( "sdf_reinit_mode" );
        if ( tReinitModeStr == "interface" )
        {
            mReinitMode = ReinitMode::INTERFACE_SDF;
        }
        else
        {
            mReinitMode = ReinitMode::FIELD_REMAP;  // Default
        }

        // Material phase for sign determination - read from OPT params
        mMaterialPhase = tOPTParameterList( 0 )( 0 ).get< sint >( "sdf_reinit_material_phase" );

        // Phases to exclude from interface facets (e.g., external material zone)
        std::string tExcludePhasesStr = tOPTParameterList( 0 )( 0 ).get< std::string >( "sdf_exclude_phases" );
        if ( !tExcludePhasesStr.empty() )
        {
            // Parse comma-separated list of phase indices
            std::stringstream ss( tExcludePhasesStr );
            std::string token;
            while ( std::getline( ss, token, ',' ) )
            {
                // Trim whitespace
                token.erase( 0, token.find_first_not_of( " \t" ) );
                token.erase( token.find_last_not_of( " \t" ) + 1 );
                if ( !token.empty() )
                {
                    mExcludePhases.push_back( std::stoi( token ) );
                }
            }
        }

        // Number of Laplacian smoothing iterations
        mSmoothingIterations = tOPTParameterList( 0 )( 0 ).get< sint >( "sdf_smoothing_iterations" );
        
        // Laplacian smoothing strength (lambda)
        mSmoothingLambda = tOPTParameterList( 0 )( 0 ).get< real >( "sdf_smoothing_lambda" );

        // Initialize the design as an SDF on the first optimization iteration
        mInitializeAtStart = tOPTParameterList( 0 )( 0 ).get< bool >( "sdf_initialize_at_start" );

        // Wean-off tolerance on the reinit's relative change to the field
        mWeanTol = tOPTParameterList( 0 )( 0 ).get< real >( "sdf_reinit_wean_tol" );

        // Redistance every iteration in-place (no GCMMA restart); smooth every reinit_freq
        mReinitEveryStep = tOPTParameterList( 0 )( 0 ).get< bool >( "sdf_reinit_every_step" );

        mADVFiledName = tMORISParameterList( 2 )( 0 ).get< std::string >( "adv_field" );
        
        // Reinitialization frequency from mapping params
        mReinitializationFrequency = tMORISParameterList( 2 )( 0 ).get< sint >( "reinitialization_frequency" );

        // get the mesh output info
        mOutputMeshFile = tMORISParameterList( 2 )( 0 ).get< std::string >( "output_mesh_file" );
        mTimeOffset     = tMORISParameterList( 2 )( 0 ).get< real >( "time_offset" );

        // dof_type is only needed for FIELD_REMAP mode  
        if ( mReinitMode == ReinitMode::FIELD_REMAP )
        {
            std::string tDofTypeStr = tMORISParameterList( 2 )( 0 ).get< std::string >( "dof_type" );
            
            if ( !tDofTypeStr.empty() )
            {
                mAdofMeshIndex = tMSIParameterList( 0 )( 0 ).get< moris::sint >( tDofTypeStr );

                // get msi string to dof type map
                moris::map< std::string, MSI::Dof_Type > tMSIDofTypeMap =
                        moris::MSI::get_msi_dof_type_map();

                // get the quantity dof type from parameter list
                string_to_vector( tDofTypeStr, mDofTypes, tMSIDofTypeMap );
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    Reinitialize_Performer::perform(
            Vector< std::shared_ptr< hmr::HMR > >&             aHMRPerformers,
            Vector< std::shared_ptr< gen::Geometry_Engine > >& aGENPerformer,
            Vector< std::shared_ptr< mtk::Mesh_Manager > >&    aMTKPerformer,
            Vector< std::shared_ptr< mdl::Model > >            aMDLPerformer )
    {
        // Tracer to trace the time
        Tracer tTracer( "WRK", "Reinitialize ADVs", "Perform Reinitialize" );

        // initialize and populate the fields
        Vector< std::shared_ptr< mtk::Field > > tGENFields;
        tGENFields.append( aGENPerformer( 0 )->get_mtk_fields() );

        // find the index of the desired adv field that will be reinitialized
        auto itr = std::find_if( tGENFields.begin(), tGENFields.end(), [ & ]( std::shared_ptr< mtk::Field > const & aFiled )    //
                { return aFiled->get_label() == mADVFiledName; } );

        // find the index of the adv field
        moris_index tADVFieldIndex = std::distance( tGENFields.begin(), itr );

        // get the the adv discretization mesh index
        uint tDiscretizationMeshIndex = ( *itr )->get_discretization_mesh_index();

        // get interpolation mesh from mesh pair
        moris::mtk::Mesh* tTargetMesh = ( *itr )->get_mesh_pair().get_interpolation_mesh();

        // get the solution field and get a matrix of the solutions
        // generate a cell containing the indices of the bspline coefficients
        // since indices are consecutive and they start from 0
        Vector< moris_index > tLocalCoeffIndices( tTargetMesh->get_num_entities( mtk::EntityRank::BSPLINE ) );
        std::iota( tLocalCoeffIndices.begin(), tLocalCoeffIndices.end(), 0 );

        moris::sol::Dist_Vector* tPartialSolutionVector = aMDLPerformer( 0 )->get_solver_interface()->get_solution_vector( mDofTypes, tLocalCoeffIndices );
        tPartialSolutionVector->extract_copy( mCoefficients );

        // delete the pointer as it is not needed anymore
        delete tPartialSolutionVector;

        // create field object for this mesh ,the discretization index is zero as there is only one discretization in the newly constructed IP mesh
        std::shared_ptr< mtk::Field_Discrete > tFieldSource = std::make_shared< mtk::Field_Discrete >( aMTKPerformer( 0 )->get_mesh_pair( 0 ), mAdofMeshIndex );

        // unlock fields and set the coeff
        tFieldSource->unlock_field();
        tFieldSource->set_coefficients( mCoefficients );

        // compute the nodal values based on the coeff
        tFieldSource->compute_nodal_values();

        // create field object for this mesh ,the discretization index is zero as there is only one discretization in the newly constructed IP mesh
        std::shared_ptr< mtk::Field_Discrete > tFieldTarget = std::make_shared< mtk::Field_Discrete >( ( *itr )->get_mesh_pair(), tDiscretizationMeshIndex );
        tFieldTarget->set_label( mADVFiledName );

        // set the nodal values
        tFieldTarget->unlock_field();
        tFieldTarget->set_values( tFieldSource->get_values() );

        // invoke the mapper and map to the target field
        mtk::Mapper tMapper;
        tFieldTarget->unlock_field();
        tMapper.map_input_field_to_output_field_2( tFieldTarget.get() );

        // compute the nodal value
        tFieldTarget->compute_nodal_values();

        // get the coefficents and store them
        mCoefficients = tFieldTarget->get_coefficients();

        // clip the values and
        this->impose_upper_lower_bound( aGENPerformer, tFieldTarget.get() );

        // replace the newly constructed field
        tGENFields( tADVFieldIndex ) = tFieldTarget;

        // store the fields
        mMTKFields = tGENFields;

        // output the fields if asked
        if ( mOutputMeshFile != "" )
        {
            this->output_fields( tFieldTarget.get(), tFieldSource.get(), mOutputMeshFile );
        }
    }

    //------------------------------------------------------------------------------

    moris::sint
    Reinitialize_Performer::get_reinitialization_frequency() const
    {
        return mReinitializationFrequency;
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat > const &
    Reinitialize_Performer::get_coefficients() const
    {
        return mCoefficients;
    }
    //------------------------------------------------------------------------------
    void
    Reinitialize_Performer::impose_upper_lower_bound( Vector< std::shared_ptr< gen::Geometry_Engine > >& aGENPerformer, mtk::Field* aField )
    {
        // lower bound and upper bound are defined on proc 0 and they need to be communicated to other
        // Note:  we make an assumption that all the lower bounds and upper bounds are equal

        // initialize  the upper and lower abound
        moris::real tLowerBound;
        moris::real tUpperBound;

        // assign the values on processor 0
        if ( par_rank() == 0 )
        {
            tLowerBound = aGENPerformer( 0 )->get_lower_bounds()( 0 );
            tUpperBound = aGENPerformer( 0 )->get_upper_bounds()( 0 );
        }    // Bcast the values to other processeors
        MPI_Bcast( &tLowerBound, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
        MPI_Bcast( &tUpperBound, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

        // clip the values of the adv
        for ( uint iADV = 0; iADV < mCoefficients.numel(); iADV++ )
        {
            mCoefficients( iADV ) = std::max( tLowerBound, std::min( mCoefficients( iADV ), tUpperBound ) );
        }

        // update the field based on the newly clipped coeff
        aField->unlock_field();
        aField->set_coefficients( mCoefficients );
        aField->compute_nodal_values();
    }

    //------------------------------------------------------------------------------

    Vector< std::shared_ptr< mtk::Field > >
    Reinitialize_Performer::get_mtk_fields() const
    {
        return mMTKFields;
    }

    //------------------------------------------------------------------------------

    void
    Reinitialize_Performer::output_fields( mtk::Field* aTarget, mtk::Field* aSource, std::string aExoFileName ) const
    {
        //
        Tracer tTracer( "WRK", "Reinitialize ADVs", "Outputting Fields" );
        // time shift
        real tTimeShift = 0.0;

        if ( mTimeOffset > 0 )
        {
            // get optimization iteration
            uint tOptIter = gLogger.get_opt_iteration();

            // set name
            std::string tOptIterStrg = std::to_string( tOptIter );
            aExoFileName += ".e-s." + std::string( 4 - tOptIterStrg.length(), '0' ) + tOptIterStrg;

            // determine time shift
            tTimeShift = tOptIter * mTimeOffset;
        }
        // call the lagrange mesh
        moris::mtk::Mesh* tMesh = aTarget->get_mesh_pair().get_interpolation_mesh();
        // Write mesh
        mtk::Writer_Exodus tWriter( tMesh );
        tWriter.write_mesh( "./", aExoFileName, "./", "gen_temp.exo" );

        // write time to file
        tWriter.set_time( tTimeShift );

        // Set nodal fields based on field names
        Vector< std::string > tNodalFieldNames = { "Mapped_Field", "Original_Field" };
        tWriter.set_nodal_fields( tNodalFieldNames );

        // Create field on mesh
        tWriter.write_nodal_field( tNodalFieldNames( 0 ), aTarget->get_values() );
        tWriter.write_nodal_field( tNodalFieldNames( 1 ), aSource->get_values() );

        // Finalize
        tWriter.close_file( true );
    }

    //------------------------------------------------------------------------------

    void
    Reinitialize_Performer::set_reinit_mode( ReinitMode aMode )
    {
        mReinitMode = aMode;
    }

    //------------------------------------------------------------------------------

    ReinitMode
    Reinitialize_Performer::get_reinit_mode() const
    {
        return mReinitMode;
    }

    //------------------------------------------------------------------------------

    bool
    Reinitialize_Performer::get_initialize_at_start() const
    {
        return mInitializeAtStart;
    }

    //------------------------------------------------------------------------------

    bool
    Reinitialize_Performer::last_reinit_below_wean_tol() const
    {
        return mLastReinitBelowWeanTol;
    }

    //------------------------------------------------------------------------------

    void
    Reinitialize_Performer::set_material_phase( moris_index aPhase )
    {
        mMaterialPhase = aPhase;
    }

    //------------------------------------------------------------------------------

    Vector< moris_index > const&
    Reinitialize_Performer::get_exclude_phases() const
    {
        return mExcludePhases;
    }

    //------------------------------------------------------------------------------

    void
    Reinitialize_Performer::set_exclude_phases( Vector< moris_index > const& aExcludePhases )
    {
        mExcludePhases = aExcludePhases;
    }

    //------------------------------------------------------------------------------

    void
    Reinitialize_Performer::compute_sdf_from_interface(
            xtk::Cut_Integration_Mesh&                         aCutMesh,
            mtk::Mesh*                                         aInterpMesh,
            Vector< std::shared_ptr< gen::Geometry_Engine > >& aGENPerformer,
            Vector< std::shared_ptr< mtk::Mesh_Manager > >&    aMTKPerformer,
            std::shared_ptr< mtk::Field >&                     aSDFField,
            bool                                               aApplySmoothing )
    {
        Tracer tTracer( "WRK", "Reinitialize ADVs", "Compute SDF from Interface" );

        // Get spatial dimension
        uint tSpatialDim = aInterpMesh->get_spatial_dim();
        uint tNodesPerFacet = ( tSpatialDim == 3 ) ? 3 : 2;  // Triangles in 3D, segments in 2D
        
        // Get facet connectivity from XTK's facet-based connectivity
        std::shared_ptr< xtk::Facet_Based_Connectivity > tFaceConn = aCutMesh.get_face_connectivity();

        MORIS_ERROR( tFaceConn != nullptr,
                "compute_sdf_from_interface - Face connectivity not available in XTK mesh" );

        // =========================================================================
        // Phase-based facet filtering using bulk-phase-to-bulk-phase interface
        // =========================================================================
        
        Vector< moris_index > tFilteredFacetIndices;
        
        if ( mExcludePhases.size() > 0 )
        {
            // Get bulk-phase-to-bulk-phase interface
            auto const& tBpToBpInterface = aCutMesh.get_bulk_phase_to_bulk_phase_dbl_side_interface();
            uint tNumPhases = tBpToBpInterface.size();
            
            // Helper to check if a phase is excluded
            auto tIsExcluded = [this]( moris_index aPhase ) {
                for ( uint i = 0; i < mExcludePhases.size(); ++i )
                {
                    if ( mExcludePhases( i ) == aPhase )
                    {
                        return true;
                    }
                }
                return false;
            };
            
            // Iterate through all phase pairs and collect facets from allowed pairs
            uint tTotalFacets = 0;
            uint tExcludedFacets = 0;
            
            for ( uint iP0 = 0; iP0 < tNumPhases; ++iP0 )
            {
                for ( uint iP1 = 0; iP1 < tBpToBpInterface( iP0 ).size(); ++iP1 )
                {
                    // Skip if no interface between these phases
                    if ( tBpToBpInterface( iP0 )( iP1 ) == nullptr )
                    {
                        continue;
                    }
                    
                    // Check if either phase is excluded
                    bool tExcludeThisPair = tIsExcluded( iP0 ) || tIsExcluded( iP1 );
                    
                    // Get the double-side group for this phase pair
                    std::shared_ptr< xtk::IG_Cell_Double_Side_Group > tDblSideGroup = tBpToBpInterface( iP0 )( iP1 );
                    uint tNumFacetsInPair = tDblSideGroup->mLeaderIgCells.size();
                    
                    if ( tExcludeThisPair )
                    {
                        tExcludedFacets += tNumFacetsInPair;
                    }
                    else
                    {
                        // This phase pair is allowed - we don't have direct facet indices here,
                        // so we'll need to use the global interface facets but filter somehow
                        // For now, track that we need these facets
                        tTotalFacets += tNumFacetsInPair;
                    }
                }
            }
            
            MORIS_LOG_INFO( "Phase exclusion: keeping %d facets from %d phases, excluded %d facets from phase pairs containing excluded phases",
                    tTotalFacets, tNumPhases, tExcludedFacets );
            
            // Since the bulk-phase interface doesn't give us direct facet indices,
            // we need to filter the global interface facets based on which phases they connect.
            // We can determine this by checking the bulk phase of cells adjacent to each facet.
            
            Vector< moris_index > const& tAllInterfaceFacetIndices = aCutMesh.get_interface_facets();
            
            // For each facet, check the bulk phases of adjacent cells
            for ( uint iFacet = 0; iFacet < tAllInterfaceFacetIndices.size(); ++iFacet )
            {
                moris_index tFacetIndex = tAllInterfaceFacetIndices( iFacet );
                
                // Skip if facet index is invalid
                if ( tFacetIndex < 0 || (uint)tFacetIndex >= tFaceConn->mFacetToCell.size() )
                {
                    continue;
                }
                
                // Get cells adjacent to this facet
                Vector< moris::mtk::Cell* > const& tAdjacentCells = tFaceConn->mFacetToCell( tFacetIndex );
                
                // Count valid (non-null) adjacent cells and check bulk phases
                uint tValidCellCount = 0;
                bool tFacetHasExcludedPhase = false;
                
                for ( uint iCell = 0; iCell < tAdjacentCells.size(); ++iCell )
                {
                    if ( tAdjacentCells( iCell ) != nullptr )
                    {
                        ++tValidCellCount;
                        moris_index tCellIndex = tAdjacentCells( iCell )->get_index();
                        moris_index tBulkPhase = aCutMesh.get_cell_bulk_phase( tCellIndex );
                        if ( tIsExcluded( tBulkPhase ) )
                        {
                            tFacetHasExcludedPhase = true;
                        }
                    }
                }
                
                // Exclude facet if:
                // 1. It has less than 2 valid adjacent cells (domain boundary facet), OR
                // 2. Any adjacent cell is in an excluded phase
                if ( tValidCellCount < 2 || tFacetHasExcludedPhase )
                {
                    continue;
                }
                
                tFilteredFacetIndices.push_back( tFacetIndex );
            }
            
            MORIS_LOG_INFO( "Phase filtering: %zu facets after filtering (from %zu total)",
                    tFilteredFacetIndices.size(), tAllInterfaceFacetIndices.size() );
        }
        else
        {
            // No phase exclusion - use all interface facets
            Vector< moris_index > const& tAllInterfaceFacetIndices = aCutMesh.get_interface_facets();
            tFilteredFacetIndices = tAllInterfaceFacetIndices;
        }

        Vector< moris_index > const& tInterfaceFacetIndices = tFilteredFacetIndices;

        if ( tInterfaceFacetIndices.size() == 0 )
        {
            MORIS_LOG_WARNING( "No interface facets found after filtering, skipping SDF reinitialization" );
            return;
        }

        // Count facets
        uint tNumFacets = tInterfaceFacetIndices.size();

        // Collect all unique facet node coordinates
        std::map< moris_index, Matrix< DDRMat > > tUniqueVertexCoords;
        Matrix< IndexMat > tFacetConn( tNumFacets, tNodesPerFacet );

        for ( uint iFacet = 0; iFacet < tNumFacets; ++iFacet )
        {
            moris_index tFacetIndex = tInterfaceFacetIndices( iFacet );

            // Bounds check for facet index
            if ( tFacetIndex < 0 || (uint)tFacetIndex >= tFaceConn->mFacetVertices.size() )
            {
                MORIS_LOG_WARNING( "Facet index %d out of bounds (max %zu), skipping",
                        tFacetIndex, tFaceConn->mFacetVertices.size() );
                continue;
            }

            // Get vertices of this facet directly from mFacetVertices member
            Vector< moris::mtk::Vertex* > const& tFacetVertices = tFaceConn->mFacetVertices( tFacetIndex );

            for ( uint iNode = 0; iNode < tNodesPerFacet && iNode < tFacetVertices.size(); ++iNode )
            {
                if ( tFacetVertices( iNode ) == nullptr )
                {
                    continue;  // Skip null vertices
                }

                moris_index tVertexIndex = tFacetVertices( iNode )->get_index();
                tFacetConn( iFacet, iNode ) = tVertexIndex;

                // Store vertex coordinates if not already present
                if ( tUniqueVertexCoords.find( tVertexIndex ) == tUniqueVertexCoords.end() )
                {
                    tUniqueVertexCoords[ tVertexIndex ] = tFacetVertices( iNode )->get_coords();
                }
            }
        }

        // Build facet node coordinate matrix
        uint tNumUniqueVerts = tUniqueVertexCoords.size();
        Matrix< DDRMat > tFacetNodeCoords( tNumUniqueVerts, tSpatialDim );

        // Create mapping from original indices to compact indices
        std::map< moris_index, moris_index > tVertexIndexMap;
        moris_index tCompactIdx = 0;
        for ( auto const& [tOrigIdx, tCoords] : tUniqueVertexCoords )
        {
            tVertexIndexMap[ tOrigIdx ] = tCompactIdx;
            for ( uint iDim = 0; iDim < tSpatialDim; ++iDim )
            {
                tFacetNodeCoords( tCompactIdx, iDim ) = tCoords( iDim );
            }
            ++tCompactIdx;
        }

        // Update facet connectivity to use compact indices
        for ( uint iFacet = 0; iFacet < tNumFacets; ++iFacet )
        {
            for ( uint iNode = 0; iNode < tNodesPerFacet; ++iNode )
            {
                tFacetConn( iFacet, iNode ) = tVertexIndexMap[ tFacetConn( iFacet, iNode ) ];
            }
        }

        // Get interpolation mesh node coordinates
        uint tNumNodes = aInterpMesh->get_num_nodes();
        Matrix< DDRMat > tNodeCoords( tNumNodes, tSpatialDim );

        for ( uint iNode = 0; iNode < tNumNodes; ++iNode )
        {
            Matrix< DDRMat > tCoords = aInterpMesh->get_mtk_vertex( iNode ).get_coords();
            for ( uint iDim = 0; iDim < tSpatialDim; ++iDim )
            {
                tNodeCoords( iNode, iDim ) = tCoords( iDim );
            }
        }

        // Determine bulk phase per node using GEN
        Vector< moris_index > tNodeBulkPhase( tNumNodes );
        for ( uint iNode = 0; iNode < tNumNodes; ++iNode )
        {
            Matrix< DDRMat > tCoords = aInterpMesh->get_mtk_vertex( iNode ).get_coords();
            tNodeBulkPhase( iNode ) = aGENPerformer( 0 )->get_phase_index( iNode, tCoords );
        }

        // Grab the EXISTING level-set field's nodal values so we can preserve its sign.
        // Re-deriving the sign from the bulk phase every reinit is inconsistent with GEN's
        // level-set convention, so the GEN->SDF->GEN round-trip inverts the field at each
        // reinit (period-2 limit cycle). Preserving the existing sign keeps the geometry fixed
        // and makes the reinit a pure redistancing (correct |grad phi|, same inside/outside).
        Matrix< DDRMat > tExistingValues;
        {
            Vector< std::shared_ptr< mtk::Field > > tCurrentFields;
            tCurrentFields.append( aGENPerformer( 0 )->get_mtk_fields() );
            auto tSrcItr = std::find_if( tCurrentFields.begin(), tCurrentFields.end(),
                    [ & ]( std::shared_ptr< mtk::Field > const& aField ) {
                        return aField->get_label() == mADVFiledName;
                    } );
            if ( tSrcItr != tCurrentFields.end() )
            {
                tExistingValues = ( *tSrcItr )->get_values();
            }
        }

        // Compute SDF using SDF_From_Interface
        Matrix< DDRMat > tSDF;
        sdf::SDF_From_Interface::compute(
                tNumNodes,
                tNodeCoords,
                tFacetNodeCoords,
                tFacetConn,
                tNodeBulkPhase,
                mMaterialPhase,
                tSDF );

        // Preserve the sign of the existing field instead of the bulk-phase-derived sign.
        // Falls back to the computed (bulk-phase) sign if the existing field is unavailable.
        if ( tExistingValues.numel() == tSDF.numel() )
        {
            for ( uint iNode = 0; iNode < tSDF.numel(); ++iNode )
            {
                real tSign = ( tExistingValues( iNode ) >= 0.0 ) ? 1.0 : -1.0;
                tSDF( iNode ) = tSign * std::abs( tSDF( iNode ) );
            }
        }

        // Get ADV bounds to clamp the signed distance into the design-variable range
        real tLowerBound = aGENPerformer( 0 )->get_lower_bounds()( 0 );
        real tUpperBound = aGENPerformer( 0 )->get_upper_bounds()( 0 );

        // Helper to check if a phase is excluded (for logging only; the facets of excluded
        // phases were already filtered out of the interface above, so excluded-phase nodes
        // already carry the true distance to the nearest *kept* facet)
        uint tPhaseExcludedNodes = 0;
        auto tIsExcludedPhase = [this]( moris_index aPhase ) {
            for ( uint i = 0; i < mExcludePhases.size(); ++i )
            {
                if ( mExcludePhases( i ) == aPhase )
                {
                    return true;
                }
            }
            return false;
        };

        // Clamp the true signed distance to the ADV bounds. Excluded-phase and domain-boundary
        // nodes are NOT saturated to +/- the bound magnitude: they keep the genuine distance to
        // the nearest (non-excluded) interface facet, so the field stays a true signed distance
        // everywhere instead of forming an artificial +/-bound wall at the frame/boundary.
        for ( uint iNode = 0; iNode < tNumNodes; ++iNode )
        {
            if ( mExcludePhases.size() > 0 && tIsExcludedPhase( tNodeBulkPhase( iNode ) ) )
            {
                ++tPhaseExcludedNodes;
            }

            tSDF( iNode ) = std::max( tLowerBound, std::min( tSDF( iNode ), tUpperBound ) );
        }

        MORIS_LOG_INFO( "SDF clamped to [%.4f, %.4f], phase-excluded nodes: %d (true distance kept, no boundary saturation)",
                tLowerBound, tUpperBound, tPhaseExcludedNodes );

        // Diagnostic only: the relative correction this reinit makes to the field, ||SDF(x)-x||/||x||.
        // This is the field's between-reinit drift (dominated by |grad phi| re-normalization) and does
        // NOT shrink at convergence, so it is logged for insight but no longer drives the wean-off
        // decision. Weaning is decided in the workflow from the actual design-variable change.
        mLastReinitBelowWeanTol = false;
        if ( tExistingValues.numel() == tSDF.numel() )
        {
            real tNumer = 0.0;
            real tDenom = 0.0;
            for ( uint iNode = 0; iNode < tSDF.numel(); ++iNode )
            {
                tNumer += std::pow( tSDF( iNode ) - tExistingValues( iNode ), 2 );
                tDenom += std::pow( tExistingValues( iNode ), 2 );
            }
            real tRelChange = ( tDenom > 0.0 ) ? std::sqrt( tNumer / tDenom ) : 0.0;
            MORIS_LOG_INFO( "SDF reinit field drift (diagnostic): %.4e", tRelChange );
        }

        // Create the MTK Field holding the computed SDF values. Build it on the ADV field's
        // OWN mesh pair + discretization_mesh_index so the resulting coefficients align with
        // the GEN ADV vector. (Previously the INTERFACE_SDF path used mesh_pair(0)/index 0,
        // whose coefficient count never matched the ADVs, so every in-place apply was skipped
        // with "coeff count != ADV count".) This mirrors the field_remap path above; falls
        // back to mesh_pair(0)/index 0 if the ADV field cannot be located by label.
        Vector< std::shared_ptr< mtk::Field > > tAdvLookupFields;
        tAdvLookupFields.append( aGENPerformer( 0 )->get_mtk_fields() );
        auto tAdvItr = std::find_if( tAdvLookupFields.begin(), tAdvLookupFields.end(),
                [ this ]( std::shared_ptr< mtk::Field > const& aField ) {
                    return aField->get_label() == mADVFiledName;
                } );

        if ( tAdvItr != tAdvLookupFields.end() )
        {
            aSDFField = std::make_shared< mtk::Field_Discrete >(
                    ( *tAdvItr )->get_mesh_pair(), ( *tAdvItr )->get_discretization_mesh_index() );
        }
        else
        {
            moris_index tMeshIndex = ( mReinitMode == ReinitMode::INTERFACE_SDF ) ? 0 : mAdofMeshIndex;
            aSDFField = std::make_shared< mtk::Field_Discrete >(
                    aMTKPerformer( 0 )->get_mesh_pair( 0 ), tMeshIndex );
        }
        aSDFField->set_label( mADVFiledName );
        aSDFField->unlock_field();
        aSDFField->set_values( tSDF );

        // Apply Laplacian smoothing if enabled (and requested this call). In every-step mode
        // smoothing is decoupled from redistancing and only requested every reinit_freq iterations.
        if ( mSmoothingIterations > 0 && aApplySmoothing )
        {
            apply_laplacian_smoothing( tSDF, aInterpMesh, mSmoothingIterations, mSmoothingLambda );
            aSDFField->unlock_field();
            aSDFField->set_values( tSDF );
        }

        // Map nodal values to coefficients
        mtk::Mapper tMapper;
        aSDFField->unlock_field();
        tMapper.map_input_field_to_output_field_2( aSDFField.get() );

        // Compute coefficients from nodal values
        mCoefficients = aSDFField->get_coefficients();

        // Clip values to bounds
        this->impose_upper_lower_bound( aGENPerformer, aSDFField.get() );

        // Store field for GEN update
        Vector< std::shared_ptr< mtk::Field > > tGENFields;
        tGENFields.append( aGENPerformer( 0 )->get_mtk_fields() );

        // Find and replace the ADV field
        auto itr = std::find_if( tGENFields.begin(), tGENFields.end(),
                [ & ]( std::shared_ptr< mtk::Field > const& aField ) {
                    return aField->get_label() == mADVFiledName;
                } );

        if ( itr != tGENFields.end() )
        {
            moris_index tADVFieldIndex = std::distance( tGENFields.begin(), itr );
            tGENFields( tADVFieldIndex ) = aSDFField;
        }

        mMTKFields = tGENFields;

        MORIS_LOG_INFO( "SDF reinitialization complete: %d nodes, %d interface facets",
                tNumNodes, tNumFacets );
    }

    //------------------------------------------------------------------------------

    void
    Reinitialize_Performer::apply_laplacian_smoothing(
            Matrix< DDRMat >&   aCoefficients,
            mtk::Mesh*          aInterpMesh,
            uint                aNumIterations,
            real                aLambda )
    {
        // Get number of nodes
        uint tNumNodes = aCoefficients.numel();

        MORIS_LOG_INFO( "Applying Laplacian smoothing: %d nodes, %d iterations, lambda=%.2f",
                tNumNodes, aNumIterations, aLambda );

        // For each smoothing iteration
        for ( uint iIter = 0; iIter < aNumIterations; ++iIter )
        {
            // Create copy of current values
            Matrix< DDRMat > tOldValues = aCoefficients;

            // For each node, compute weighted average with neighbors
            for ( uint iNode = 0; iNode < tNumNodes; ++iNode )
            {
                // Get elements connected to this node
                Matrix< IndexMat > tConnectedElements = aInterpMesh->get_elements_connected_to_node_loc_inds( iNode );

                // Collect unique neighbor nodes
                std::set< moris_index > tNeighborNodes;

                for ( uint iElem = 0; iElem < tConnectedElements.numel(); ++iElem )
                {
                    moris_index tElemIndex = tConnectedElements( iElem );

                    // Get nodes of this element
                    Matrix< IndexMat > tElemNodes = aInterpMesh->get_nodes_connected_to_element_loc_inds( tElemIndex );

                    for ( uint iElemNode = 0; iElemNode < tElemNodes.numel(); ++iElemNode )
                    {
                        moris_index tNeighborNode = tElemNodes( iElemNode );
                        if ( tNeighborNode != (moris_index)iNode && tNeighborNode < (moris_index)tNumNodes )
                        {
                            tNeighborNodes.insert( tNeighborNode );
                        }
                    }
                }

                // Compute average of neighbor values
                if ( !tNeighborNodes.empty() )
                {
                    real tNeighborSum = 0.0;
                    for ( moris_index tNeighbor : tNeighborNodes )
                    {
                        tNeighborSum += tOldValues( tNeighbor );
                    }
                    real tNeighborAvg = tNeighborSum / tNeighborNodes.size();

                    // Apply Laplacian smoothing: new = (1-lambda)*old + lambda*neighbor_avg
                    aCoefficients( iNode ) = ( 1.0 - aLambda ) * tOldValues( iNode ) + aLambda * tNeighborAvg;
                }
                // If no neighbors (shouldn't happen), keep original value
            }
        }

        MORIS_LOG_INFO( "Laplacian smoothing complete. ADV range: [%.6f, %.6f]",
                aCoefficients.min(), aCoefficients.max() );
    }

}    // namespace moris::wrk
