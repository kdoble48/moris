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
#include <unordered_map>
#include <vector>
#include "cl_Matrix.hpp"
#include "cl_Tracer.hpp"
#include "cl_Logger.hpp"
#include "cl_MTK_Mesh_Pair.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_Communication_Tools.hpp"

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

    namespace
    {
        // The helpers below make the interface redistance collective. On one processor each of them returns its input.

        /** concatenation of every processor's vector (in processor order) on processor 0; empty on the others */
        template< typename T >
        std::vector< T >
        gather_on_proc_0( const std::vector< T >& aLocal, MPI_Datatype aType )
        {
            if ( par_size() == 1 )
            {
                return aLocal;
            }
            int                tNumLocal = (int)aLocal.size();
            std::vector< int > tCounts( par_size(), 0 );
            MPI_Gather( &tNumLocal, 1, MPI_INT, tCounts.data(), 1, MPI_INT, 0, get_comm() );

            std::vector< int > tOffsets( par_size(), 0 );
            std::vector< T >   tAll;
            if ( par_rank() == 0 )
            {
                int tTotal = 0;
                for ( int iProc = 0; iProc < (int)par_size(); iProc++ )
                {
                    tOffsets[ iProc ] = tTotal;
                    tTotal += tCounts[ iProc ];
                }
                tAll.resize( tTotal );
            }
            MPI_Gatherv( aLocal.data(), tNumLocal, aType, tAll.data(), tCounts.data(), tOffsets.data(), aType, 0, get_comm() );
            return tAll;
        }

        /** concatenation of every processor's vector (in processor order), on every processor */
        std::vector< real >
        gather_on_all_procs( const std::vector< real >& aLocal )
        {
            if ( par_size() == 1 )
            {
                return aLocal;
            }
            int                tNumLocal = (int)aLocal.size();
            std::vector< int > tCounts( par_size(), 0 );
            MPI_Allgather( &tNumLocal, 1, MPI_INT, tCounts.data(), 1, MPI_INT, get_comm() );

            std::vector< int > tOffsets( par_size(), 0 );
            int                tTotal = 0;
            for ( int iProc = 0; iProc < (int)par_size(); iProc++ )
            {
                tOffsets[ iProc ] = tTotal;
                tTotal += tCounts[ iProc ];
            }
            std::vector< real > tAll( tTotal );
            MPI_Allgatherv( aLocal.data(), tNumLocal, MPI_DOUBLE, tAll.data(), tCounts.data(), tOffsets.data(), MPI_DOUBLE, get_comm() );
            return tAll;
        }

        /** processor 0's vector on every processor */
        template< typename T >
        void
        broadcast_from_proc_0( std::vector< T >& aVector, MPI_Datatype aType )
        {
            if ( par_size() == 1 )
            {
                return;
            }
            int tSize = (int)aVector.size();
            MPI_Bcast( &tSize, 1, MPI_INT, 0, get_comm() );
            aVector.resize( tSize );
            MPI_Bcast( aVector.data(), tSize, aType, 0, get_comm() );
        }

        /**
         * The interface facets of all processors as one facet soup (every facet with its own nodes). A processor holds
         * the facets of its own and its aura cells only, but the distance is to the whole interface. A facet that
         * several processors hold is listed several times, which does not change the minimum distance.
         */
        void
        gather_interface_facets(
                Matrix< DDRMat >&   aFacetNodeCoords,
                Matrix< IndexMat >& aFacetConn,
                uint                aSpatialDim )
        {
            uint tNumLocalFacets = aFacetConn.n_rows();
            uint tNodesPerFacet  = aFacetConn.n_cols();
            uint tStride         = tNodesPerFacet * aSpatialDim;

            std::vector< real > tLocal( tNumLocalFacets * tStride, 0.0 );
            for ( uint iFacet = 0; iFacet < tNumLocalFacets; ++iFacet )
            {
                for ( uint iNode = 0; iNode < tNodesPerFacet; ++iNode )
                {
                    for ( uint iDim = 0; iDim < aSpatialDim; ++iDim )
                    {
                        tLocal[ iFacet * tStride + iNode * aSpatialDim + iDim ] = aFacetNodeCoords( aFacetConn( iFacet, iNode ), iDim );
                    }
                }
            }

            std::vector< real > tAll       = gather_on_all_procs( tLocal );
            uint                tNumFacets = tAll.size() / tStride;

            aFacetNodeCoords.set_size( tNumFacets * tNodesPerFacet, aSpatialDim );
            aFacetConn.set_size( tNumFacets, tNodesPerFacet );
            for ( uint iFacet = 0; iFacet < tNumFacets; ++iFacet )
            {
                for ( uint iNode = 0; iNode < tNodesPerFacet; ++iNode )
                {
                    aFacetConn( iFacet, iNode ) = (moris_index)( iFacet * tNodesPerFacet + iNode );
                    for ( uint iDim = 0; iDim < aSpatialDim; ++iDim )
                    {
                        aFacetNodeCoords( iFacet * tNodesPerFacet + iNode, iDim ) = tAll[ iFacet * tStride + iNode * aSpatialDim + iDim ];
                    }
                }
            }
        }

        /** nodal values of the owning processor on every processor that holds the node (shared and aura nodes) */
        void
        synchronize_nodal_values( Matrix< DDRMat >& aNodalValues, mtk::Mesh* aMesh )
        {
            if ( par_size() == 1 )
            {
                return;
            }
            uint     tNumNodes = aNodalValues.numel();
            moris_id tMaxId    = 0;
            for ( uint iNode = 0; iNode < tNumNodes; ++iNode )
            {
                tMaxId = std::max( tMaxId, aMesh->get_glb_entity_id_from_entity_loc_index( (moris_index)iNode, mtk::EntityRank::NODE ) );
            }
            tMaxId = max_all( tMaxId );

            std::vector< real > tOwned( tMaxId + 1, 0.0 );
            std::vector< real > tGlobal( tMaxId + 1, 0.0 );
            for ( uint iNode = 0; iNode < tNumNodes; ++iNode )
            {
                if ( aMesh->get_entity_owner( (moris_index)iNode, mtk::EntityRank::NODE ) == (uint)par_rank() )
                {
                    tOwned[ aMesh->get_glb_entity_id_from_entity_loc_index( (moris_index)iNode, mtk::EntityRank::NODE ) ] = aNodalValues( iNode );
                }
            }
            MPI_Allreduce( tOwned.data(), tGlobal.data(), (int)tOwned.size(), MPI_DOUBLE, MPI_SUM, get_comm() );
            for ( uint iNode = 0; iNode < tNumNodes; ++iNode )
            {
                aNodalValues( iNode ) = tGlobal[ aMesh->get_glb_entity_id_from_entity_loc_index( (moris_index)iNode, mtk::EntityRank::NODE ) ];
            }
        }

        /** owning processor of every (local) design coefficient; -1 for a coefficient that no node uses */
        Vector< moris_id >
        coefficient_owners( mtk::Interpolation_Mesh* aMesh, uint aDiscretizationIndex, uint aNumCoefficients )
        {
            Vector< moris_id > tOwners( aNumCoefficients, -1 );
            uint               tNumNodes = aMesh->get_num_nodes();
            for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
            {
                if ( not aMesh->get_mtk_vertex( iNode ).has_interpolation( aDiscretizationIndex ) )
                {
                    continue;
                }
                Matrix< IndexMat > tIndices   = aMesh->get_coefficient_indices_of_node( iNode, aDiscretizationIndex );
                Matrix< IdMat >    tNodeOwner = aMesh->get_coefficient_owners_of_node( iNode, aDiscretizationIndex );
                for ( uint iC = 0; iC < tIndices.numel(); iC++ )
                {
                    if ( tIndices( iC ) >= 0 && (uint)tIndices( iC ) < aNumCoefficients )
                    {
                        tOwners( tIndices( iC ) ) = tNodeOwner( iC );
                    }
                }
            }
            return tOwners;
        }
    }    // namespace

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

        // in-place taper: patience (consecutive accepted iterates below the tolerance) and hard cutoff
        mWeanPatience = (uint)std::max( 1, (int)tOPTParameterList( 0 )( 0 ).get< sint >( "sdf_reinit_wean_patience" ) );
        mStopIter     = (uint)std::max( 0, (int)tOPTParameterList( 0 )( 0 ).get< sint >( "sdf_reinit_stop_iter" ) );

        // Redistance every iteration in-place (no GCMMA restart); smooth every reinit_freq
        mReinitEveryStep = tOPTParameterList( 0 )( 0 ).get< bool >( "sdf_reinit_every_step" );

        // footprint diagnostics of the in-place redistance
        mDiagnostics = tOPTParameterList( 0 )( 0 ).get< bool >( "sdf_reinit_diagnostics" );

        // interface-preserving in-place redistance (freeze every coefficient that can reach a cut cell)
        mPreserveInterface = tOPTParameterList( 0 )( 0 ).get< bool >( "sdf_preserve_interface" );

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

        MORIS_ERROR( itr != tGENFields.end(),
                "Reinitialize_Performer: no GEN field is labelled '%s' (the mapping sublist's adv_field). Check the design "
                "geometry's \"name\" parameter.", mADVFiledName.c_str() );

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

        // the interface is the union of every processor's facets: all processors take this branch together
        if ( sum_all( (uint)tInterfaceFacetIndices.size() ) == 0 )
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

        // on several processors: the distance is to the whole interface, not to this processor's part of it (MB-113)
        if ( par_size() > 1 )
        {
            gather_interface_facets( tFacetNodeCoords, tFacetConn, tSpatialDim );
            tNumFacets = tFacetConn.n_rows();
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
        // (the full bound vectors live on processor 0; as in impose_upper_lower_bound() all bounds are taken as equal)
        real tLowerBound = 0.0;
        real tUpperBound = 0.0;
        if ( par_rank() == 0 )
        {
            tLowerBound = aGENPerformer( 0 )->get_lower_bounds()( 0 );
            tUpperBound = aGENPerformer( 0 )->get_upper_bounds()( 0 );
        }
        if ( par_size() > 1 )
        {
            MPI_Bcast( &tLowerBound, 1, MPI_DOUBLE, 0, get_comm() );
            MPI_Bcast( &tUpperBound, 1, MPI_DOUBLE, 0, get_comm() );
        }

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
                // every node once: on its owning processor
                if ( par_size() > 1 && aInterpMesh->get_entity_owner( (moris_index)iNode, mtk::EntityRank::NODE ) != (uint)par_rank() )
                {
                    continue;
                }
                tNumer += std::pow( tSDF( iNode ) - tExistingValues( iNode ), 2 );
                tDenom += std::pow( tExistingValues( iNode ), 2 );
            }
            tNumer = sum_all( tNumer );
            tDenom = sum_all( tDenom );
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

        // a missing label used to fall back silently to mesh pair 0 / index 0 (sdf_reinit_sweep layer 3, MB-019)
        MORIS_ERROR( tAdvItr != tAdvLookupFields.end(),
                "Reinitialize_Performer (interface SDF): no GEN field is labelled '%s' (the mapping sublist's adv_field). "
                "Check the design geometry's \"name\" parameter.", mADVFiledName.c_str() );
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

    namespace
    {
        /** the ADV field (by label) among GEN's MTK fields; MORIS_ERROR if absent */
        std::shared_ptr< mtk::Field >
        find_adv_field( Vector< std::shared_ptr< gen::Geometry_Engine > >& aGENPerformer, const std::string& aLabel )
        {
            Vector< std::shared_ptr< mtk::Field > > tFields;
            tFields.append( aGENPerformer( 0 )->get_mtk_fields() );
            for ( uint iField = 0; iField < tFields.size(); iField++ )
            {
                if ( tFields( iField )->get_label() == aLabel )
                {
                    return tFields( iField );
                }
            }
            MORIS_ERROR( false, "Reinitialize_Performer: no GEN field is labelled '%s'.", aLabel.c_str() );
            return nullptr;
        }
    }    // namespace

    //------------------------------------------------------------------------------

    Vector< Vector< moris_index > >
    Reinitialize_Performer::cell_coefficients( mtk::Interpolation_Mesh* aMesh, uint aDiscretizationIndex )
    {
        uint                            tNumCells = aMesh->get_num_elems();
        Vector< Vector< moris_index > > tCellCoefficients( tNumCells );
        for ( uint iCell = 0; iCell < tNumCells; iCell++ )
        {
            Matrix< IndexMat >      tNodes = aMesh->get_nodes_connected_to_element_loc_inds( (moris_index)iCell );
            std::set< moris_index > tUnion;
            for ( uint iNode = 0; iNode < tNodes.numel(); iNode++ )
            {
                // on several processors the outer nodes of the aura carry no interpolation; the cells of this processor do
                if ( not aMesh->get_mtk_vertex( tNodes( iNode ) ).has_interpolation( aDiscretizationIndex ) )
                {
                    continue;
                }
                Matrix< IndexMat > tCoefficients = aMesh->get_coefficient_indices_of_node( (uint)tNodes( iNode ), aDiscretizationIndex );
                for ( uint iC = 0; iC < tCoefficients.numel(); iC++ )
                {
                    tUnion.insert( tCoefficients( iC ) );
                }
            }
            tCellCoefficients( iCell ).reserve( tUnion.size() );
            for ( moris_index tIndex : tUnion )
            {
                tCellCoefficients( iCell ).push_back( tIndex );
            }
        }
        return tCellCoefficients;
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat >
    Reinitialize_Performer::nodal_values( mtk::Interpolation_Mesh* aMesh, uint aDiscretizationIndex, const Matrix< DDRMat >& aCoeffs )
    {
        uint             tNumNodes = aMesh->get_num_nodes();
        Matrix< DDRMat > tValues( tNumNodes, 1, 0.0 );
        for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
        {
            // on several processors the outer nodes of the aura carry no interpolation (their value stays 0)
            if ( not aMesh->get_mtk_vertex( iNode ).has_interpolation( aDiscretizationIndex ) )
            {
                continue;
            }
            const Matrix< DDRMat >& tT            = aMesh->get_t_matrix_of_node_loc_ind( iNode, aDiscretizationIndex );
            Matrix< IndexMat >      tCoefficients = aMesh->get_coefficient_indices_of_node( iNode, aDiscretizationIndex );
            real                    tValue        = 0.0;
            for ( uint iC = 0; iC < tCoefficients.numel(); iC++ )
            {
                tValue += tT( iC ) * aCoeffs( tCoefficients( iC ) );
            }
            tValues( iNode ) = tValue;
        }
        return tValues;
    }

    //------------------------------------------------------------------------------

    uint
    Reinitialize_Performer::freeze_interface_band(
            Matrix< DDRMat >&                      aNewCoeffs,
            const Matrix< DDRMat >&                aOldCoeffs,
            const Vector< Vector< moris_index > >& aCellCoefficients,
            real                                   aTolerance,
            uint&                                  aRings )
    {
        MORIS_ERROR( aNewCoeffs.numel() == aOldCoeffs.numel(),
                "Reinitialize_Performer::freeze_interface_band: %zu new vs %zu old coefficients.",
                (size_t)aNewCoeffs.numel(), (size_t)aOldCoeffs.numel() );

        uint                tNumCells = aCellCoefficients.size();
        std::vector< bool > tFrozen( aOldCoeffs.numel(), false );
        std::vector< int >  tCellSign( tNumCells, 0 );    // +1 / -1: psi strictly signed on the cell, 0: possibly cut

        // 1. possibly cut cells (psi coefficients not all of one strict sign): freeze all their coefficients
        for ( uint iCell = 0; iCell < tNumCells; iCell++ )
        {
            const Vector< moris_index >& tCoeffs = aCellCoefficients( iCell );
            bool                         tAllPos = tCoeffs.size() > 0;
            bool                         tAllNeg = tCoeffs.size() > 0;
            for ( uint iC = 0; iC < tCoeffs.size(); iC++ )
            {
                real tValue = aOldCoeffs( tCoeffs( iC ) );
                tAllPos     = tAllPos && tValue > aTolerance;
                tAllNeg     = tAllNeg && tValue < -aTolerance;
            }
            tCellSign[ iCell ] = tAllPos ? 1 : ( tAllNeg ? -1 : 0 );
            if ( tCellSign[ iCell ] == 0 )
            {
                for ( uint iC = 0; iC < tCoeffs.size(); iC++ )
                {
                    tFrozen[ tCoeffs( iC ) ] = true;
                }
            }
        }

        // 2. strictly signed cells must stay strictly signed with the same sign; grow the frozen set until they do
        aRings = 0;
        bool tChanged = true;
        while ( tChanged )
        {
            tChanged = false;
            for ( uint iC = 0; iC < aOldCoeffs.numel(); iC++ )
            {
                if ( tFrozen[ iC ] )
                {
                    aNewCoeffs( iC ) = aOldCoeffs( iC );
                }
            }
            for ( uint iCell = 0; iCell < tNumCells; iCell++ )
            {
                if ( tCellSign[ iCell ] == 0 )
                {
                    continue;
                }
                const Vector< moris_index >& tCoeffs = aCellCoefficients( iCell );
                bool                         tOk     = true;
                for ( uint iC = 0; iC < tCoeffs.size() && tOk; iC++ )
                {
                    real tValue = aNewCoeffs( tCoeffs( iC ) );
                    tOk         = tCellSign[ iCell ] > 0 ? tValue > aTolerance : tValue < -aTolerance;
                }
                if ( !tOk )
                {
                    for ( uint iC = 0; iC < tCoeffs.size(); iC++ )
                    {
                        if ( !tFrozen[ tCoeffs( iC ) ] )
                        {
                            tFrozen[ tCoeffs( iC ) ] = true;
                            tChanged                 = true;
                        }
                    }
                }
            }
            aRings += tChanged ? 1 : 0;
        }

        uint tNumFrozen = 0;
        for ( uint iC = 0; iC < aOldCoeffs.numel(); iC++ )
        {
            tNumFrozen += tFrozen[ iC ] ? 1 : 0;
        }
        return tNumFrozen;
    }

    //------------------------------------------------------------------------------

    bool
    Reinitialize_Performer::apply_redistance_to_advs(
            Vector< std::shared_ptr< gen::Geometry_Engine > >& aGENPerformer,
            Vector< real >&                                    aADVs,
            uint&                                              aNumWritten,
            uint&                                              aNumCoefficients ) const
    {
        aNumWritten      = 0;
        aNumCoefficients = 0;

        std::shared_ptr< mtk::Field > tField     = find_adv_field( aGENPerformer, mADVFiledName );
        mtk::Mesh_Pair                tMeshPair  = tField->get_mesh_pair();
        mtk::Interpolation_Mesh*      tMesh      = tMeshPair.get_interpolation_mesh();
        uint                          tMeshIndex = (uint)tField->get_discretization_mesh_index();

        // global ADV ID of every coefficient this processor holds (owned, shared and aura)
        Vector< sint > tLocalADVIds = aGENPerformer( 0 )->get_design_determining_adv_ids( mADVFiledName );

        MORIS_ERROR( tLocalADVIds.size() > 0,
                "SDF in-place redistance: GEN has no design named '%s' (the mapping sublist's adv_field), so the "
                "redistance cannot be applied (MB-019: this used to be skipped silently every iteration).",
                mADVFiledName.c_str() );

        // the redistanced coefficients must be those of the design on every processor
        uint tNumLocal = mCoefficients.n_rows();
        uint tMatch    = ( tLocalADVIds.size() == tNumLocal && tNumLocal > 0 ) ? 1 : 0;
        if ( min_all( tMatch ) == 0 )
        {
            aNumCoefficients = tNumLocal;
            return false;
        }

        // 1. every coefficient once, from its owning processor, onto processor 0: ADV ID and redistanced value
        Vector< moris_id >  tOwners = coefficient_owners( tMesh, tMeshIndex, tNumLocal );
        std::vector< sint > tOwnedIds;
        std::vector< real > tOwnedValues;
        tOwnedIds.reserve( tNumLocal );
        tOwnedValues.reserve( tNumLocal );
        for ( uint iC = 0; iC < tNumLocal; iC++ )
        {
            if ( tOwners( iC ) == par_rank() )
            {
                tOwnedIds.push_back( tLocalADVIds( iC ) );
                tOwnedValues.push_back( mCoefficients( iC ) );
            }
        }
        std::vector< sint > tIds    = gather_on_proc_0( tOwnedIds, MPI_INT );
        std::vector< real > tValues = gather_on_proc_0( tOwnedValues, MPI_DOUBLE );

        // 2. for the freeze band: the coefficients of every cell, from the processor that owns the cell (as ADV IDs)
        std::vector< sint > tCellSizes;
        std::vector< sint > tCellIds;
        if ( mPreserveInterface )
        {
            Vector< Vector< moris_index > > tLocalCells = cell_coefficients( tMesh, tMeshIndex );
            for ( uint iCell = 0; iCell < tLocalCells.size(); iCell++ )
            {
                if ( par_size() > 1 && tMesh->get_mtk_cell( (moris_index)iCell ).get_owner() != par_rank() )
                {
                    continue;
                }
                tCellSizes.push_back( (sint)tLocalCells( iCell ).size() );
                for ( uint iC = 0; iC < tLocalCells( iCell ).size(); iC++ )
                {
                    tCellIds.push_back( tLocalADVIds( tLocalCells( iCell )( iC ) ) );
                }
            }
            tCellSizes = gather_on_proc_0( tCellSizes, MPI_INT );
            tCellIds   = gather_on_proc_0( tCellIds, MPI_INT );
        }

        // 3. processor 0 holds the optimizer's ADV vector: psi from the design's slice, freeze, write phi back
        std::vector< real > tOld( tIds.size(), 0.0 );
        if ( par_rank() == 0 )
        {
            uint tNumCoefficients = tIds.size();

            Vector< sint > tIdVector( tNumCoefficients );
            for ( uint iC = 0; iC < tNumCoefficients; iC++ )
            {
                tIdVector( iC ) = tIds[ iC ];
            }
            // positions in aADVs, not global ADV IDs (MB-017: using the IDs shifted the design by one index)
            Vector< sint > tPositions = aGENPerformer( 0 )->get_adv_positions( tIdVector );

            // psi: the design coefficients this evaluation was called with
            Matrix< DDRMat > tOldCoeffs( tNumCoefficients, 1, 0.0 );
            Matrix< DDRMat > tNewCoeffs( tNumCoefficients, 1, 0.0 );
            for ( uint iC = 0; iC < tNumCoefficients; iC++ )
            {
                sint tPosition = tPositions( iC );
                if ( tPosition >= 0 && tPosition < (sint)aADVs.size() )
                {
                    tOldCoeffs( iC ) = aADVs( tPosition );
                }
                tNewCoeffs( iC ) = tValues[ iC ];
            }

            // sdf_preserve_interface (MB-027): keep every coefficient that can reach a cut cell, so the cut, the criteria
            // and their sensitivities are those of psi bit for bit
            if ( mPreserveInterface )
            {
                std::unordered_map< sint, moris_index > tIndexOfId;
                tIndexOfId.reserve( tNumCoefficients );
                for ( uint iC = 0; iC < tNumCoefficients; iC++ )
                {
                    tIndexOfId[ tIds[ iC ] ] = (moris_index)iC;
                }

                Vector< Vector< moris_index > > tCells( tCellSizes.size() );
                uint                            tCursor = 0;
                for ( uint iCell = 0; iCell < tCellSizes.size(); iCell++ )
                {
                    tCells( iCell ).reserve( tCellSizes[ iCell ] );
                    for ( sint iC = 0; iC < tCellSizes[ iCell ]; iC++ )
                    {
                        std::unordered_map< sint, moris_index >::const_iterator tFound = tIndexOfId.find( tCellIds[ tCursor++ ] );
                        MORIS_ERROR( tFound != tIndexOfId.end(),
                                "sdf_preserve_interface: a cell uses a design coefficient that no processor owns." );
                        tCells( iCell ).push_back( tFound->second );
                    }
                }

                // strict-sign tolerance: GEN's default isocontour_tolerance, so no free coefficient lands where GEN
                // would treat a node as on the interface
                uint tRings  = 0;
                uint tFrozen = freeze_interface_band( tNewCoeffs, tOldCoeffs, tCells, 1.0e-12, tRings );

                MORIS_LOG_INFO( "SDF redistance preserves the interface: %u of %u coefficients frozen (%u growth passes)",
                        tFrozen, tNumCoefficients, tRings );
            }

            for ( uint iC = 0; iC < tNumCoefficients; iC++ )
            {
                sint tPosition = tPositions( iC );
                if ( tPosition >= 0 && tPosition < (sint)aADVs.size() )
                {
                    aADVs( tPosition ) = tNewCoeffs( iC );
                    ++aNumWritten;
                }
                tOld[ iC ]    = tOldCoeffs( iC );
                tValues[ iC ] = tNewCoeffs( iC );
            }
            aNumCoefficients = tNumCoefficients;
        }

        // every processor continues with the same decision (re-cut or not)
        broadcast( aNumWritten );
        broadcast( aNumCoefficients );

        // 4. sdf_reinit_diagnostics: psi and phi on this processor's coefficients
        if ( mDiagnostics )
        {
            broadcast_from_proc_0( tIds, MPI_INT );
            broadcast_from_proc_0( tOld, MPI_DOUBLE );
            broadcast_from_proc_0( tValues, MPI_DOUBLE );

            std::unordered_map< sint, uint > tIndexOfId;
            tIndexOfId.reserve( tIds.size() );
            for ( uint iC = 0; iC < tIds.size(); iC++ )
            {
                tIndexOfId[ tIds[ iC ] ] = iC;
            }
            Matrix< DDRMat > tOldLocal( tNumLocal, 1, 0.0 );
            Matrix< DDRMat > tNewLocal( tNumLocal, 1, 0.0 );
            for ( uint iC = 0; iC < tNumLocal; iC++ )
            {
                std::unordered_map< sint, uint >::const_iterator tFound = tIndexOfId.find( tLocalADVIds( iC ) );
                if ( tFound != tIndexOfId.end() )
                {
                    tOldLocal( iC ) = tOld[ tFound->second ];
                    tNewLocal( iC ) = tValues[ tFound->second ];
                }
            }
            this->log_redistance_footprint( aGENPerformer, tOldLocal, tNewLocal );
        }

        return true;
    }

    //------------------------------------------------------------------------------

    Vector< Vector< moris_index > >
    Reinitialize_Performer::get_cell_coefficients( Vector< std::shared_ptr< gen::Geometry_Engine > >& aGENPerformer ) const
    {
        std::shared_ptr< mtk::Field > tField    = find_adv_field( aGENPerformer, mADVFiledName );
        mtk::Mesh_Pair                tMeshPair = tField->get_mesh_pair();
        return cell_coefficients( tMeshPair.get_interpolation_mesh(), (uint)tField->get_discretization_mesh_index() );
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat >
    Reinitialize_Performer::get_nodal_values(
            Vector< std::shared_ptr< gen::Geometry_Engine > >& aGENPerformer,
            const Matrix< DDRMat >&                            aCoeffs ) const
    {
        std::shared_ptr< mtk::Field > tField    = find_adv_field( aGENPerformer, mADVFiledName );
        mtk::Mesh_Pair                tMeshPair = tField->get_mesh_pair();
        return nodal_values( tMeshPair.get_interpolation_mesh(), (uint)tField->get_discretization_mesh_index(), aCoeffs );
    }

    //------------------------------------------------------------------------------

    void
    Reinitialize_Performer::log_redistance_footprint(
            Vector< std::shared_ptr< gen::Geometry_Engine > >& aGENPerformer,
            const Matrix< DDRMat >&                            aOldCoeffs,
            const Matrix< DDRMat >&                            aNewCoeffs ) const
    {
        std::shared_ptr< mtk::Field > tField     = find_adv_field( aGENPerformer, mADVFiledName );
        mtk::Mesh_Pair                tMeshPair  = tField->get_mesh_pair();
        mtk::Interpolation_Mesh*      tMesh      = tMeshPair.get_interpolation_mesh();
        uint                          tMeshIndex = (uint)tField->get_discretization_mesh_index();

        Vector< Vector< moris_index > > tCellCoefficients = this->get_cell_coefficients( aGENPerformer );
        Matrix< DDRMat >                tOldNodal         = this->get_nodal_values( aGENPerformer, aOldCoeffs );
        Matrix< DDRMat >                tNewNodal         = this->get_nodal_values( aGENPerformer, aNewCoeffs );

        uint tNumCells = tCellCoefficients.size();
        uint tNumNodes = tOldNodal.numel();
        uint tNumCoeff = aOldCoeffs.numel();

        // cells, nodes and coefficients are counted once, on their owning processor; the band (everything that touches
        // a cell cut in psi) is marked by the owner of the cell and merged over the processors below
        bool               tParallel = par_size() > 1;
        Vector< moris_id > tOwners   = coefficient_owners( tMesh, tMeshIndex, tNumCoeff );

        std::vector< int > tBandCoefficient( tNumCoeff, 0 );
        std::vector< int > tBandNode( tNumNodes, 0 );
        uint                tCutOld = 0, tCutNew = 0, tCutChanged = 0;
        for ( uint iCell = 0; iCell < tNumCells; iCell++ )
        {
            if ( tParallel && tMesh->get_mtk_cell( (moris_index)iCell ).get_owner() != par_rank() )
            {
                continue;
            }
            Matrix< IndexMat > tNodes = tMesh->get_nodes_connected_to_element_loc_inds( (moris_index)iCell );
            uint               tPosOld = 0, tPosNew = 0;
            for ( uint iNode = 0; iNode < tNodes.numel(); iNode++ )
            {
                tPosOld += tOldNodal( tNodes( iNode ) ) >= 0.0 ? 1 : 0;
                tPosNew += tNewNodal( tNodes( iNode ) ) >= 0.0 ? 1 : 0;
            }
            bool tIsCutOld = tPosOld > 0 && tPosOld < tNodes.numel();
            bool tIsCutNew = tPosNew > 0 && tPosNew < tNodes.numel();
            tCutOld += tIsCutOld ? 1 : 0;
            tCutNew += tIsCutNew ? 1 : 0;
            tCutChanged += tIsCutOld != tIsCutNew ? 1 : 0;
            if ( tIsCutOld )
            {
                for ( uint iC = 0; iC < tCellCoefficients( iCell ).size(); iC++ )
                {
                    tBandCoefficient[ tCellCoefficients( iCell )( iC ) ] = 1;
                }
                for ( uint iNode = 0; iNode < tNodes.numel(); iNode++ )
                {
                    tBandNode[ tNodes( iNode ) ] = 1;
                }
            }
        }

        // a cut cell of another processor can touch a coefficient or a node of this one: merge the band flags by ADV ID
        // and by node ID
        if ( tParallel )
        {
            moris_id tMaxNodeId = 0;
            for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
            {
                tMaxNodeId = std::max( tMaxNodeId, tMesh->get_glb_entity_id_from_entity_loc_index( (moris_index)iNode, mtk::EntityRank::NODE ) );
            }
            tMaxNodeId = max_all( tMaxNodeId );
            std::vector< int > tLocalNodeFlags( tMaxNodeId + 1, 0 );
            std::vector< int > tGlobalNodeFlags( tMaxNodeId + 1, 0 );
            for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
            {
                if ( tBandNode[ iNode ] == 1 )
                {
                    tLocalNodeFlags[ tMesh->get_glb_entity_id_from_entity_loc_index( (moris_index)iNode, mtk::EntityRank::NODE ) ] = 1;
                }
            }
            MPI_Allreduce( tLocalNodeFlags.data(), tGlobalNodeFlags.data(), (int)tLocalNodeFlags.size(), MPI_INT, MPI_MAX, get_comm() );
            for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
            {
                tBandNode[ iNode ] = tGlobalNodeFlags[ tMesh->get_glb_entity_id_from_entity_loc_index( (moris_index)iNode, mtk::EntityRank::NODE ) ];
            }

            Vector< sint > tADVIds = aGENPerformer( 0 )->get_design_determining_adv_ids( mADVFiledName );
            sint           tMinId  = MORIS_SINT_MAX;
            sint           tMaxId  = -MORIS_SINT_MAX;
            for ( uint iC = 0; iC < tNumCoeff; iC++ )
            {
                if ( tOwners( iC ) >= 0 )
                {
                    tMinId = std::min( tMinId, tADVIds( iC ) );
                    tMaxId = std::max( tMaxId, tADVIds( iC ) );
                }
            }
            tMinId = min_all( tMinId );
            tMaxId = max_all( tMaxId );

            std::vector< int > tLocalFlags( tMaxId - tMinId + 1, 0 );
            std::vector< int > tGlobalFlags( tMaxId - tMinId + 1, 0 );
            for ( uint iC = 0; iC < tNumCoeff; iC++ )
            {
                if ( tOwners( iC ) >= 0 && tBandCoefficient[ iC ] == 1 )
                {
                    tLocalFlags[ tADVIds( iC ) - tMinId ] = 1;
                }
            }
            MPI_Allreduce( tLocalFlags.data(), tGlobalFlags.data(), (int)tLocalFlags.size(), MPI_INT, MPI_MAX, get_comm() );
            for ( uint iC = 0; iC < tNumCoeff; iC++ )
            {
                if ( tOwners( iC ) >= 0 )
                {
                    tBandCoefficient[ iC ] = tGlobalFlags[ tADVIds( iC ) - tMinId ];
                }
            }
        }

        real tMaxBand = 0.0, tMaxFar = 0.0;
        uint tNumBand = 0, tNumOwnedCoeff = 0;
        for ( uint iC = 0; iC < tNumCoeff; iC++ )
        {
            if ( tParallel && tOwners( iC ) != par_rank() )
            {
                continue;
            }
            tNumOwnedCoeff++;
            real tDelta = std::abs( aNewCoeffs( iC ) - aOldCoeffs( iC ) );
            if ( tBandCoefficient[ iC ] == 1 )
            {
                tNumBand++;
                tMaxBand = std::max( tMaxBand, tDelta );
            }
            else
            {
                tMaxFar = std::max( tMaxFar, tDelta );
            }
        }

        uint tFlips = 0, tFlipsBand = 0, tNumOwnedNodes = 0;
        real tMaxNodalBand = 0.0;
        for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
        {
            if ( tParallel && tMesh->get_entity_owner( (moris_index)iNode, mtk::EntityRank::NODE ) != (uint)par_rank() )
            {
                continue;
            }
            tNumOwnedNodes++;
            bool tFlip = ( tOldNodal( iNode ) >= 0.0 ) != ( tNewNodal( iNode ) >= 0.0 );
            tFlips += tFlip ? 1 : 0;
            if ( tBandNode[ iNode ] == 1 )
            {
                tFlipsBand += tFlip ? 1 : 0;
                tMaxNodalBand = std::max( tMaxNodalBand, std::abs( tNewNodal( iNode ) - tOldNodal( iNode ) ) );
            }
        }

        MORIS_LOG_INFO( "SDF-FOOTPRINT coeffs=%u band_coeffs=%u max|dc|_band=%.6e max|dc|_far=%.6e max|dphi_node|_band=%.6e "
                        "nodes=%u sign_flips=%u (band %u) cut_cells psi=%u phi=%u changed=%u",
                sum_all( tNumOwnedCoeff ), sum_all( tNumBand ), max_all( tMaxBand ), max_all( tMaxFar ), max_all( tMaxNodalBand ),
                sum_all( tNumOwnedNodes ), sum_all( tFlips ), sum_all( tFlipsBand ), sum_all( tCutOld ), sum_all( tCutNew ), sum_all( tCutChanged ) );
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

            // a processor does not hold every neighbour of its shared and aura nodes: take their values from the owner
            synchronize_nodal_values( aCoefficients, aInterpMesh );
        }

        MORIS_LOG_INFO( "Laplacian smoothing complete. ADV range: [%.6f, %.6f]",
                aCoefficients.min(), aCoefficients.max() );
    }

}    // namespace moris::wrk
