/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_WRK_Workflow_HMR_XTK.cpp
 *
 */

#include "cl_WRK_Performer_Manager.hpp"
#include "cl_WRK_Workflow_HMR_XTK.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_HMR.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Mesh_Checker.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "cl_MTK_Cell.hpp"
#include "fn_WRK_Feature_Redetect.hpp"
#include "cl_MDL_Model.hpp"
#include "cl_WRK_GEN_Performer.hpp"

#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

#include "cl_Stopwatch.hpp"
#include "cl_WRK_perform_refinement.hpp"
#include "cl_WRK_perform_remeshing.hpp"

#include "fn_norm.hpp"

#include "cl_WRK_DataBase_Performer.hpp"
#include "cl_MIG.hpp"

#include "cl_WRK_Reinitialize_Performer.hpp"
#include "HDF5_Tools.hpp"

namespace moris::wrk
{
    namespace
    {
        // MB-089: whether two ADV vectors hold the same design (bit for bit)
        bool
        same_advs( const Vector< real >& aFirst, const Vector< real >& aSecond )
        {
            if ( aFirst.size() != aSecond.size() )
            {
                return false;
            }
            for ( uint iADV = 0; iADV < aFirst.size(); iADV++ )
            {
                if ( aFirst( iADV ) != aSecond( iADV ) )
                {
                    return false;
                }
            }
            return true;
        }

        // MB-097: processor 0's ADVs on every processor. The full ADV vector exists on processor 0 only: the optimizer
        // runs there, the other processors evaluate through opt::Algorithm::dummy_solve() with an empty vector, and
        // GEN's get_advs() is empty on them. Collective: every processor calls it.
        void
        broadcast_advs( Vector< real >& aADVs )
        {
            if ( par_size() == 1 )
            {
                return;
            }
            uint tNumADVs = aADVs.size();
            broadcast( tNumADVs );
            aADVs.resize( tNumADVs );
            MPI_Bcast( aADVs.memptr(), (int)tNumADVs, MPI_DOUBLE, 0, moris::get_comm() );
        }
    }    // namespace

    //--------------------------------------------------------------------------------------------------------------

    Workflow_HMR_XTK::Workflow_HMR_XTK( wrk::Performer_Manager* aPerformerManager )
            : Workflow( aPerformerManager )
    {
        // log & trace this function
        Tracer tTracer( "WRK", "HMR-XTK Workflow", "Create" );

        // Performer set for this workflow
        mPerformerManager->mHMRPerformer.resize( 1 );
        mPerformerManager->mGENPerformer.resize( 1 );
        mPerformerManager->mXTKPerformer.resize( 1 );
        mPerformerManager->mMTKPerformer.resize( 2 );
        mPerformerManager->mMDLPerformer.resize( 1 );
        mPerformerManager->mRemeshingMiniPerformer.resize( 1 );
        mPerformerManager->mDataBasePerformer.resize( 1 );

        // load parameter lists
        Module_Parameter_Lists tHMRParameterList   = aPerformerManager->mLibrary->get_parameters_for_module( Module_Type::HMR );
        Module_Parameter_Lists tGENParameterList   = aPerformerManager->mLibrary->get_parameters_for_module( Module_Type::GEN );
        Module_Parameter_Lists tMORISParameterList = aPerformerManager->mLibrary->get_parameters_for_module( Module_Type::MORISGENERAL );

        // create HMR performer
        mPerformerManager->mHMRPerformer( 0 ) = std::make_shared< hmr::HMR >( tHMRParameterList, mPerformerManager->mLibrary );

        // create MTK performer - will be used for HMR mesh
        mPerformerManager->mMTKPerformer( 0 ) = std::make_shared< mtk::Mesh_Manager >();

        // Create GE performer
        mPerformerManager->mGENPerformer( 0 ) = std::make_shared< gen::Geometry_Engine >(
                tGENParameterList,
                mPerformerManager->mLibrary );

        // create MTK performer - will be used for XTK mesh
        mPerformerManager->mMTKPerformer( 1 ) = std::make_shared< mtk::Mesh_Manager >();

        // create MDL performer
        mPerformerManager->mMDLPerformer( 0 ) = std::make_shared< mdl::Model >( mPerformerManager->mLibrary, 0 );

        // Set MTK performer to HMR
        mPerformerManager->mHMRPerformer( 0 )->set_performer( mPerformerManager->mMTKPerformer( 0 ) );

        // Set performer to MDL
        mPerformerManager->mMDLPerformer( 0 )->set_performer( mPerformerManager->mMTKPerformer( 1 ) );

        if ( tMORISParameterList.size() != 0 )
        {
            if ( !tMORISParameterList( 0 ).empty() )
            {
                // create re-meshing performer - will be used for HMR mesh
                mPerformerManager->mRemeshingMiniPerformer( 0 ) =
                        std::make_shared< wrk::Remeshing_Mini_Performer >( tMORISParameterList( 0 )( 0 ), mPerformerManager->mLibrary );
            }

            if ( !tMORISParameterList( 2 ).empty() )
            {
                // allocate size
                mPerformerManager->mReinitializePerformer.resize( 1 );

                // construct the parameter list
                mPerformerManager->mReinitializePerformer( 0 ) = std::make_shared< wrk::Reinitialize_Performer >( mPerformerManager->mLibrary );
            }
        }

        // Additive discretization (field + B-spline correction) cannot be remapped by remeshing or reinitialization
        this->check_additive_discretization_compatibility();
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Workflow_HMR_XTK::gen_refines_mesh()
    {
        std::shared_ptr< gen::Geometry_Engine > tGEN = mPerformerManager->mGENPerformer( 0 );
        for ( uint iField = 0; iField < tGEN->get_num_refinement_fields(); iField++ )
        {
            for ( uint iNum : tGEN->get_num_refinements( iField ) )
            {
                if ( iNum > 0 )
                {
                    return true;
                }
            }
        }
        return false;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Workflow_HMR_XTK::refine_hmr_on_current_design()
    {
        Tracer tTracer( "WRK", "HMR-XTK Workflow", "Rebuild HMR with GEN refinement on current ADVs" );

        // a new background mesh from the parameter list, so that cells refined for the earlier design are coarse again
        Module_Parameter_Lists tHMRParameterList = mPerformerManager->mLibrary->get_parameters_for_module( Module_Type::HMR );
        mPerformerManager->mHMRPerformer( 0 ) = std::make_shared< hmr::HMR >( tHMRParameterList, mPerformerManager->mLibrary );
        mPerformerManager->mMTKPerformer( 0 ) = std::make_shared< mtk::Mesh_Manager >();
        mPerformerManager->mHMRPerformer( 0 )->set_performer( mPerformerManager->mMTKPerformer( 0 ) );

        // a new GEN from its parameter list, at the current ADVs: after the first distribution, fields stored on the old
        // mesh (discretization_mesh_index = -1, e.g. fixed pads) can no longer be evaluated on the new mesh, and the
        // MB-085 equalisation only sees fields that are still analytic; a new GEN is in the state of the first refinement.
        // Every processor refines with the full ADV vector (MB-097)
        Vector< real >                          tCurrentADVs      = this->get_current_advs_on_all_procs();
        Module_Parameter_Lists tGENParameterList = mPerformerManager->mLibrary->get_parameters_for_module( Module_Type::GEN );
        std::shared_ptr< gen::Geometry_Engine > tNewGEN = std::make_shared< gen::Geometry_Engine >( tGENParameterList, mPerformerManager->mLibrary );
        MORIS_ERROR( tCurrentADVs.size() == tNewGEN->get_num_primitive_advs(),
                "Workflow_HMR_XTK::refine_hmr_on_current_design() - %zu ADVs, %u primitive ADVs: designs discretized on the "
                "remeshed B-spline mesh need a remeshing performer.",
                tCurrentADVs.size(),
                tNewGEN->get_num_primitive_advs() );
        tNewGEN->set_initial_primitive_advs( tCurrentADVs );
        mPerformerManager->mGENPerformer( 0 ) = tNewGEN;

        // uniform initial refinement, then the GEN-driven refinement with GEN's current ADVs
        mPerformerManager->mHMRPerformer( 0 )->perform_initial_refinement();

        Refinement_Mini_Performer tRefinementPerformer;
        std::shared_ptr< Performer > tGenInterfacePerformer = std::make_shared< wrk::Gen_Performer >( mPerformerManager->mGENPerformer( 0 ) );
        tRefinementPerformer.perform_refinement_old( mPerformerManager->mHMRPerformer( 0 ), { tGenInterfacePerformer } );

        mPerformerManager->mHMRPerformer( 0 )->perform();

        // the MDL fields belong to the replaced mesh
        mMDLFieldsOnCurrentMesh = false;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Workflow_HMR_XTK::apply_restart_advs_before_refinement()
    {
        std::shared_ptr< gen::Geometry_Engine > tGEN = mPerformerManager->mGENPerformer( 0 );
        if ( tGEN->advs_distributed() or tGEN->get_num_primitive_advs() == 0 or not this->gen_refines_mesh() )
        {
            return;
        }

        Module_Parameter_Lists tOPTParameterList = mPerformerManager->mLibrary->get_parameters_for_module( Module_Type::OPT );
        if ( tOPTParameterList.size() == 0 or tOPTParameterList( 0 ).size() == 0
                or not tOPTParameterList( 0 )( 0 ).get< bool >( "is_optimization_problem" ) )
        {
            return;
        }
        std::string tRestartFile = tOPTParameterList( 0 )( 0 ).get< std::string >( "restart_file" );
        if ( tRestartFile.empty() )
        {
            return;
        }

        // every rank reads the file (read only, no parallel extension), as GEN refines on every rank
        hid_t          tFileID = open_hdf5_file( tRestartFile, false, true );
        Vector< real > tRestartADVs;
        herr_t         tStatus = 0;
        load_vector_from_hdf5_file( tFileID, "ADVs", tRestartADVs.data(), tStatus );
        close_hdf5_file( tFileID );

        // the restart file holds the full ADV vector; the primitive ADVs come first. ADVs on a B-spline mesh do not
        // exist before the refinement (their designs are refined on the initial field, unchanged here)
        uint tNumPrimitive = tGEN->get_num_primitive_advs();
        MORIS_ERROR( tRestartADVs.size() >= tNumPrimitive,
                "Workflow_HMR_XTK - restart file %s holds %zu ADVs, the deck %u primitive ADVs.",
                tRestartFile.c_str(),
                tRestartADVs.size(),
                tNumPrimitive );

        Vector< real > tPrimitive( tNumPrimitive );
        for ( uint iADV = 0; iADV < tNumPrimitive; iADV++ )
        {
            tPrimitive( iADV ) = tRestartADVs( iADV );
        }
        real tMaxChange = tGEN->set_initial_primitive_advs( tPrimitive );

        MORIS_LOG_INFO( "Workflow_HMR_XTK - GEN refinement of the initial mesh uses the %u primitive ADVs of the restart file "
                        "%s (max change from the deck's initial values %g; MB-086)",
                tNumPrimitive,
                tRestartFile.c_str(),
                tMaxChange );
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< sint >
    Workflow_HMR_XTK::get_feature_redetect_thresholds( uint aNumADVs )
    {
        std::shared_ptr< gen::Geometry_Engine > tGEN = mPerformerManager->mGENPerformer( 0 );

        Vector< sint > tThresholds( tGEN->get_number_of_geometries(), -1 );
        if ( mPerformerManager->mRemeshingMiniPerformer( 0 ) or not this->gen_refines_mesh()
                or tGEN->get_num_primitive_advs() != aNumADVs )
        {
            return tThresholds;
        }
        for ( uint iGeom = 0; iGeom < tThresholds.size(); iGeom++ )
        {
            tThresholds( iGeom ) = tGEN->get_feature_refinement_redetect_edges( iGeom );
        }
        return tThresholds;
    }

    //--------------------------------------------------------------------------------------------------------------

    uint
    Workflow_HMR_XTK::get_num_advs_on_all_procs( const Vector< real >& aADVs ) const
    {
        uint tNumADVs = aADVs.size();
        broadcast( tNumADVs );
        return tNumADVs;
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< real >
    Workflow_HMR_XTK::get_current_advs_on_all_procs()
    {
        Vector< real > tADVs = mPerformerManager->mGENPerformer( 0 )->get_advs();
        broadcast_advs( tADVs );
        return tADVs;
    }

    //--------------------------------------------------------------------------------------------------------------

    std::pair< uint, uint >
    Workflow_HMR_XTK::count_refined_background_cells()
    {
        mtk::Interpolation_Mesh* tMesh = mPerformerManager->mMTKPerformer( 0 )->get_interpolation_mesh( 0 );

        uint tRefined = 0;
        uint tTotal   = 0;
        for ( uint iCell = 0; iCell < tMesh->get_num_elems(); iCell++ )
        {
            const mtk::Cell& tCell = tMesh->get_mtk_cell( (moris_index)iCell );
            if ( tCell.get_owner() != (moris_id)par_rank() )
            {
                continue;
            }
            tTotal++;
            tRefined += tCell.get_level() > 0 ? 1 : 0;
        }
        return { sum_all( tRefined ), sum_all( tTotal ) };
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Workflow_HMR_XTK::check_additive_discretization_compatibility()
    {
        // Same activity conditions as perform(): remeshing performer built / reinitialize performer built
        bool tRemeshingActive      = mPerformerManager->mRemeshingMiniPerformer( 0 ) != nullptr;
        bool tReinitializingActive = mPerformerManager->mReinitializePerformer.size() > 0;

        if ( not( tRemeshingActive or tReinitializingActive ) )
        {
            return;
        }

        Vector< std::string > tAdditiveDesigns = mPerformerManager->mGENPerformer( 0 )->get_additive_discretized_design_names();

        if ( tAdditiveDesigns.size() > 0 )
        {
            std::string tDesignList;
            for ( const std::string& iName : tAdditiveDesigns )
            {
                tDesignList += ( tDesignList.empty() ? "" : ", " ) + iName;
            }

            MORIS_ERROR( false,
                    "Workflow_HMR_XTK - design(s) [%s] use discretization_mode = \"additive\", which is incompatible with %s%s%s: "
                    "the analytic field + B-spline correction split is not preserved when the field is remapped. "
                    "Use discretization_mode = \"replace\" or disable remeshing / reinitialization.",
                    tDesignList.c_str(),
                    tRemeshingActive ? "HMR remeshing" : "",
                    tRemeshingActive and tReinitializingActive ? " and " : "",
                    tReinitializingActive ? "level-set reinitialization" : "" );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Workflow_HMR_XTK::initialize(
            Vector< real >&  aADVs,
            Vector< real >&  aLowerBounds,
            Vector< real >&  aUpperBounds,
            Matrix< IdMat >& aIjklIDs )
    {
        // Trace & log this function
        Tracer tTracer( "WRK", "HMR-XTK Workflow", "Initialize" );

        mInitializeOptimizationRestart = false;
        
        // If mSDFReinitPending was true, we're coming from an SDF reinit restart.
        // Set skip flag to prevent re-triggering SDF reinit on the restart pass.
        // Also apply the computed SDF coefficients to GEN BEFORE distribute_advs().
        if ( mSDFReinitPending )
        {
            mSkipNextSDFReinit = true;
            
            // Apply the stored SDF coefficients to GEN
            if ( mPerformerManager->mReinitializePerformer.size() > 0 )
            {
                Reinitialize_Performer* tReinitPerformer = mPerformerManager->mReinitializePerformer( 0 ).get();
                Matrix< DDRMat > const& tCoeffs = tReinitPerformer->get_coefficients();
                
                if ( tCoeffs.n_rows() > 0 )
                {
                    // Get current ADV count from GEN to verify size matches
                    Vector< real > const& tCurrentADVs = mPerformerManager->mGENPerformer( 0 )->get_advs();
                    uint tExpectedSize = tCurrentADVs.size();
                    uint tStoredSize = tCoeffs.n_rows();
                    
                    if ( tStoredSize != tExpectedSize )
                    {
                        MORIS_LOG_WARNING( "SDF coefficient count (%d) does not match current ADV count (%d). "
                                "This can happen if mesh topology changed during reinit. Skipping SDF application.",
                                (int)tStoredSize, (int)tExpectedSize );
                    }
                    else
                    {
                        // Convert to Vector for ADVs
                        Vector< real > tNewADVs( tCoeffs.n_rows() );
                        for ( uint i = 0; i < tCoeffs.n_rows(); ++i )
                        {
                            tNewADVs( i ) = tCoeffs( i );
                        }
                        
                        // Set new ADVs in GEN - these will be picked up by distribute_advs() later
                        mPerformerManager->mGENPerformer( 0 )->set_advs( tNewADVs );
                        
                        MORIS_LOG_INFO( "Applied %d SDF coefficients to GEN in initialize()", (int)tCoeffs.n_rows() );
                    }
                }
            }
        }
        
        // Clear SDF reinit pending flag (coefficients now applied above)
        mSDFReinitPending = false;

        // sdf_reinit_mode "user_function": this restart was requested by the reinitialization schedule (perform()). GEN still
        // holds the ADVs of the last evaluated design (the trial point that triggered the restart was never set); replace the
        // targeted design's ADVs by the deck function's result BEFORE the mesh and GEN are rebuilt below, so that the GEN
        // refinement, distribute_advs() and the ADVs handed to the optimizer all belong to the new design. The optimizer then
        // re-initializes its state at it (asymptotes at the bounds, move limits at their initial value, objective scale
        // re-derived) and evaluates it as its first iterate, so criteria and sensitivities are those of the new design.
        if ( mUserReinitPending )
        {
            mUserReinitPending = false;

            MORIS_ERROR( mPerformerManager->mReinitializePerformer.size() > 0
                                 and mPerformerManager->mReinitializePerformer( 0 )->get_reinit_mode() == ReinitMode::USER_FUNCTION,
                    "Workflow_HMR_XTK::initialize() - user_function reinitialization pending without a user_function performer." );

            mPerformerManager->mReinitializePerformer( 0 )->apply_user_function( mPerformerManager->mGENPerformer, mUserReinitIteration );

            mUserReinitCount++;
            MORIS_LOG_INFO(
                    "Workflow_HMR_XTK::initialize() - user_function reinitialization %u applied (optimization iteration %u); "
                    "the optimizer restarts at the new design",
                    mUserReinitCount,
                    mUserReinitIteration );
        }

        mIter = 0;

        Vector< std::shared_ptr< mtk::Field > > tFieldsIn;
        Vector< std::shared_ptr< mtk::Field > > tFieldsOut;

        // perform in first optimization iteration only
        if ( tIsFirstOptSolve )
        {
            // Step 1: HMR refinement -------------------------------------------------------------------

            // mPerformerManager->mHMRPerformer( 0 )->reset_HMR();

            // check whether HMR should be NOT initialized from restart file
            if ( not mPerformerManager->mHMRPerformer( 0 )->get_restarted_from_file() )
            {
                // uniform initial refinement
                mPerformerManager->mHMRPerformer( 0 )->perform_initial_refinement();

                // MB-086: refine around the restart design, not the deck's initial ADVs
                this->apply_restart_advs_before_refinement();

                // HMR refined by GE
                Refinement_Mini_Performer tRefinementPerformer;

                // get the GEN interface performer
                std::shared_ptr< gen::Geometry_Engine > tGeometryEngine = mPerformerManager->mGENPerformer( 0 );

                std::shared_ptr< Performer > tGenInterfacePerformer = std::make_shared< wrk::Gen_Performer >( tGeometryEngine );

                // perform initial refinement for GEN and user defined refinement functions
                tRefinementPerformer.perform_refinement_old( mPerformerManager->mHMRPerformer( 0 ), { tGenInterfacePerformer } );
            }

            // HMR finalize
            mPerformerManager->mHMRPerformer( 0 )->perform();

            // switch flag for first solve in optimization process
            tIsFirstOptSolve = false;
        }
        // after first solve in optimization process, no remeshing performer, GEN-driven refinement: rebuild the mesh with
        // the GEN refinement of the current ADVs (e.g. after an XTK failure because the design left the refined band).
        // MB-097: the ADV count of processor 0 (the only one holding the full vector), so that every processor takes the
        // same branch of these collective rebuilds
        else if ( not mPerformerManager->mRemeshingMiniPerformer( 0 ) and this->gen_refines_mesh()
                  and mPerformerManager->mGENPerformer( 0 )->get_num_primitive_advs()
                              == this->get_num_advs_on_all_procs( mPerformerManager->mGENPerformer( 0 )->get_advs() ) )
        {
            this->refine_hmr_on_current_design();
        }
        // after the first solve, no remeshing performer and no GEN refinement, primitive ADVs only: the background mesh does
        // not depend on the design, so it is kept; GEN is rebuilt at the current ADVs as in refine_hmr_on_current_design()
        // (without the HMR rebuild). Any optimizer restart of such a deck (reinitialize_interface_iter, a user_function
        // reinitialization, an XTK failure, a non-finite MMA step) used to stop here: "remeshing performer has not been built".
        else if ( not mPerformerManager->mRemeshingMiniPerformer( 0 )
                  and mPerformerManager->mGENPerformer( 0 )->get_num_primitive_advs()
                              == this->get_num_advs_on_all_procs( mPerformerManager->mGENPerformer( 0 )->get_advs() ) )
        {
            Tracer tRebuildTracer( "WRK", "HMR-XTK Workflow", "Rebuild GEN on current ADVs" );

            Vector< real >         tCurrentADVs      = this->get_current_advs_on_all_procs();
            Module_Parameter_Lists tGENParameterList = mPerformerManager->mLibrary->get_parameters_for_module( Module_Type::GEN );

            std::shared_ptr< gen::Geometry_Engine > tNewGEN =
                    std::make_shared< gen::Geometry_Engine >( tGENParameterList, mPerformerManager->mLibrary );
            tNewGEN->set_initial_primitive_advs( tCurrentADVs );
            mPerformerManager->mGENPerformer( 0 ) = tNewGEN;
        }
        // after first solve in optimization process
        else
        {
            // get access to the performers
            std::shared_ptr< gen::Geometry_Engine >         tGeometryEngine     = mPerformerManager->mGENPerformer( 0 );
            std::shared_ptr< mdl::Model >                   tModelPerformer     = mPerformerManager->mMDLPerformer( 0 );
            std::shared_ptr< Remeshing_Mini_Performer >     tRemeshingPerformer = mPerformerManager->mRemeshingMiniPerformer( 0 );
            Vector< std::shared_ptr< hmr::HMR > >&          tHmrPerformers      = mPerformerManager->mHMRPerformer;
            Vector< std::shared_ptr< mtk::Mesh_Manager > >& tMtkPerformers      = mPerformerManager->mMTKPerformer;

            // get refinement fields from GEN and MDL performers; the MDL fields only if a solve ran on this mesh
            tFieldsIn.append( tGeometryEngine->get_mtk_fields() );
            if ( mMDLFieldsOnCurrentMesh )
            {
                tFieldsIn.append( tModelPerformer->get_mtk_fields() );
            }
            else
            {
                MORIS_LOG_WARNING( "Workflow_HMR_XTK::initialize() - no forward solve on the current mesh since the last remesh; "
                                   "remeshing on the geometry fields only" );
            }

            // check remeshing mini-performer has been built
            MORIS_ERROR(
                    tRemeshingPerformer,
                    "Workflow_HMR_XTK::initialize() - remeshing performer has not been built." );

            // refine meshes
            tRemeshingPerformer->perform_remeshing(
                    tFieldsIn,
                    tHmrPerformers,
                    tMtkPerformers,
                    tFieldsOut );

            // the MDL fields now live on the replaced mesh
            mMDLFieldsOnCurrentMesh = false;

            // get access to the library
            std::shared_ptr< Library_IO > tLibrary = mPerformerManager->mLibrary;

            // re-initialize GEN
            Module_Parameter_Lists tGENParameterList = tLibrary->get_parameters_for_module( Module_Type::GEN );
            tGeometryEngine                       = std::make_shared< gen::Geometry_Engine >( tGENParameterList, tLibrary );
        }

        // Step 2: Initialize Level set field in GEN -----------------------------------------------
        {
            // retrieve the mesh pair
            const mtk::Mesh_Pair& tMeshPair = mPerformerManager->mMTKPerformer( 0 )->get_mesh_pair( 0 );

            // initialize GEN
            std::shared_ptr< gen::Geometry_Engine > tGeometryEngine = mPerformerManager->mGENPerformer( 0 );
            tGeometryEngine->distribute_advs( tMeshPair, tFieldsOut );

            // Get ADVs
            aADVs        = tGeometryEngine->get_advs();
            aLowerBounds = tGeometryEngine->get_lower_bounds();
            aUpperBounds = tGeometryEngine->get_upper_bounds();
            aIjklIDs     = tGeometryEngine->get_IjklIDs();
        }

        // MB-089: the mesh is refined for these ADVs; the first cut of this design gives the resolution baseline
        mRefinedForADVs = aADVs;
        mRedetectBaselineEdges.clear();

    }    // end function: Workflow_HMR_XTK::initialize()

    //--------------------------------------------------------------------------------------------------------------

    Vector< real >
    Workflow_HMR_XTK::perform( Vector< real >& aNewADVs )
    {
        // get optimization iteration
        sint tOptIter = gLogger.get_iteration( "OPT", "Manager", "Perform" );

        // compute total optimization iteration accounting for iterations performed //
        // in previous instances of optimization algorithm
        tOptIter = tOptIter + mIter;

        // return vector of design criteria with NANs causing the optimization algorithm to restart
        // Skip this block for INTERFACE_SDF mode - we handle NANs in the SDF reinit block (line ~420)
        bool tIsInterfaceSDFMode = mPerformerManager->mReinitializePerformer.size() > 0 &&
                mPerformerManager->mReinitializePerformer( 0 )->get_reinit_mode() == ReinitMode::INTERFACE_SDF;
        
        if ( !tIsInterfaceSDFMode && ( mIter >= mReinitializeIterInterval or ( uint ) tOptIter == mReinitializeIterFirst ) )
        {
            mInitializeOptimizationRestart = true;

            // user_function mode: initialize() applies the deck function on this restart (to the last evaluated design)
            if ( mPerformerManager->mReinitializePerformer.size() > 0
                    and mPerformerManager->mReinitializePerformer( 0 )->get_reinit_mode() == ReinitMode::USER_FUNCTION )
            {
                mUserReinitPending   = true;
                mUserReinitIteration = gLogger.get_opt_iteration();
                MORIS_LOG_INFO(
                        "Workflow_HMR_XTK::perform() - user_function reinitialization requested after %u evaluations "
                        "since the last (re)initialization (optimization iteration %u); returning NaN criteria to "
                        "restart the optimizer",
                        mIter,
                        mUserReinitIteration );
            }

            Vector< real > tVector( mNumCriteria, std::numeric_limits< real >::quiet_NaN() );

            return tVector;
        }

        // Stage *: Re-initialization of the adv field
        if ( mPerformerManager->mReinitializePerformer.size() > 0 )
        {
            Reinitialize_Performer* tReinitPerformer = mPerformerManager->mReinitializePerformer( 0 ).get();
            sint tReinitFreq = tReinitPerformer->get_reinitialization_frequency();

            // FIELD_REMAP mode: reinitialize ADVs from solution field before XTK
            if ( tReinitPerformer->get_reinit_mode() == ReinitMode::FIELD_REMAP 
                 && tOptIter > 0 && tOptIter % tReinitFreq == 0 )
            {
                // Set new advs in GE
                Tracer tTracer( "GEN", "Levelset", "Re-InitializeADVs" );

                mPerformerManager->mGENPerformer( 0 )->distribute_advs(
                        mPerformerManager->mMTKPerformer( 0 )->get_mesh_pair( 0 ),
                        mPerformerManager->mReinitializePerformer( 0 )->get_mtk_fields() );

                // get advs from GE and overwrite them
                aNewADVs = mPerformerManager->mGENPerformer( 0 )->get_advs();
            }
            // INTERFACE_SDF mode: Handle restart after SDF reinit triggered NANs
            // Apply the stored SDF coefficients when mSDFReinitPending is set
            else if ( tReinitPerformer->get_reinit_mode() == ReinitMode::INTERFACE_SDF
                      && mSDFReinitPending )
            {
                Tracer tTracer( "GEN", "Levelset", "Interface-SDF-Apply-ADVs" );
                
                // Get the stored SDF coefficients and use them as new ADVs
                Matrix< DDRMat > const& tCoeffs = tReinitPerformer->get_coefficients();
                
                if ( tCoeffs.n_rows() > 0 )
                {
                    // Convert to Vector for ADVs
                    Vector< real > tNewADVs( tCoeffs.n_rows() );
                    for ( uint i = 0; i < tCoeffs.n_rows(); ++i )
                    {
                        tNewADVs( i ) = tCoeffs( i );
                    }
                    
                    // Set new ADVs in GEN AND update aNewADVs
                    // This is safe on restart - optimizer hasn't computed criteria yet
                    mPerformerManager->mGENPerformer( 0 )->set_advs( tNewADVs );
                    aNewADVs = tNewADVs;
                    
                    MORIS_LOG_INFO( "Applied %zu SDF coefficients as new ADVs", tCoeffs.n_rows() );
                }
                
                // Clear the pending flag - ADVs have been applied
                mSDFReinitPending = false;
            }
            else
            {
                mPerformerManager->mGENPerformer( 0 )->set_advs( aNewADVs );
            }
        }
        else
        {
            // Set new advs in GE
            mPerformerManager->mGENPerformer( 0 )->set_advs( aNewADVs );
        }

        // Stage 1: HMR refinement
        if ( mPerformerManager->mRemeshingMiniPerformer( 0 ) )
        {
            sint tReinitFreq = mPerformerManager->mRemeshingMiniPerformer( 0 )->get_reinitialization_frequency();

            if ( tOptIter > 0 and tOptIter % tReinitFreq == 0 )
            {
                // remeshing rebuilds the B-spline basis, but the optimizer's ADVs are the unfiltered variables psi of
                // the old basis; there is no mapping of psi to the new basis
                MORIS_ERROR( not mPerformerManager->mGENPerformer( 0 )->has_adv_filter(),
                        "Workflow_HMR_XTK: remeshing during the optimization is not supported with adv_filter_radius > 0" );

                // allocate auxiliary arrays
                Vector< real >  tADVs;
                Vector< real >  tLowerBounds;
                Vector< real >  tUpperBounds;
                Matrix< IdMat > tIjklIDs;

                // initialize HMR and GEN
                this->initialize( tADVs, tLowerBounds, tUpperBounds, tIjklIDs );

                // Set new advs in GE
                mPerformerManager->mGENPerformer( 0 )->set_advs( aNewADVs );
            }
        }

        // Stage 2: XTK -----------------------------------------------------------------------------

        // Read parameter list from shared object
        Module_Parameter_Lists tXTKParameterList = mPerformerManager->mLibrary->get_parameters_for_module( Module_Type::XTK );

        // Create XTK
        xtk::Model* tXTKPerformer = new xtk::Model( tXTKParameterList( 0 )( 0 ) );

        std::shared_ptr< mtk::Mesh_Manager > tMTKPerformer = std::make_shared< mtk::Mesh_Manager >();

        // Set performers
        tXTKPerformer->set_geometry_engine( mPerformerManager->mGENPerformer( 0 ).get() );
        tXTKPerformer->set_input_performer( mPerformerManager->mMTKPerformer( 0 ) );
        tXTKPerformer->set_output_performer( tMTKPerformer );

        // Compute level set data in GEN
        // FIXME: HMR stores mesh with aura on 0
        mPerformerManager->mGENPerformer( 0 )->reset_mesh_information(
                mPerformerManager->mMTKPerformer( 0 )->get_interpolation_mesh( 0 ) );

        // Output GEN fields, if requested
        mPerformerManager->mGENPerformer( 0 )->output_fields(
                mPerformerManager->mMTKPerformer( 0 )->get_interpolation_mesh( 0 ) );

        // mtk::Mesh_Checker tMeshCheckerHMR(
        //         0,
        //         mPerformerManager->mMTKPerformer( 0 )->get_interpolation_mesh(0),
        //         mPerformerManager->mMTKPerformer( 0 )->get_integration_mesh(0));
        // tMeshCheckerHMR.perform();
        // tMeshCheckerHMR.print_diagnostics();

        bool tDeleteXTK = tXTKPerformer->delete_xtk_after_generation();
        // XTK perform - decompose - enrich - ghost - multigrid
        bool tFlag = tXTKPerformer->perform_decomposition();

        // MB-089: re-detection of sub-cell features. The GEN refinement (MB-084) was computed for an earlier design; when
        // the XTK geometry resolution check of this cut counts more hidden edges than the deck threshold
        // (feature_refinement_redetect_edges) and than the first cut after the last refinement, redo the refinement on
        // the current ADVs through the MB-087 path below, as if the decomposition had failed. The counts are those of
        // the check that already ran in this decomposition (no extra field evaluations).
        // MB-097: the redo below is collective (HMR and GEN are rebuilt on every processor), so every input of the
        // decision is the same on every processor: aNewADVs is complete on processor 0 only (the other processors
        // evaluate through opt::Algorithm::dummy_solve() with an empty vector), hence processor 0's ADV count and its
        // comparison with the design the mesh was refined for; tFlag and the hidden-edge counts are global already
        uint                      tNumADVs            = this->get_num_advs_on_all_procs( aNewADVs );
        Vector< sint >            tRedetectThresholds = this->get_feature_redetect_thresholds( tNumADVs );
        Vector< uint >            tHiddenEdges;
        uint                      tRefinedForThisDesignMsg = same_advs( aNewADVs, mRefinedForADVs ) ? 1 : 0;
        broadcast( tRefinedForThisDesignMsg );
        bool                      tMeshRefinedForThisDesign = tRefinedForThisDesignMsg == 1;
        Feature_Redetect_Decision tRedetect;
        if ( tFlag and tXTKPerformer->get_cut_integration_mesh() )
        {
            tHiddenEdges = tXTKPerformer->get_cut_integration_mesh()->mResolutionHiddenEdges;
            tRedetect    = decide_feature_redetect( tRedetectThresholds, tHiddenEdges, mRedetectBaselineEdges, tMeshRefinedForThisDesign );
        }

        // The design left the band GEN refined for an earlier design: cut cells on several levels, decomposition fails.
        // Without a remeshing performer, redo the GEN refinement (incl. the MB-084 sub-cell feature refinement) on the
        // current ADVs and cut again, in this evaluation. The number of ADVs does not change for designs that are not
        // discretized on the remeshed B-spline mesh, so the optimizer is not restarted.
        // MB-097: processor 0 decides, every processor follows (the inputs agree; this keeps any future rank-local input
        // from splitting the processors at the collectives of the redo)
        bool tRedo = ( not tFlag or tRedetect.mRedo ) and not mPerformerManager->mRemeshingMiniPerformer( 0 ) and this->gen_refines_mesh()
                 and mPerformerManager->mGENPerformer( 0 )->get_num_primitive_advs() == tNumADVs;
        uint tRedoMsg = tRedo ? 1 : 0;
        broadcast( tRedoMsg );
        if ( tRedoMsg == 1 )
        {
            Tracer tRefineTracer( "WRK", "HMR-XTK Workflow", "Refine on current design" );

            mNumRefinementRedos++;

            std::pair< uint, uint > tCellsBefore = this->count_refined_background_cells();

            if ( not tFlag )
            {
                MORIS_LOG_INFO( "Workflow_HMR_XTK::perform() - decomposition failed on the mesh refined for an earlier design; "
                                "redoing the GEN refinement on the current ADVs and cutting again" );
            }
            else
            {
                std::string tBaseline = tRedetect.mGeometry < mRedetectBaselineEdges.size()
                                              ? std::to_string( mRedetectBaselineEdges( tRedetect.mGeometry ) )
                                              : std::string( "none yet" );
                MORIS_LOG_INFO( "Workflow_HMR_XTK::perform() - optimization iteration %d: geometry %u under-resolved on the mesh "
                                "refined for an earlier design (MB-089): %u background edges cross the interface more often "
                                "than the cut, threshold %u (feature_refinement_redetect_edges %d; first cut after the last "
                                "refinement: %s); redoing the GEN refinement on the current ADVs and cutting again",
                        (int)tOptIter,
                        tRedetect.mGeometry,
                        tRedetect.mHiddenEdges,
                        tRedetect.mThreshold,
                        (int)tRedetectThresholds( tRedetect.mGeometry ),
                        tBaseline.c_str() );
            }

            delete tXTKPerformer;

            Vector< real >  tADVs;
            Vector< real >  tLowerBounds;
            Vector< real >  tUpperBounds;
            Matrix< IdMat > tIjklIDs;

            // initialize() restarts the count of evaluations since the last re-initialization (mIter). Neither a
            // re-detection (MB-089) nor a redo after a failed decomposition (MB-087) re-initializes the optimization, so
            // both keep the count: resetting it on the MB-087 path postponed the OPT reinitialization schedule
            // (reinitialize_interface_iter) indefinitely while decompositions kept failing, e.g. as topology forms (MB-093)
            uint tIterBeforeRedo = mIter;
            this->initialize( tADVs, tLowerBounds, tUpperBounds, tIjklIDs );
            mIter           = tIterBeforeRedo;
            mRefinedForADVs = aNewADVs;

            std::pair< uint, uint > tCellsAfter = this->count_refined_background_cells();
            MORIS_LOG_INFO( "Workflow_HMR_XTK::perform() - optimization iteration %d: GEN refinement redone on the current "
                            "ADVs: %u of %u background (leaf) cells on a refined level, before %u of %u",
                    (int)tOptIter,
                    tCellsAfter.first,
                    tCellsAfter.second,
                    tCellsBefore.first,
                    tCellsBefore.second );

            mPerformerManager->mGENPerformer( 0 )->set_advs( aNewADVs );

            tXTKPerformer = new xtk::Model( tXTKParameterList( 0 )( 0 ) );
            tMTKPerformer = std::make_shared< mtk::Mesh_Manager >();

            tXTKPerformer->set_geometry_engine( mPerformerManager->mGENPerformer( 0 ).get() );
            tXTKPerformer->set_input_performer( mPerformerManager->mMTKPerformer( 0 ) );
            tXTKPerformer->set_output_performer( tMTKPerformer );

            mPerformerManager->mGENPerformer( 0 )->reset_mesh_information(
                    mPerformerManager->mMTKPerformer( 0 )->get_interpolation_mesh( 0 ) );

            mPerformerManager->mGENPerformer( 0 )->output_fields(
                    mPerformerManager->mMTKPerformer( 0 )->get_interpolation_mesh( 0 ) );

            tDeleteXTK = tXTKPerformer->delete_xtk_after_generation();
            tFlag      = tXTKPerformer->perform_decomposition();

            if ( tFlag and tXTKPerformer->get_cut_integration_mesh() )
            {
                tHiddenEdges = tXTKPerformer->get_cut_integration_mesh()->mResolutionHiddenEdges;
                if ( tRedetect.mRedo and tRedetect.mGeometry < tHiddenEdges.size() )
                {
                    MORIS_LOG_INFO( "Workflow_HMR_XTK::perform() - optimization iteration %d: geometry %u after the redo "
                                    "(MB-089): %u background edges cross the interface more often than the cut (%u before)",
                            (int)tOptIter,
                            tRedetect.mGeometry,
                            tHiddenEdges( tRedetect.mGeometry ),
                            tRedetect.mHiddenEdges );
                }
            }
            tMeshRefinedForThisDesign = true;
        }

        // MB-089: the first cut on a mesh refined for the design it cuts sets the resolution baseline
        if ( tFlag and tMeshRefinedForThisDesign and mRedetectBaselineEdges.size() == 0 )
        {
            mRedetectBaselineEdges = tHiddenEdges;
        }

        if ( not tFlag )
        {
            mInitializeOptimizationRestart = true;

            MORIS_ERROR( mNumCriteria != MORIS_UINT_MAX,
                    "Workflow_HMR_XTK::perform() problem with mNumCriteria. "
                    "This can happen if the xtk interface interfaces different refinement level in the first optimization iteration" );

            Vector< real > tVector( mNumCriteria, std::numeric_limits< real >::quiet_NaN() );

            if ( tDeleteXTK )
            {
                // delete the xtk
                delete tXTKPerformer;
            }

            return tVector;
        }

        // store whether the new ghost has been used
        bool tUseNewGhostSets = tXTKPerformer->uses_SPG_based_enrichment();

        // XTK perform - enrich - ghost - multigrid
        tXTKPerformer->perform_enrichment();

        // set by the every-step in-place redistance below: the design was replaced by its SDF after this cut
        bool tRecut = false;

        // Stage *: Interface-based SDF reinitialization (after XTK has interface facets)
        if ( mPerformerManager->mReinitializePerformer.size() > 0 )
        {
            Reinitialize_Performer* tReinitPerformer = mPerformerManager->mReinitializePerformer( 0 ).get();

            if ( tReinitPerformer->get_reinit_mode() == ReinitMode::INTERFACE_SDF )
            {
                sint tReinitFreq = tReinitPerformer->get_reinitialization_frequency();

                // One-time SDF initialization at the first iteration: start the optimizer
                // from a true signed distance field rather than the raw initial level set.
                bool tInitAtStart = tReinitPerformer->get_initialize_at_start()
                                 && !mSDFInitialized && tOptIter > 0;

                // Regular scheduled reinit on the configured frequency
                bool tScheduled = tOptIter > 0 && tOptIter % tReinitFreq == 0;

                // In-place mode: redistance the field in-place (no restart) every tReinitFreq
                // iterations (reinitialization_frequency; =1 means every step). Keeps the field near
                // a true SDF on that cadence; smoothing fires on the same schedule.
                bool tEveryStep  = tReinitPerformer->get_reinit_every_step();
                bool tRedistance = tEveryStep ? ( tOptIter > 0 && tOptIter % tReinitFreq == 0 ) : tScheduled;
                bool tDoSmooth   = tEveryStep ? ( tOptIter % tReinitFreq == 0 ) : true;

                // The native MMA in reinit_on_accept_only mode suppresses the in-place redistance for its inner
                // conservatism trials (they must be evaluated at the design the approximation was built for). The flag is
                // set on rank 0 only; broadcast it so that every rank takes the same branch.
                uint tSuppressMsg = 0;
                if ( tEveryStep )
                {
                    tSuppressMsg = this->is_reinit_suppressed() ? 1 : 0;
                    broadcast( tSuppressMsg );

                    // the taper decision is made on rank 0 (notify_accepted_iterate); share it
                    uint tWeanedMsg = mSDFWeanedOff ? 1 : 0;
                    broadcast( tWeanedMsg );
                    mSDFWeanedOff = tWeanedMsg == 1;
                }
                bool tSuppressed = tSuppressMsg == 1;

                // Skip SDF reinit on restart pass (first iteration after restart)
                if ( mSkipNextSDFReinit )
                {
                    MORIS_LOG_INFO( "Skipping SDF reinit on restart pass" );
                    mSkipNextSDFReinit = false;
                }
                else if ( tSuppressed && ( tInitAtStart || tRedistance ) )
                {
                    MORIS_LOG_INFO( "SDF redistance suppressed for this evaluation (inner GCMMA trial, reinit_on_accept_only)" );
                }
                // Trigger SDF reinit: the one-time start-of-run init, every-step redistance, or schedule
                else if ( ( tInitAtStart || tRedistance ) && !mSDFReinitPending && !mSDFWeanedOff )
                {
                    if ( tInitAtStart )
                    {
                        MORIS_LOG_INFO( "Initializing design as SDF on first iteration (sdf_initialize_at_start)" );
                    }
                    Tracer tTracer( "GEN", "Levelset", "Interface-SDF-Reinitialize" );

                    // the redistance writes B-spline coefficients into the ADV vector, but with an ADV filter the ADVs are the
                    // unfiltered variables psi; writing the redistanced field into psi would filter it again at every
                    // redistance, which acts as a curvature flow and moves the interface
                    MORIS_ERROR( not mPerformerManager->mGENPerformer( 0 )->has_adv_filter(),
                            "Workflow_HMR_XTK: in-place SDF reinitialization is not supported with adv_filter_radius > 0; "
                            "turn the redistance off (reinitialization_frequency) or the filter off" );

                    // Get XTK cut integration mesh
                    xtk::Cut_Integration_Mesh* tCutMeshPtr = tXTKPerformer->get_cut_integration_mesh();

                    MORIS_ERROR( tCutMeshPtr != nullptr,
                            "Workflow_HMR_XTK::perform - Cut integration mesh not available for SDF reinit" );

                    xtk::Cut_Integration_Mesh& tCutMesh = *tCutMeshPtr;

                    // Get interpolation mesh
                    mtk::Mesh* tInterpMesh = mPerformerManager->mMTKPerformer( 0 )->get_interpolation_mesh( 0 );

                    // Compute SDF from interface and update field
                    std::shared_ptr< mtk::Field > tSDFField;
                    tReinitPerformer->compute_sdf_from_interface(
                            tCutMesh,
                            tInterpMesh,
                            mPerformerManager->mGENPerformer,
                            mPerformerManager->mMTKPerformer,
                            tSDFField,
                            tDoSmooth );

                    // Every-step in-place redistance: overwrite the ADVs with the redistanced
                    // coefficients and continue WITHOUT a NaN/restart. The native MMA optimizer
                    // (cl_OPT_Algorithm_MMA) owns its loop and writes these back into its design
                    // vector, so the redistanced design persists and the criteria/gradient stay
                    // consistent. (The TPL GCMMA cannot do this -- it would hit the gradient
                    // consistency assert -- so every-step mode requires algorithm = "mma".)
                    // Skipped on iterations where the basis size changed (remesh handled elsewhere).
                    if ( tSDFField != nullptr && tEveryStep )
                    {
                        // KNOWN LIMITATION (2026-06): this every-step in-place apply modifies the
                        // design vector BETWEEN the criteria and gradient evaluations, which trips
                        // Problem::compute_design_criteria_gradients' consistency assert
                        // (mADVNormTolerance) after a few hundred iterations. The slice-write below
                        // is correct (reinit applies, PDVs preserved, converges near baseline) but
                        // the in-place/no-restart timing is not. The consistency-safe path is
                        // scheduled-restart mode (every_step=false): apply the SAME slice-write in
                        // Workflow_HMR_XTK::initialize() so the restart recomputes criteria+gradients
                        // at the redistanced design. See runs/studies/sdf_reinit_sweep/FINDINGS.md.
                        //
                        // Write the redistanced level-set coefficients into ONLY that geometry's
                        // slice of the global ADV vector, identified by its determining ADV ids.
                        // This preserves every other design's ADVs and ALL PDVs (no FEM-coupled
                        // DENSITY-PDV desync -- unlike distribute_advs, which rebuilds GEN and
                        // segfaults the PDV host on the no-restart path), and removes the broken
                        // requirement that the single-field coefficient count equal the whole
                        // composite ADV vector. The next optimizer iteration's normal set_advs()
                        // re-syncs all fields and PDVs consistently.
                        // a copy: the interface-preserving option below edits it
                        Matrix< DDRMat > tCoeffs = tReinitPerformer->get_coefficients();
                        // positions in aNewADVs, not global ADV IDs (MB-017: using the IDs shifted the design by one index)
                        Vector< sint >          tLSAdvIds = mPerformerManager->mGENPerformer( 0 )
                                ->get_design_adv_positions( tReinitPerformer->get_adv_field_name() );

                        MORIS_ERROR( tLSAdvIds.size() > 0,
                                "SDF in-place redistance: GEN has no design named '%s' (the mapping sublist's adv_field), so the "
                                "redistance cannot be applied (MB-019: this used to be skipped silently every iteration).",
                                tReinitPerformer->get_adv_field_name().c_str() );

                        if ( tLSAdvIds.size() == tCoeffs.n_rows() && tCoeffs.n_rows() > 0 )
                        {
                            if ( tReinitPerformer->get_diagnostics() || tReinitPerformer->get_preserve_interface() )
                            {
                                // psi: the design coefficients this evaluation was called with, gathered from the same slice
                                Matrix< DDRMat > tOldCoeffs( tCoeffs.n_rows(), 1, 0.0 );
                                for ( uint iC = 0; iC < tCoeffs.n_rows(); ++iC )
                                {
                                    sint tAdvId = tLSAdvIds( iC );
                                    if ( tAdvId >= 0 && tAdvId < (sint)aNewADVs.size() )
                                    {
                                        tOldCoeffs( iC ) = aNewADVs( tAdvId );
                                    }
                                }

                                // sdf_preserve_interface (MB-027): keep every coefficient that can reach a cut cell, so the
                                // cut, the criteria and their sensitivities are those of psi bit for bit
                                if ( tReinitPerformer->get_preserve_interface() )
                                {
                                    tReinitPerformer->preserve_interface( mPerformerManager->mGENPerformer, tOldCoeffs, tCoeffs );
                                }

                                if ( tReinitPerformer->get_diagnostics() )
                                {
                                    tReinitPerformer->log_redistance_footprint( mPerformerManager->mGENPerformer, tOldCoeffs, tCoeffs );
                                }
                            }

                            uint tWritten = 0;
                            for ( uint iC = 0; iC < tCoeffs.n_rows(); ++iC )
                            {
                                sint tAdvId = tLSAdvIds( iC );
                                if ( tAdvId >= 0 && tAdvId < (sint)aNewADVs.size() )
                                {
                                    aNewADVs( tAdvId ) = tCoeffs( iC );
                                    ++tWritten;
                                }
                            }
                            mSDFInitialized = true;
                            tRecut          = tWritten > 0;
                            MORIS_LOG_INFO( "SDF redistanced in-place every step (no restart): wrote %u of %u level-set coeffs into the %u-ADV vector slice%s",
                                    tWritten, (uint)tCoeffs.n_rows(), (uint)aNewADVs.size(), tDoSmooth ? ", smoothed" : "" );
                        }
                        else
                        {
                            MORIS_LOG_INFO( "SDF redistance skipped this step: %u level-set ADV ids vs %u coeffs (mesh changed)",
                                    (uint)tLSAdvIds.size(), (uint)tCoeffs.n_rows() );
                        }
                        // no NaN, no restart -- fall through to FEM with the redistanced design
                    }
                    // Distribute new SDF to ADVs (restart-based reinit; scheduled mode every reinit_freq)
                    else if ( tSDFField != nullptr )
                    {
                        // Wean-off decision driven by the actual DESIGN-VARIABLE change between
                        // reinits, ||x_k - x_{k-1}|| / ||x_{k-1}||. Unlike the reinit's field drift
                        // (which is dominated by |grad phi| re-normalization and stays large), the
                        // design-variable change genuinely shrinks as GCMMA converges. When it falls
                        // below the tolerance the design has settled, so applying further reinits
                        // would only perturb a converged design: disable reinit for the rest of the
                        // run. Skipped for the one-time start init and the first scheduled reinit
                        // (no previous design to compare against yet).
                        real tWeanTol = tReinitPerformer->get_wean_tol();
                        bool tWeanOff = false;
                        if ( !tInitAtStart && tWeanTol > 0.0
                                && mPrevReinitADVs.size() == aNewADVs.size() && aNewADVs.size() > 0 )
                        {
                            real tNumer = 0.0;
                            real tDenom = 0.0;
                            for ( uint iADV = 0; iADV < aNewADVs.size(); ++iADV )
                            {
                                tNumer += std::pow( aNewADVs( iADV ) - mPrevReinitADVs( iADV ), 2 );
                                tDenom += std::pow( mPrevReinitADVs( iADV ), 2 );
                            }
                            real tDVChange = ( tDenom > 0.0 ) ? std::sqrt( tNumer / tDenom ) : 0.0;
                            tWeanOff = ( tDVChange < tWeanTol );
                            MORIS_LOG_INFO( "SDF reinit design-variable change: %.4e (wean tol %.4e)%s",
                                    tDVChange, tWeanTol, tWeanOff ? " -> wean off" : "" );
                        }

                        // Record the current design as the reference for the next reinit
                        mPrevReinitADVs = aNewADVs;

                        if ( !tInitAtStart && tWeanOff )
                        {
                            MORIS_LOG_INFO( "SDF reinit weaned off: design-variable change below tolerance, no further reinits" );
                            mSDFWeanedOff = true;
                        }
                        else
                        {
                            // The design is now a true SDF; the one-time start-of-run
                            // initialization is satisfied and will not fire again.
                            mSDFInitialized = true;

                            // Mark that SDF coefficients are pending for NEXT iteration
                            // They will be applied at the start of the next get_criteria call
                            // (at lines 247-276) BEFORE XTK runs
                            mSDFReinitPending = true;

                            MORIS_LOG_INFO( "SDF reinitialization computed, returning NANs to trigger restart" );

                            // Return NaNs to trigger optimizer restart
                            // On restart, ADVs will be applied BEFORE XTK runs
                            mInitializeOptimizationRestart = true;

                            // Build and return vector of NANs
                            Vector< real > tNaNVector( mNumCriteria, std::numeric_limits< real >::quiet_NaN() );
                            return tNaNVector;  // Exit early - don't continue to FEM
                        }
                    }
                }
            }
        }

        // Every-step in-place redistance: the cut above was made on the pre-redistance design psi, but the optimizer
        // records phi = R(psi). Re-cut on phi so that FEM, the criteria, the adjoint and the sensitivities all belong to
        // the design the optimizer holds (otherwise each step uses psi's gradients for phi: D-015, the MBB detached from
        // its load pad at eval 13). Mirrors the per-iteration GEN -> XTK sequence above; the psi cut has only been used
        // by the redistance itself, nothing else holds it yet.
        if ( tRecut )
        {
            Tracer tRecutTracer( "WRK", "HMR-XTK Workflow", "Re-cut on redistanced design" );

            delete tXTKPerformer;

            mPerformerManager->mGENPerformer( 0 )->set_advs( aNewADVs );

            tXTKPerformer = new xtk::Model( tXTKParameterList( 0 )( 0 ) );
            tMTKPerformer = std::make_shared< mtk::Mesh_Manager >();

            tXTKPerformer->set_geometry_engine( mPerformerManager->mGENPerformer( 0 ).get() );
            tXTKPerformer->set_input_performer( mPerformerManager->mMTKPerformer( 0 ) );
            tXTKPerformer->set_output_performer( tMTKPerformer );

            mPerformerManager->mGENPerformer( 0 )->reset_mesh_information(
                    mPerformerManager->mMTKPerformer( 0 )->get_interpolation_mesh( 0 ) );

            tDeleteXTK = tXTKPerformer->delete_xtk_after_generation();

            if ( not tXTKPerformer->perform_decomposition() )
            {
                mInitializeOptimizationRestart = true;

                Vector< real > tVector( mNumCriteria, std::numeric_limits< real >::quiet_NaN() );

                if ( tDeleteXTK )
                {
                    delete tXTKPerformer;
                }

                return tVector;
            }

            tUseNewGhostSets = tXTKPerformer->uses_SPG_based_enrichment();

            tXTKPerformer->perform_enrichment();

            MORIS_LOG_INFO( "Re-cut on the redistanced design before FEM" );
        }

        // NOTE: mSDFReinitPending is NOT cleared here.
        // It persists until the next iteration, where it's checked at lines 247-276
        // and the ADVs are applied BEFORE XTK runs (then the flag is cleared there)

        if ( tDeleteXTK )
        {
            // construct the data base with the mtk performer from xtk
            mPerformerManager->mDataBasePerformer( 0 ) = std::make_shared< DataBase_Performer >( tMTKPerformer );

            // create the mtk performer that will hold the data base mesh pair and set it
            std::shared_ptr< mtk::Mesh_Manager > tMTKDataBasePerformer = std::make_shared< mtk::Mesh_Manager >();
            mPerformerManager->mDataBasePerformer( 0 )->set_output_performer( tMTKDataBasePerformer );

            // turn off the mesh check if no FEM-model will be constructed on the mesh
            mPerformerManager->mDataBasePerformer( 0 )->set_mesh_check( !tXTKPerformer->kill_workflow_flag() );

            // perform the mtk data base
            mPerformerManager->mDataBasePerformer( 0 )->perform();

            // set the mtk performer
            mPerformerManager->mMTKPerformer( 1 ) = tMTKDataBasePerformer;
        }

        // stop workflow if only pre-processing output is requested
        if ( tXTKPerformer->only_generate_xtk_temp() )
        {
            MORIS_LOG( "------------------------------------------------------------------------------" );
            MORIS_LOG( "Only output of the foreground mesh requested. Stopping workflow after XTK/MTK." );
            MORIS_LOG( "------------------------------------------------------------------------------" );

            // delete data that needs to be deleted explicitly to prevent memory leaks
            delete tXTKPerformer;
            mPerformerManager->mDataBasePerformer( 0 )->free_memory();

            // return empty vector for design criteria
            Vector< real > tVector( 1, std::numeric_limits< real >::quiet_NaN() );
            return tVector;
        }

        // output T-matrices and/or MPCs if requested
        this->output_T_matrices( tMTKPerformer, tXTKPerformer );

        // stop workflow if T-Matrices have been outputted
        if ( tXTKPerformer->kill_workflow_flag() )
        {
            MORIS_LOG( "----------------------------------------------------------------------------------------------------" );
            MORIS_LOG( "T-Matrix output or triangulation of all elements in post requested. Stopping workflow after XTK/MTK." );
            MORIS_LOG( "----------------------------------------------------------------------------------------------------" );

            // delete data that needs to be deleted explicitly to prevent memory leaks
            delete tXTKPerformer;
            mPerformerManager->mDataBasePerformer( 0 )->free_memory();

            // return empty vector for design criteria
            Vector< real > tVector( 1, std::numeric_limits< real >::quiet_NaN() );
            return tVector;
        }

        if ( tDeleteXTK )
        {
            // delete the xtk-performer
            delete tXTKPerformer;
        }
        else
        {
            // IMPORTANT!!! do not overwrite previous XTK and MTK performer before we know if this XTK performer triggers a restart.
            // otherwise the fem::field meshes are deleted and cannot be used anymore.
            mPerformerManager->mXTKPerformer( 0 ) = std::shared_ptr< xtk::Model >( tXTKPerformer );
            mPerformerManager->mMTKPerformer( 1 ) = tMTKPerformer;
        }

        // mtk::Mesh_Checker tMeshCheckerXTK(
        //         0,
        //         mPerformerManager->mMTKPerformer( 1 )->get_mesh_pair(0).get_interpolation_mesh(),
        //         mPerformerManager->mMTKPerformer( 1 )->get_mesh_pair(0).get_integration_mesh() );
        // tMeshCheckerXTK.perform();
        // tMeshCheckerXTK.print_diagnostics();

        // mPerformerManager->mMTKPerformer( 1 )->get_mesh_pair(0).get_integration_mesh()->save_MPC_to_hdf5();
        // mPerformerManager->mMTKPerformer( 1 )->get_mesh_pair(0).get_integration_mesh()->save_IG_global_T_matrix_to_file();

        // get the MIG parameter list
        Module_Parameter_Lists tMIGParameterList = mPerformerManager->mLibrary->get_parameters_for_module( Module_Type::MIG );

        // check if there are MIG parameters specified
        if ( tMIGParameterList.size() > 0 )
        {
            moris::mig::MIG tMIGPerformer = moris::mig::MIG(
                    mPerformerManager->mMTKPerformer( 1 ),
                    tMIGParameterList( 0 )( 0 ),
                    mPerformerManager->mGENPerformer( 0 ).get() );

            tMIGPerformer.perform();
        }

        if ( tDeleteXTK )
        {
            // free the memory and delete the unused data
            mPerformerManager->mDataBasePerformer( 0 )->free_memory();
        }

        mPerformerManager->mMDLPerformer( 0 )->set_performer( mPerformerManager->mMTKPerformer( 1 ) );

        // Assign PDVs
        mPerformerManager->mGENPerformer( 0 )->create_pdvs( mPerformerManager->mMTKPerformer( 1 )->get_mesh_pair( 0 ) );

        // Stage 3: MDL perform ---------------------------------------------------------------------

        mPerformerManager->mMDLPerformer( 0 )->set_design_variable_interface(
                mPerformerManager->mGENPerformer( 0 )->get_design_variable_interface() );

        // store whether the new ghost sets are being used
        mPerformerManager->mMDLPerformer( 0 )->set_use_new_ghost_mesh_sets( tUseNewGhostSets );

        mPerformerManager->mMDLPerformer( 0 )->initialize();

        mPerformerManager->mGENPerformer( 0 )->communicate_requested_IQIs();

        // Build MDL components and solve
        mPerformerManager->mMDLPerformer( 0 )->perform();
        mMDLFieldsOnCurrentMesh = true;

        // perform mapping at this stage between solution field and adv field as some data will be deleted
        // This is FIELD_REMAP mode only - INTERFACE_SDF mode handles reinitialization earlier (after XTK)
        if ( mPerformerManager->mReinitializePerformer.size() > 0 )
        {
            Reinitialize_Performer* tReinitPerformer = mPerformerManager->mReinitializePerformer( 0 ).get();
            
            // Only run perform() for FIELD_REMAP mode
            if ( tReinitPerformer->get_reinit_mode() == ReinitMode::FIELD_REMAP )
            {
                // decide if we need to perform mapping
                if ( tOptIter % tReinitPerformer->get_reinitialization_frequency() ==
                        tReinitPerformer->get_reinitialization_frequency() - 1 )
                {
                    // perform mapping of the solution to the adv
                    mPerformerManager->mReinitializePerformer( 0 )->perform( mPerformerManager->mHMRPerformer,
                            mPerformerManager->mGENPerformer,
                            mPerformerManager->mMTKPerformer,
                            mPerformerManager->mMDLPerformer );
                }
            }
        }

        // evaluate IQIs
        Vector< moris::Matrix< DDRMat > > tVal = mPerformerManager->mMDLPerformer( 0 )->get_IQI_values();

        // get number of design criteria
        mNumCriteria = tVal.size();

        // Communicate IQIs
        for ( uint iIQIIndex = 0; iIQIIndex < mNumCriteria; iIQIIndex++ )
        {
            tVal( iIQIIndex )( 0 ) = sum_all( tVal( iIQIIndex )( 0 ) );
        }

        // build vector of design criteria
        Vector< real > tVector( mNumCriteria, 0.0 );
        for ( uint iCriteriaIndex = 0; iCriteriaIndex < mNumCriteria; iCriteriaIndex++ )
        {
            tVector( iCriteriaIndex ) = tVal( iCriteriaIndex )( 0 );
        }

        // increment optimization iteration counter
        mIter++;

        // return vector of design criteria
        return tVector;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Workflow_HMR_XTK::notify_accepted_iterate( const Vector< real >& aPrevious, const Vector< real >& aAccepted )
    {
        if ( mPerformerManager->mReinitializePerformer.size() == 0 || mSDFWeanedOff )
        {
            return;
        }
        Reinitialize_Performer* tReinitPerformer = mPerformerManager->mReinitializePerformer( 0 ).get();
        if ( tReinitPerformer->get_reinit_mode() != ReinitMode::INTERFACE_SDF || !tReinitPerformer->get_reinit_every_step() )
        {
            return;
        }
        if ( !mWeanMonitor )
        {
            mWeanMonitor = std::make_unique< SDF_Wean_Monitor >( tReinitPerformer->create_wean_monitor() );
        }
        if ( mWeanMonitor->mTol <= 0.0 && mWeanMonitor->mStopIter == 0 )
        {
            return;    // taper not configured
        }

        // an empty accepted vector: accepted at the inner-loop cap, not a convergence sample
        if ( aAccepted.size() == 0 )
        {
            bool tWeanedCapped = mWeanMonitor->update_capped();
            MORIS_LOG_INFO( "SDF taper: accepted iterate %u at the inner-loop cap (not counted as a small change)%s",
                    mWeanMonitor->mAccepted, tWeanedCapped ? " -> weaned off" : "" );
            if ( tWeanedCapped )
            {
                MORIS_LOG_INFO( "SDF in-place redistance weaned off: %s; no further redistances", mWeanMonitor->mReason.c_str() );
                mSDFWeanedOff = true;
            }
            return;
        }

        real tNumer = 0.0;
        real tDenom = 0.0;
        for ( uint iADV = 0; iADV < aPrevious.size() && iADV < aAccepted.size(); iADV++ )
        {
            tNumer += std::pow( aAccepted( iADV ) - aPrevious( iADV ), 2 );
            tDenom += std::pow( aPrevious( iADV ), 2 );
        }
        real tChange = tDenom > 0.0 ? std::sqrt( tNumer / tDenom ) : 0.0;

        bool tWeaned = mWeanMonitor->update( tChange );
        MORIS_LOG_INFO( "SDF taper: accepted iterate %u, relative design change %.4e (wean tol %.4e, %u/%u below, stop iter %u)%s",
                mWeanMonitor->mAccepted, tChange, mWeanMonitor->mTol, mWeanMonitor->mBelow, mWeanMonitor->mPatience,
                mWeanMonitor->mStopIter, tWeaned ? " -> weaned off" : "" );
        if ( tWeaned )
        {
            MORIS_LOG_INFO( "SDF in-place redistance weaned off: %s; no further redistances", mWeanMonitor->mReason.c_str() );
            mSDFWeanedOff = true;
        }
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat >
    Workflow_HMR_XTK::compute_dcriteria_dadv()
    {
        mPerformerManager->mGENPerformer( 0 )->communicate_requested_IQIs();

        mPerformerManager->mMDLPerformer( 0 )->perform( 1 );

        Matrix< DDRMat > tDCriteriaDAdv = mPerformerManager->mGENPerformer( 0 )->get_dcriteria_dadv();

        if ( par_rank() == 0 )
        {
            MORIS_LOG_INFO( "--------------------------------------------------------------------------------" );
            MORIS_LOG_INFO( "Gradients of design criteria wrt ADVs:" );

            for ( uint i = 0; i < tDCriteriaDAdv.n_rows(); ++i )
            {
                Matrix< DDRMat > tdIQIdAdv = tDCriteriaDAdv.get_row( i );

                auto tItrMin = std::min_element( tdIQIdAdv.data(), tdIQIdAdv.data() + tdIQIdAdv.numel() );
                auto tIndMin = std::distance( tdIQIdAdv.data(), tItrMin );

                auto tItrMax = std::max_element( tdIQIdAdv.data(), tdIQIdAdv.data() + tdIQIdAdv.numel() );
                auto tIndMax = std::distance( tdIQIdAdv.data(), tItrMax );

                MORIS_LOG_INFO( "Criteria(%i): norm = %e   min = %e  (index = %li)   max = %e  (index = %li)",
                        i,
                        norm( tdIQIdAdv ),
                        tdIQIdAdv.min(),
                        tIndMin,
                        tdIQIdAdv.max(),
                        tIndMax );
            }

            MORIS_LOG_INFO( "--------------------------------------------------------------------------------" );
        }

        mPerformerManager->mMDLPerformer( 0 )->free_memory();

        return tDCriteriaDAdv;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Workflow_HMR_XTK::output_T_matrices(
            const std::shared_ptr< mtk::Mesh_Manager >& aMTKPerformer,
            xtk::Model* const &                         aXTKPerformer )
    {
        // Output T-matrices if requested
        std::string tTmatrixFileName = aXTKPerformer->get_global_T_matrix_output_file_name();
        if ( tTmatrixFileName != "" )
        {
            mPerformerManager->mMTKPerformer( 1 )->get_mesh_pair( 0 ).get_integration_mesh()->save_IG_global_T_matrix_to_file( tTmatrixFileName );

            // return flag stopping the workflow after the T-Matrix output
            return true;
        }

        // Output T-matrices if requested
        std::string tElementalTmatrixFileName = aXTKPerformer->get_elemental_T_matrix_output_file_name();
        if ( tElementalTmatrixFileName != "" )
        {
            uint tNumBsplineMeshes = mPerformerManager->mMTKPerformer( 1 )->get_mesh_pair( 0 ).get_interpolation_mesh()->get_num_interpolations();
            mPerformerManager->mMTKPerformer( 1 )->get_mesh_pair( 0 ).get_integration_mesh()->save_elemental_T_matrices_to_file( tElementalTmatrixFileName, tNumBsplineMeshes );

            // return flag stopping the workflow after the T-Matrix output
            return true;
        }

        // Output MPCs if requested
        std::string tMpcFileName = aXTKPerformer->get_MPC_output_file_name();
        if ( tMpcFileName != "" )
        {
            mPerformerManager->mMTKPerformer( 1 )->get_mesh_pair( 0 ).get_integration_mesh()->save_MPC_to_hdf5( tMpcFileName );
        }

        return false;
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::wrk
