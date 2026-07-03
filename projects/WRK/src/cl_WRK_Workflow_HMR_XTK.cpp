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
#include "cl_HMR.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Mesh_Checker.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_XTK_Model.hpp"
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

namespace moris::wrk
{

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
        // after first solve in optimization process
        else
        {
            // get access to the performers
            std::shared_ptr< gen::Geometry_Engine >         tGeometryEngine     = mPerformerManager->mGENPerformer( 0 );
            std::shared_ptr< mdl::Model >                   tModelPerformer     = mPerformerManager->mMDLPerformer( 0 );
            std::shared_ptr< Remeshing_Mini_Performer >     tRemeshingPerformer = mPerformerManager->mRemeshingMiniPerformer( 0 );
            Vector< std::shared_ptr< hmr::HMR > >&          tHmrPerformers      = mPerformerManager->mHMRPerformer;
            Vector< std::shared_ptr< mtk::Mesh_Manager > >& tMtkPerformers      = mPerformerManager->mMTKPerformer;

            // get refinement fields from GEN and MDL performers
            tFieldsIn.append( tGeometryEngine->get_mtk_fields() );
            tFieldsIn.append( tModelPerformer->get_mtk_fields() );

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

                // Skip SDF reinit on restart pass (first iteration after restart)
                if ( mSkipNextSDFReinit )
                {
                    MORIS_LOG_INFO( "Skipping SDF reinit on restart pass" );
                    mSkipNextSDFReinit = false;
                }
                // Trigger SDF reinit: the one-time start-of-run init, every-step redistance, or schedule
                else if ( ( tInitAtStart || tRedistance ) && !mSDFReinitPending && !mSDFWeanedOff )
                {
                    if ( tInitAtStart )
                    {
                        MORIS_LOG_INFO( "Initializing design as SDF on first iteration (sdf_initialize_at_start)" );
                    }
                    Tracer tTracer( "GEN", "Levelset", "Interface-SDF-Reinitialize" );

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
                        Matrix< DDRMat > const& tCoeffs = tReinitPerformer->get_coefficients();
                        if ( tCoeffs.n_rows() == aNewADVs.size() && aNewADVs.size() > 0 )
                        {
                            for ( uint iADV = 0; iADV < tCoeffs.n_rows(); ++iADV )
                            {
                                aNewADVs( iADV ) = tCoeffs( iADV );
                            }
                            mSDFInitialized = true;
                            MORIS_LOG_INFO( "SDF redistanced in-place every step (no restart): %u coeffs%s",
                                    (uint)tCoeffs.n_rows(), tDoSmooth ? ", smoothed" : "" );
                        }
                        else
                        {
                            MORIS_LOG_INFO( "SDF redistance skipped this step: coeff count %u != ADV count %u (mesh changed)",
                                    (uint)tCoeffs.n_rows(), (uint)aNewADVs.size() );
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
