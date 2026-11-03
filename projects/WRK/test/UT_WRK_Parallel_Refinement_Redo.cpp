/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_WRK_Parallel_Refinement_Redo.cpp
 *
 * MB-097: the redo of the GEN refinement on the current design (after a failed decomposition, MB-087, and on a sub-cell
 * feature re-detection, MB-089) and the restart paths of Workflow_HMR_XTK::initialize() (MB-092) are collective, so every
 * processor must take the same branch. The ADVs are complete on processor 0 only: the optimizer runs there and the other
 * processors evaluate through opt::Algorithm::dummy_solve() with an empty ADV vector. This test calls the workflow the
 * same way (the ADVs on processor 0, an empty vector elsewhere) on the WRK_Input_3.so deck; before the fix processor 0
 * entered the redo alone and the run hung at its first collective, or the other processors stopped with "remeshing
 * performer has not been built" on a restart. Run it under a timeout: a hang is the failure mode.
 *
 * The criteria are compared with the serial values (the 1-processor run of this test), so the redo on several
 * processors gives the serial mesh and solution.
 *
 */

#include "catch.hpp"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include "cl_Communication_Tools.hpp"
#include "paths.hpp"

#include "cl_WRK_Performer_Manager.hpp"
#include "cl_WRK_Workflow_HMR_XTK.hpp"
#include "cl_Library_Factory.hpp"

namespace moris::wrk
{
    namespace
    {
        // the ADVs as the optimizer hands them to the workflow: complete on processor 0, empty on the other processors
        Vector< real >
        advs_on_proc_0( const Vector< real >& aADVs )
        {
            return par_rank() == 0 ? aADVs : Vector< real >( 0 );
        }

        std::shared_ptr< Library_IO >
        load_input_3( bool aGenRefines )
        {
            setenv( "WRK_INPUT_3_GEN_REFINE", aGenRefines ? "1" : "0", 1 );

            std::string tInputFilePath = moris::get_moris_bin_dir() + "/lib/WRK_Input_3.so";

            Library_Factory               tLibraryFactory;
            std::shared_ptr< Library_IO > tLibrary = tLibraryFactory.create_Library( Library_Type::STANDARD );
            tLibrary->load_parameter_list( tInputFilePath, File_Type::SO_FILE );
            tLibrary->finalize();
            return tLibrary;
        }

        // one evaluation; the criteria (stored temperature of the solid, disk volume) must be finite on every processor
        Vector< real >
        evaluate( Workflow_HMR_XTK& aWorkflow, const Vector< real >& aADVs, const std::string& aLabel )
        {
            Vector< real > tADVs     = advs_on_proc_0( aADVs );
            Vector< real > tCriteria = aWorkflow.get_criteria( tADVs );

            REQUIRE( tCriteria.size() == 2 );
            CHECK( std::isfinite( tCriteria( 0 ) ) );
            CHECK( std::isfinite( tCriteria( 1 ) ) );

            if ( par_rank() == 0 )
            {
                std::cout << std::setprecision( 15 ) << "MB-097 " << aLabel << " on " << par_size() << " processor(s): criteria "
                          << tCriteria( 0 ) << " " << tCriteria( 1 ) << ", refinement redos " << aWorkflow.get_num_refinement_redos()
                          << '\n';
            }
            return tCriteria;
        }

        // the serial values (1-processor run of this test, ee7785d6a + the MB-097 fix)
        void
        check_serial_value( const Vector< real >& aCriteria, real aStoredSolid, real aVolumeDisk )
        {
            if ( par_rank() == 0 )
            {
                CHECK( aCriteria( 0 ) == Approx( aStoredSolid ).epsilon( 1.0e-8 ) );
                CHECK( aCriteria( 1 ) == Approx( aVolumeDisk ).epsilon( 1.0e-8 ) );
            }
        }
    }    // namespace

    TEST_CASE( "MB-097 refinement redo and restarts are collective", "[WRK],[MB-097]" )
    {
        Vector< real >  tADVs;
        Vector< real >  tLowerBounds;
        Vector< real >  tUpperBounds;
        Matrix< IdMat > tIjklIDs;

        // the GS-LS decks' path: GEN refines the mesh around the design
        SECTION( "GEN-refined mesh: redo after a failed decomposition, re-detection, restart" )
        {
            std::shared_ptr< Library_IO > tLibrary = load_input_3( true );
            Performer_Manager             tPerformerManager( tLibrary );
            Workflow_HMR_XTK              tWorkflow( &tPerformerManager );

            tWorkflow.initialize( tADVs, tLowerBounds, tUpperBounds, tIjklIDs );
            REQUIRE( tADVs.size() == ( par_rank() == 0 ? 3 : 0 ) );

            // 1: the design the mesh is refined for. The small disk crosses a coarse background edge twice: GEN refines
            // the two cells around it (MB-084), and the cut cells of the large disk to the same level (MB-085)
            Vector< real > tFirst = evaluate( tWorkflow, { 0.31, 0.18, 0.03 }, "initial design" );
            CHECK( tWorkflow.get_num_refinement_redos() == 0 );

            // 2: the large disk moves out of the refined band and the small disk vanishes: cut cells on levels 0 and 1,
            // the decomposition fails, the refinement is redone on this design (MB-087: nothing to refine, a uniform
            // mesh) and the cut succeeds in the same evaluation
            Vector< real > tMoved = evaluate( tWorkflow, { 0.66, 0.18, -0.05 }, "moved disk" );
            CHECK( tWorkflow.get_num_refinement_redos() == 1 );

            // 3: the small disk appears again: no cut cell sees it, the resolution check counts a hidden edge (more than
            // the threshold 0 and than the first cut after the last refinement), the refinement is redone on this design
            // (MB-089) and the disk is cut
            Vector< real > tFeature = evaluate( tWorkflow, { 0.66, 0.18, 0.03 }, "small disk" );
            CHECK( tWorkflow.get_num_refinement_redos() == 2 );

            // 4: the same design again: the mesh was refined for it, no redo
            Vector< real > tAgain = evaluate( tWorkflow, { 0.66, 0.18, 0.03 }, "small disk again" );
            CHECK( tWorkflow.get_num_refinement_redos() == 2 );

            // 5: an optimizer restart rebuilds HMR and GEN on the current design (refine_hmr_on_current_design())
            tWorkflow.initialize( tADVs, tLowerBounds, tUpperBounds, tIjklIDs );
            if ( par_rank() == 0 )
            {
                REQUIRE( tADVs.size() == 3 );
                CHECK( tADVs( 0 ) == 0.66 );
                CHECK( tADVs( 2 ) == 0.03 );
            }
            Vector< real > tRestarted = evaluate( tWorkflow, { 0.66, 0.18, 0.03 }, "after restart" );
            CHECK( tWorkflow.get_num_refinement_redos() == 2 );

            if ( par_rank() == 0 )
            {
                CHECK( tAgain( 0 ) == Approx( tFeature( 0 ) ).epsilon( 1.0e-12 ) );
                CHECK( tAgain( 1 ) == Approx( tFeature( 1 ) ).epsilon( 1.0e-12 ) );
                CHECK( tRestarted( 0 ) == Approx( tFeature( 0 ) ).epsilon( 1.0e-12 ) );
                CHECK( tRestarted( 1 ) == Approx( tFeature( 1 ) ).epsilon( 1.0e-12 ) );
            }

            check_serial_value( tFirst, 0.226679118762088, 0.103362722606113 );
            check_serial_value( tMoved, 0.272286316555731, 0.0977523559699494 );
            check_serial_value( tFeature, 0.266963529115395, 0.103345810335808 );
        }

        // a design-independent mesh: a restart keeps the mesh and rebuilds GEN on the current ADVs
        SECTION( "design-independent mesh: restart rebuilds GEN on the current ADVs" )
        {
            std::shared_ptr< Library_IO > tLibrary = load_input_3( false );
            Performer_Manager             tPerformerManager( tLibrary );
            Workflow_HMR_XTK              tWorkflow( &tPerformerManager );

            tWorkflow.initialize( tADVs, tLowerBounds, tUpperBounds, tIjklIDs );

            Vector< real > tFirst  = evaluate( tWorkflow, { 0.31, 0.18, -0.05 }, "fixed mesh, initial design" );
            Vector< real > tSecond = evaluate( tWorkflow, { 0.43, 0.21, -0.05 }, "fixed mesh, second design" );

            tWorkflow.initialize( tADVs, tLowerBounds, tUpperBounds, tIjklIDs );
            if ( par_rank() == 0 )
            {
                REQUIRE( tADVs.size() == 3 );
                CHECK( tADVs( 0 ) == 0.43 );
                CHECK( tADVs( 1 ) == 0.21 );
            }
            Vector< real > tRestarted = evaluate( tWorkflow, { 0.43, 0.21, -0.05 }, "fixed mesh, after restart" );
            CHECK( tWorkflow.get_num_refinement_redos() == 0 );

            if ( par_rank() == 0 )
            {
                CHECK( tRestarted( 0 ) == Approx( tSecond( 0 ) ).epsilon( 1.0e-12 ) );
                CHECK( tRestarted( 1 ) == Approx( tSecond( 1 ) ).epsilon( 1.0e-12 ) );
            }

            check_serial_value( tFirst, 0.232830772494089, 0.0978569069187233 );
            check_serial_value( tSecond, 0.226280854756853, 0.134165596907745 );
        }
    }
}    // namespace moris::wrk
