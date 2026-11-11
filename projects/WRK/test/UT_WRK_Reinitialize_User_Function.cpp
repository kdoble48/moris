/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_WRK_Reinitialize_User_Function.cpp
 *
 * sdf_reinit_mode = "user_function": a deck function replaces the ADVs of one design on the optimizer restart that the
 * OPT reinitialization schedule (reinitialize_interface_iter) requests.
 *  - the pure core Reinitialize_Performer::apply_user_reinitialization(): slice, clipping, failure handling;
 *  - an optimization (WRK_Input_2.so, native MMA): the function sees the last evaluated design, the next evaluation is
 *    exactly the function's result, and the run continues to its iteration limit. On several processors the function is
 *    called once, on processor 0, and every evaluated design equals the one of the 1-processor run.
 *
 */

#include "catch.hpp"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include "cl_Communication_Tools.hpp"
#include "paths.hpp"

#include "cl_WRK_Reinitialize_Performer.hpp"
#include "cl_WRK_Performer_Manager.hpp"
#include "cl_WRK_Workflow_Factory.hpp"
#include "cl_OPT_Manager.hpp"
#include "cl_OPT_Criteria_Interface.hpp"
#include "cl_Library_Factory.hpp"

namespace moris::wrk
{
    namespace
    {
        // a known transform of a design's ADVs: new_i = 2 old_i + iteration
        Vector< real >
        double_plus_iteration( const Vector< real >& aADVs, uint aIteration )
        {
            Vector< real > tNew( aADVs.size() );
            for ( uint i = 0; i < aADVs.size(); i++ )
            {
                tNew( i ) = 2.0 * aADVs( i ) + (real)aIteration;
            }
            return tNew;
        }

        Vector< real >
        returns_nothing( const Vector< real >& aADVs, uint aIteration )
        {
            return Vector< real >( 0 );
        }

        Vector< real >
        returns_nan( const Vector< real >& aADVs, uint aIteration )
        {
            Vector< real > tNew = aADVs;
            tNew( 0 )           = std::nan( "" );
            return tNew;
        }

        typedef uint ( *Count_Getter )();
        typedef uint ( *Indexed_Count_Getter )( uint );
        typedef Vector< real > ( *ADV_Getter )();
        typedef Vector< real > ( *Indexed_ADV_Getter )( uint );
    }    // namespace

    TEST_CASE( "user_function reinitialization: the pure core", "[WRK],[Reinitialize_User_Function]" )
    {
        // 5 ADVs; the design owns positions 3 and 1 (in its coefficient order)
        Vector< real > tADVs      = { 0.0, 1.0, 2.0, 3.0, 4.0 };
        Vector< sint > tPositions = { 3, 1 };
        Vector< real > tLower( 5, -10.0 );
        Vector< real > tUpper( 5, 10.0 );

        SECTION( "replaces the design's slice only" )
        {
            User_Reinit_Report tReport = Reinitialize_Performer::apply_user_reinitialization(
                    tADVs, tPositions, tLower, tUpper, &double_plus_iteration, 2, "double_plus_iteration" );

            // the function received ( 3, 1 ): new = ( 8, 4 ) written to positions 3 and 1
            CHECK( tADVs( 3 ) == 8.0 );
            CHECK( tADVs( 1 ) == 4.0 );
            CHECK( tADVs( 0 ) == 0.0 );
            CHECK( tADVs( 2 ) == 2.0 );
            CHECK( tADVs( 4 ) == 4.0 );
            CHECK( tReport.mNumADVs == 2 );
            CHECK( tReport.mNumClipped == 0 );
            CHECK( tReport.mMaxChange == Approx( 5.0 ) );
            CHECK( tReport.mRelChange == Approx( std::sqrt( 25.0 + 9.0 ) / std::sqrt( 9.0 + 1.0 ) ) );
        }

        SECTION( "clips each returned ADV to its own bounds" )
        {
            tUpper( 3 )                = 7.5;
            User_Reinit_Report tReport = Reinitialize_Performer::apply_user_reinitialization(
                    tADVs, tPositions, tLower, tUpper, &double_plus_iteration, 2, "double_plus_iteration" );
            CHECK( tADVs( 3 ) == 7.5 );
            CHECK( tADVs( 1 ) == 4.0 );
            CHECK( tReport.mNumClipped == 1 );
            CHECK( tReport.mMaxClip == Approx( 0.5 ) );
        }

        SECTION( "a failed function stops the run, never keeps the old ADVs silently" )
        {
            CHECK_THROWS( Reinitialize_Performer::apply_user_reinitialization(
                    tADVs, tPositions, tLower, tUpper, &returns_nothing, 2, "returns_nothing" ) );
            CHECK_THROWS( Reinitialize_Performer::apply_user_reinitialization(
                    tADVs, tPositions, tLower, tUpper, &returns_nan, 2, "returns_nan" ) );
            CHECK_THROWS( Reinitialize_Performer::apply_user_reinitialization(
                    tADVs, Vector< sint >{ 3, -1 }, tLower, tUpper, &double_plus_iteration, 2, "unresolved position" ) );
            CHECK_THROWS( Reinitialize_Performer::apply_user_reinitialization(
                    tADVs, Vector< sint >( 0 ), tLower, tUpper, &double_plus_iteration, 2, "empty design" ) );
        }
    }

    namespace
    {
        void
        check_optimization_with_user_reinitialization( bool aGenRefines, const Vector< Vector< real > >& aSerialDesigns )
        {
            setenv( "WRK_INPUT_2_GEN_REFINE", aGenRefines ? "1" : "0", 1 );

            std::string tInputFilePath = moris::get_moris_bin_dir() + "/lib/WRK_Input_2.so";

            Library_Factory               tLibraryFactory;
            std::shared_ptr< Library_IO > tLibrary = tLibraryFactory.create_Library( Library_Type::STANDARD );
            tLibrary->load_parameter_list( tInputFilePath, File_Type::SO_FILE );
            tLibrary->finalize();

            Module_Parameter_Lists tOPTParameterList = tLibrary->get_parameters_for_module( Module_Type::OPT );

            Performer_Manager tPerformerManager( tLibrary );

            Vector< std::shared_ptr< opt::Criteria_Interface > > tWorkflows = { create_workflow( "HMR_XTK", &tPerformerManager ) };

            typedef void ( *Reset_Function )();
            Reset_Function tReset = tLibrary->load_function< Reset_Function >( "reset_record" );
            tReset();

            opt::Manager tManager( tOPTParameterList, tWorkflows, tLibrary );
            REQUIRE_NOTHROW( tManager.perform() );

            // the optimizer, the objective and the deck function run on processor 0; the other processors only take part in
            // the evaluations and must not have called the function
            if ( par_rank() != 0 )
            {
                Count_Getter tNumCallsElsewhere = tLibrary->load_function< Count_Getter >( "get_num_reinit_calls" );
                CHECK( tNumCallsElsewhere() == 0 );
                return;
            }

            Count_Getter         tNumCalls       = tLibrary->load_function< Count_Getter >( "get_num_reinit_calls" );
            Count_Getter         tNumEvaluations = tLibrary->load_function< Count_Getter >( "get_num_evaluations" );
            Indexed_Count_Getter tCallAfter      = tLibrary->load_function< Indexed_Count_Getter >( "get_reinit_after_evaluation" );
            Indexed_ADV_Getter   tEvaluated      = tLibrary->load_function< Indexed_ADV_Getter >( "get_evaluated_advs" );
            ADV_Getter           tLastInput      = tLibrary->load_function< ADV_Getter >( "get_last_reinit_input" );

            // reinitialize_interface_iter = 3: the first call follows the third evaluation
            REQUIRE( tNumCalls() >= 1 );
            CHECK( tCallAfter( 0 ) == 3 );

            for ( uint iCall = 0; iCall < tNumCalls(); iCall++ )
            {
                uint tAfter = tCallAfter( iCall );
                REQUIRE( tAfter >= 1 );

                // the run continues after every call: the next evaluation is exactly the known transform (centre x 0.43,
                // radius 0.27), i.e. the optimizer restarted at the function's result
                REQUIRE( tNumEvaluations() > tAfter );
                Vector< real > tNext = tEvaluated( tAfter );
                REQUIRE( tNext.size() == 2 );
                CHECK( tNext( 0 ) == 0.43 );
                CHECK( tNext( 1 ) == 0.27 );

                // and differs from the design before the call (it was moved)
                Vector< real > tBefore = tEvaluated( tAfter - 1 );
                CHECK( std::abs( tBefore( 0 ) - 0.43 ) + std::abs( tBefore( 1 ) - 0.27 ) > 1.0e-6 );
            }

            // the function received the ADVs of the last evaluated design (not the unevaluated trial point)
            uint           tLastAfter = tCallAfter( tNumCalls() - 1 );
            Vector< real > tInput     = tLastInput();
            Vector< real > tLastEval  = tEvaluated( tLastAfter - 1 );
            REQUIRE( tInput.size() == 2 );
            CHECK( tInput( 0 ) == tLastEval( 0 ) );
            CHECK( tInput( 1 ) == tLastEval( 1 ) );

            // the optimizer kept iterating after the last call
            CHECK( tNumEvaluations() >= tLastAfter + 1 );

            // every evaluated design is the one of the 1-processor run (the serial values of this test)
            std::cout << std::setprecision( 15 ) << "user_function reinitialization, " << ( aGenRefines ? "GEN-refined" : "fixed" )
                      << " mesh, " << par_size() << " processor(s): " << tNumCalls() << " calls, " << tNumEvaluations()
                      << " evaluations\n";
            for ( uint iEvaluation = 0; iEvaluation < tNumEvaluations(); iEvaluation++ )
            {
                Vector< real > tDesign = tEvaluated( iEvaluation );
                std::cout << "  evaluation " << iEvaluation << ": " << tDesign( 0 ) << " " << tDesign( 1 ) << '\n';
            }
            REQUIRE( tNumEvaluations() == aSerialDesigns.size() );
            for ( uint iEvaluation = 0; iEvaluation < tNumEvaluations(); iEvaluation++ )
            {
                Vector< real > tDesign = tEvaluated( iEvaluation );
                CHECK( tDesign( 0 ) == Approx( aSerialDesigns( iEvaluation )( 0 ) ).epsilon( 1.0e-7 ) );
                CHECK( tDesign( 1 ) == Approx( aSerialDesigns( iEvaluation )( 1 ) ).epsilon( 1.0e-7 ) );
            }
        }
    }    // namespace

    TEST_CASE( "user_function reinitialization: an optimization continues from the function's result",
            "[WRK],[Reinitialize_User_Function]" )
    {
        // restart path 1: GEN refines the mesh, which is rebuilt on the new design (the GS-LS decks)
        SECTION( "GEN-refined mesh" )
        {
            check_optimization_with_user_reinitialization( true,
                    { { 0.5, 0.31 },
                            { 0.480333157041999, 0.308577965563353 },
                            { 0.460333243298349, 0.309410936861192 },
                            { 0.43, 0.27 },
                            { 0.410000094139711, 0.282499984838826 },
                            { 0.390000157920133, 0.294999969424 },
                            { 0.43, 0.27 } } );
        }
        // restart path 2: design-independent mesh, kept; GEN rebuilt at the new ADVs
        SECTION( "design-independent mesh" )
        {
            check_optimization_with_user_reinitialization( false,
                    { { 0.5, 0.31 },
                            { 0.480996891841292, 0.309900988164436 },
                            { 0.460996963946742, 0.310897962809957 },
                            { 0.43, 0.27 },
                            { 0.410000142370524, 0.282499983720318 },
                            { 0.390000220889688, 0.294999969030376 },
                            { 0.43, 0.27 } } );
        }
    }
}    // namespace moris::wrk
