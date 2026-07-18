/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_IOS_Input_Deck_Typed_Layer.cpp
 *
 * Tests for the Deck API v2 typed layer: the mesh-set vocabulary (pinned against
 * set-name strings from existing EXA decks) and the criteria-expression machinery
 * (values and reverse-mode gradients against hand-computed results).
 */

#include <catch.hpp>

#include <map>

#include "cl_Input_Deck_Vocabulary.hpp"
#include "cl_Input_Deck_Expressions.hpp"
#include "cl_Vector.hpp"

namespace moris
{
    //------------------------------------------------------------------------------------------------------------------

    TEST_CASE( "deck vocabulary generates the XTK set-name grammar", "[IOS],[deck_semantics],[vocabulary]" )
    {
        deck::Phase tVoid( 0 );
        deck::Phase tMaterial( 1 );
        deck::Phase tFrame( 2 );

        // strings pinned against hand-written EXA decks
        CHECK( tMaterial.bulk() == "HMR_dummy_n_p1,HMR_dummy_c_p1" );            // Shape_Sensitivity_Sweep
        CHECK( tMaterial.side( deck::Side::Left ) == "SideSet_4_n_p1,SideSet_4_c_p1" );
        CHECK( tMaterial.side( deck::Side::Right ) == "SideSet_2_n_p1,SideSet_2_c_p1" );
        CHECK( tMaterial.ghost() == "ghost_p1" );
        CHECK( deck::interface( tMaterial, tVoid ) == "iside_b0_1_b1_0" );       // boxbeam material/void interface
        CHECK( deck::between( tMaterial, tFrame ) == "dbl_iside_p0_1_p1_2" );    // boxbeam interior/frame Nitsche
        CHECK( deck::join( tFrame.bulk(), tMaterial.bulk() )
                == "HMR_dummy_n_p2,HMR_dummy_c_p2,HMR_dummy_n_p1,HMR_dummy_c_p1" );

        CHECK( deck::Dofs::Displacement2D == "UX,UY" );
        CHECK( deck::Dofs::Displacement3D == "UX,UY,UZ" );
    }

    //------------------------------------------------------------------------------------------------------------------

    TEST_CASE( "deck expressions: values and reverse-mode gradients", "[IOS],[deck_semantics],[expressions]" )
    {
        // criteria layout: c0 = strain energy, c1 = volume, c2 = perimeter
        std::map< std::string, uint > tIndexOf = { { "SE", 0 }, { "Vol", 1 }, { "Per", 2 } };
        Vector< real >                tCriteria = { 3.0, 0.5, 8.0 };

        SECTION( "pass-through objective (the Shape_Sensitivity_Sweep form)" )
        {
            deck::Expression tObjective = deck::criterion( "SE" );

            CHECK( tObjective.evaluate( tCriteria, tIndexOf ) == 3.0 );

            Vector< real > tGradient( 3, 0.0 );
            tObjective.accumulate_gradient( tCriteria, tIndexOf, tGradient );
            CHECK( tGradient( 0 ) == 1.0 );
            CHECK( tGradient( 1 ) == 0.0 );
            CHECK( tGradient( 2 ) == 0.0 );
        }

        SECTION( "weighted normalized sum (the box-beam objective form)" )
        {
            const real tSE0 = 2.0;
            const real tP0  = 4.0;
            const real tW   = 0.2;

            deck::Expression tObjective =
                    deck::criterion( "SE" ) / tSE0 + tW * deck::criterion( "Per" ) / tP0;

            CHECK( tObjective.evaluate( tCriteria, tIndexOf ) == Approx( 3.0 / 2.0 + 0.2 * 8.0 / 4.0 ) );

            Vector< real > tGradient( 3, 0.0 );
            tObjective.accumulate_gradient( tCriteria, tIndexOf, tGradient );
            CHECK( tGradient( 0 ) == Approx( 1.0 / tSE0 ) );
            CHECK( tGradient( 1 ) == 0.0 );
            CHECK( tGradient( 2 ) == Approx( tW / tP0 ) );
        }

        SECTION( "product, quotient, and negation rules" )
        {
            // f = -( SE * Vol ) / Per : df/dSE = -Vol/Per, df/dVol = -SE/Per, df/dPer = SE*Vol/Per^2
            deck::Expression tExpression = -( deck::criterion( "SE" ) * deck::criterion( "Vol" ) ) / deck::criterion( "Per" );

            CHECK( tExpression.evaluate( tCriteria, tIndexOf ) == Approx( -3.0 * 0.5 / 8.0 ) );

            Vector< real > tGradient( 3, 0.0 );
            tExpression.accumulate_gradient( tCriteria, tIndexOf, tGradient );
            CHECK( tGradient( 0 ) == Approx( -0.5 / 8.0 ) );
            CHECK( tGradient( 1 ) == Approx( -3.0 / 8.0 ) );
            CHECK( tGradient( 2 ) == Approx( 3.0 * 0.5 / 64.0 ) );
        }

        SECTION( "constraint sugar sets type and shifts the bound" )
        {
            deck::Constraint tInequality = deck::criterion( "Vol" ) / 0.4 - 1.0 <= 0.0;
            CHECK( tInequality.mType == 1 );
            CHECK( tInequality.mExpression.evaluate( tCriteria, tIndexOf ) == Approx( 0.5 / 0.4 - 1.0 ) );

            deck::Constraint tEquality = deck::criterion( "Vol" ) == 0.5;
            CHECK( tEquality.mType == 0 );
            CHECK( tEquality.mExpression.evaluate( tCriteria, tIndexOf ) == Approx( 0.0 ).margin( 1e-14 ) );
        }

        SECTION( "criteria collection is first-appearance ordered and deduplicated" )
        {
            deck::Expression tExpression =
                    deck::criterion( "Per" ) + deck::criterion( "SE" ) * deck::criterion( "Per" );

            Vector< std::string > tOrder;
            tExpression.collect_criteria( tOrder );

            REQUIRE( tOrder.size() == 2 );
            CHECK( tOrder( 0 ) == "Per" );
            CHECK( tOrder( 1 ) == "SE" );
        }
    }

    //------------------------------------------------------------------------------------------------------------------

}    // namespace moris
