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
#include "fn_Library_Interlink_Checks.hpp"
#include "parameters.hpp"
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

    TEST_CASE( "deck interlink findings", "[IOS],[deck_semantics],[interlinks]" )
    {
        // build a minimal in-memory deck: one property, one SP, one IWG referencing both
        Vector< Module_Parameter_Lists > tAll;
        for ( uint iModule = 0; iModule < (uint)Module_Type::END_ENUM; iModule++ )
        {
            tAll.push_back( Module_Parameter_Lists( (Module_Type)iModule ) );
        }

        Module_Parameter_Lists& tFem = tAll( (uint)Module_Type::FEM );

        tFem( FEM::PROPERTIES ).add_parameter_list();
        tFem.set( "property_name", "PropYoungs" );

        tFem( FEM::STABILIZATION ).add_parameter_list();
        tFem.set( "stabilization_name", "SPNitscheDirichlet" );
        tFem.set( "leader_properties", "PropYoungs,Material" );

        tFem( FEM::IWG ).add_parameter_list();
        tFem.set( "IWG_name", "IWGDirichlet" );
        tFem.set( "leader_properties", "PropYoungs,Dirichlet" );
        tFem.set( "stabilization_parameters", "SPNitscheDirichlet,DirichletNitsche" );

        tFem( FEM::IQI ).add_parameter_list();
        tFem.set( "IQI_name", "IQIBulkStrainEnergy" );

        SECTION( "a consistent deck has no findings" )
        {
            CHECK( collect_deck_interlink_findings( tAll ).size() == 0 );
        }

        SECTION( "a misspelled SP reference is found, with a did-you-mean suggestion" )
        {
            // parameters lock once set, so the broken reference goes on a fresh IWG
            tFem( FEM::IWG ).add_parameter_list();
            tFem.set( "IWG_name", "IWGTraction" );
            tFem.set( "stabilization_parameters", "SPNitcheDirichlet,DirichletNitsche" );

            Vector< std::string > tFindings = collect_deck_interlink_findings( tAll );
            REQUIRE( tFindings.size() == 1 );
            CHECK( tFindings( 0 ).find( "SPNitcheDirichlet" ) != std::string::npos );
            CHECK( tFindings( 0 ).find( "did you mean 'SPNitscheDirichlet'" ) != std::string::npos );
        }

        SECTION( "GEN IQI_types referencing an unknown IQI is found" )
        {
            tAll( (uint)Module_Type::GEN )( 0 )( 0 ).set( "IQI_types", Vector< std::string >{ "IQIBulkStrainEnergy", "IQIBulkVolume" } );

            Vector< std::string > tFindings = collect_deck_interlink_findings( tAll );
            REQUIRE( tFindings.size() == 1 );
            CHECK( tFindings( 0 ).find( "IQIBulkVolume" ) != std::string::npos );
        }

        SECTION( "ghost sets without XTK ghost_stab are found" )
        {
            tFem( FEM::IWG )( 0 ).set( "mesh_set_names", "ghost_p1" );

            Vector< std::string > tFindings = collect_deck_interlink_findings( tAll );
            REQUIRE( tFindings.size() == 1 );
            CHECK( tFindings( 0 ).find( "ghost_stab" ) != std::string::npos );

            // enabling ghost stabilization clears the finding
            tAll( (uint)Module_Type::XTK )( 0 )( 0 ).set( "ghost_stab", true );
            CHECK( collect_deck_interlink_findings( tAll ).size() == 0 );
        }

        SECTION( "VIS list-length mismatches are found" )
        {
            Module_Parameter_Lists& tVis = tAll( (uint)Module_Type::VIS );
            tVis( 0 )( 0 ).set( "Field_Names", "UX,UY" );
            tVis( 0 )( 0 ).set( "Field_Type", "NODAL" );
            tVis( 0 )( 0 ).set( "IQI_Names", "IQIBulkStrainEnergy" );

            Vector< std::string > tFindings = collect_deck_interlink_findings( tAll );
            REQUIRE( tFindings.size() == 1 );
            CHECK( tFindings( 0 ).find( "matching lengths" ) != std::string::npos );
        }
    }

    //------------------------------------------------------------------------------------------------------------------

}    // namespace moris
