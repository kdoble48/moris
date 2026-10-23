/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_GEN_ADV_Filter.cpp
 *
 * Kreissl-Maute linear filter of level-set B-spline coefficients (level-set key adv_filter_radius): the fields use
 * F * psi with F_ij = max( 0, r - |x_i - x_j| ) / sum_j, on the Greville points of the coefficients; the optimizer
 * keeps psi, and the criteria sensitivities are F^T times the coefficient sensitivities.
 *
 */

#include "catch.hpp"
#include <cmath>
#include <random>

#include "cl_GEN_Geometry_Engine_Test.hpp"
#include "cl_GEN_Design_Factory.hpp"
#include "cl_GEN_Derived_Node.hpp"
#include "fn_GEN_create_simple_mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"

#include "fn_PRM_GEN_Parameters.hpp"

namespace moris::gen
{
    //--------------------------------------------------------------------------------------------------------------

    // 8 x 8 cells on [-1, 1]^2: coefficient spacing h = 0.25 (uniform B-splines, Greville points on a lattice)
    static const uint tNumElements = 8;
    static const real tSpacing     = 0.25;

    //--------------------------------------------------------------------------------------------------------------

    static Submodule_Parameter_Lists
    create_level_set_parameters(
            bool aLinear,
            real aFilterRadius )
    {
        Submodule_Parameter_Lists tFieldParameterLists( "FIELDS" );
        if ( aLinear )
        {
            // phi = x - 0.1: its B-spline coefficients are the Greville abscissae minus 0.1
            tFieldParameterLists.add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
            tFieldParameterLists.set( "center_x", 0.1 );
            tFieldParameterLists.set( "center_y", 0.0 );
            tFieldParameterLists.set( "normal_x", 1.0 );
            tFieldParameterLists.set( "normal_y", 0.0 );
        }
        else
        {
            tFieldParameterLists.add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
            tFieldParameterLists.set( "center_x", 0.1 );
            tFieldParameterLists.set( "center_y", -0.05 );
            tFieldParameterLists.set( "radius", 0.63 );
        }
        tFieldParameterLists.set( "discretization_mesh_index", 0 );
        tFieldParameterLists.set( "discretization_lower_bound", -10.0 );
        tFieldParameterLists.set( "discretization_upper_bound", 10.0 );
        tFieldParameterLists.set( "adv_filter_radius", aFilterRadius );
        return tFieldParameterLists;
    }

    //--------------------------------------------------------------------------------------------------------------

    static Vector< real >
    random_vector(
            uint aSize,
            uint aSeed )
    {
        std::mt19937                     tGenerator( aSeed );
        std::uniform_real_distribution<> tDistribution( -1.0, 1.0 );
        Vector< real >                   tVector( aSize );
        for ( uint iEntry = 0; iEntry < aSize; iEntry++ )
        {
            tVector( iEntry ) = tDistribution( tGenerator );
        }
        return tVector;
    }

    //--------------------------------------------------------------------------------------------------------------

    static real
    dot(
            const Vector< real >& aFirst,
            const Vector< real >& aSecond )
    {
        real tSum = 0.0;
        for ( uint iEntry = 0; iEntry < aFirst.size(); iEntry++ )
        {
            tSum += aFirst( iEntry ) * aSecond( iEntry );
        }
        return tSum;
    }

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "ADV Filter", "[gen], [adv_filter]" )
    {
        if ( par_size() == 1 )
        {
            Derived_Node::clear_override_linear_interpolation();
            mtk::Interpolation_Mesh* tMesh = create_simple_mesh( tNumElements, tNumElements, 2, 2 );

            SECTION( "Off by default and for a radius below the coefficient spacing: F = I" )
            {
                for ( real tRadius : { 0.0, 0.9 * tSpacing } )
                {
                    Submodule_Parameter_Lists  tFieldParameterLists = create_level_set_parameters( false, tRadius );
                    Geometry_Engine_Parameters tGeometryEngineParameters;
                    ADV_Manager                tADVManager;
                    Design_Factory             tDesignFactory( tFieldParameterLists, tADVManager );
                    tGeometryEngineParameters.mGeometries = tDesignFactory.get_geometries();
                    Geometry_Engine_Test tGeometryEngine( tMesh, tGeometryEngineParameters );

                    CHECK( tGeometryEngine.has_adv_filter() == ( tRadius > 0.0 ) );
                    Vector< real > tADVs     = random_vector( tGeometryEngine.get_advs().size(), 1 );
                    Vector< real > tFiltered = tGeometryEngine.apply_adv_filter( tADVs );
                    Vector< real > tAdjoint  = tGeometryEngine.apply_adv_filter_transpose( tADVs );
                    for ( uint iADV = 0; iADV < tADVs.size(); iADV++ )
                    {
                        CHECK( tFiltered( iADV ) == tADVs( iADV ) );
                        CHECK( tAdjoint( iADV ) == tADVs( iADV ) );
                    }
                }
            }

            SECTION( "Radius 1.5 h: constants and interior linear fields are reproduced, F^T is the adjoint" )
            {
                // circle: generic coefficients
                Submodule_Parameter_Lists  tFieldParameterLists = create_level_set_parameters( false, 1.5 * tSpacing );
                Geometry_Engine_Parameters tGeometryEngineParameters;
                ADV_Manager                tADVManager;
                Design_Factory             tDesignFactory( tFieldParameterLists, tADVManager );
                tGeometryEngineParameters.mGeometries = tDesignFactory.get_geometries();
                Geometry_Engine_Test tGeometryEngine( tMesh, tGeometryEngineParameters );
                REQUIRE( tGeometryEngine.has_adv_filter() );
                uint tNumADVs = tGeometryEngine.get_advs().size();

                // F 1 = 1 (row normalization, also at the edges)
                Vector< real > tConstant( tNumADVs, 0.37 );
                Vector< real > tFilteredConstant = tGeometryEngine.apply_adv_filter( tConstant );
                for ( uint iADV = 0; iADV < tNumADVs; iADV++ )
                {
                    CHECK( tFilteredConstant( iADV ) == Approx( 0.37 ).margin( 1e-14 ) );
                }

                // < F a, b > = < a, F^T b >
                Vector< real > tA     = random_vector( tNumADVs, 2 );
                Vector< real > tB     = random_vector( tNumADVs, 3 );
                real           tLeft  = dot( tGeometryEngine.apply_adv_filter( tA ), tB );
                real           tRight = dot( tA, tGeometryEngine.apply_adv_filter_transpose( tB ) );
                CHECK( tLeft == Approx( tRight ).epsilon( 1e-13 ) );

                // the filter changes a generic vector (it is not the identity at this radius)
                Vector< real > tFilteredA = tGeometryEngine.apply_adv_filter( tA );
                real           tChange    = 0.0;
                for ( uint iADV = 0; iADV < tNumADVs; iADV++ )
                {
                    tChange = std::max( tChange, std::abs( tFilteredA( iADV ) - tA( iADV ) ) );
                }
                CHECK( tChange > 0.1 );

                // set_advs keeps psi for get_advs, and the field evaluates F psi: compare with an unfiltered engine
                // given F psi directly
                Submodule_Parameter_Lists  tPlainParameterLists = create_level_set_parameters( false, 0.0 );
                Geometry_Engine_Parameters tPlainParameters;
                ADV_Manager                tPlainADVManager;
                Design_Factory             tPlainDesignFactory( tPlainParameterLists, tPlainADVManager );
                tPlainParameters.mGeometries = tPlainDesignFactory.get_geometries();
                Geometry_Engine_Test tPlainEngine( tMesh, tPlainParameters );

                Vector< real > tPsi = tGeometryEngine.get_advs();
                for ( uint iADV = 0; iADV < tNumADVs; iADV++ )
                {
                    tPsi( iADV ) += 0.05 * tA( iADV );
                }
                tGeometryEngine.set_advs( tPsi );
                tPlainEngine.set_advs( tGeometryEngine.apply_adv_filter( tPsi ) );

                Vector< real > tReturned = tGeometryEngine.get_advs();
                for ( uint iADV = 0; iADV < tNumADVs; iADV++ )
                {
                    CHECK( tReturned( iADV ) == tPsi( iADV ) );
                }
                for ( uint iNode = 0; iNode < tMesh->get_num_nodes(); iNode++ )
                {
                    Matrix< DDRMat > tCoordinates = tMesh->get_node_coordinate( iNode );
                    real             tFilteredValue =
                            tGeometryEngine.get_geometry( 0 )->get_field()->get_field_value( iNode, tCoordinates );
                    real tPlainValue = tPlainEngine.get_geometry( 0 )->get_field()->get_field_value( iNode, tCoordinates );
                    CHECK( tFilteredValue == Approx( tPlainValue ).margin( 1e-12 ) );
                }
            }

            SECTION( "Radius 1.5 h: a linear field is reproduced wherever the filter support is complete in x" )
            {
                Submodule_Parameter_Lists  tFieldParameterLists = create_level_set_parameters( true, 1.5 * tSpacing );
                Geometry_Engine_Parameters tGeometryEngineParameters;
                ADV_Manager                tADVManager;
                Design_Factory             tDesignFactory( tFieldParameterLists, tADVManager );
                tGeometryEngineParameters.mGeometries = tDesignFactory.get_geometries();
                Geometry_Engine_Test tGeometryEngine( tMesh, tGeometryEngineParameters );

                // initial coefficients = Greville x - 0.1; a row whose cone is symmetric in x reproduces them
                Vector< real > tPsi           = tGeometryEngine.get_advs();
                Vector< real > tFiltered      = tGeometryEngine.apply_adv_filter( tPsi );
                uint           tNumReproduced = 0;
                uint           tNumInterior   = 0;
                for ( uint iADV = 0; iADV < tPsi.size(); iADV++ )
                {
                    real tX = tPsi( iADV ) + 0.1;
                    if ( std::abs( tX ) < 1.0 + 0.5 * tSpacing - 1e-9 )    // not on the outermost Greville column
                    {
                        tNumInterior++;
                        if ( std::abs( tFiltered( iADV ) - tPsi( iADV ) ) < 1e-10 )
                        {
                            tNumReproduced++;
                        }
                    }
                }
                std::printf( "ADV filter: linear field reproduced on %u of %u coefficients off the outermost columns\n", tNumReproduced, tNumInterior );
                CHECK( tNumInterior > 0 );
                CHECK( tNumReproduced == tNumInterior );
            }

            delete tMesh;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::gen
