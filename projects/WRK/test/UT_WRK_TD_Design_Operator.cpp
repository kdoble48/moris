/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_WRK_TD_Design_Operator.cpp
 *
 * wrk::TD_Design_Operator (MORISGENERAL sublist TD_OPERATOR): b_A = sum_e g_e int_e B_A and the Helmholtz filter
 * ( gamma K + M ) g~ = b on the design B-spline space (sandbox-topopt IMPL_PLAN class B, tests U6-U11).
 *
 * Operator cases: HMR on D = [0,2] x [0,1], 8 x 8 cells (h_x = 2 h_y, the 2:1 cells of the cantilever), Lagrange /
 * B-spline orders (1,1), (2,2), (3,3), (2,1); GEN level sets discretized on the B-spline mesh: "Const" (-1/sqrt(2)) and
 * "LineX" (phi = x, whose coefficients are the Greville abscissae in x).
 *   U6  sum M = |D| = 2; row sums of M = int B_A; K 1 = 0; c^T M c = int x^2 = 8/3 and c^T M 1 = int x = 2 (c = LineX)
 *   U7  straight interface x0 = 0.9 (cut ratio 0.6 in the cut column): sum_A b_A = int g exactly; away from the cut
 *       column b_A = g_side int B_A
 *   U8  g = const -> g~ = const; gamma = 0 -> M g~ = b; gamma > 0 -> ( gamma K + M ) g~ = b
 *   U9  constant start: |c_A + 1/sqrt(2)| < 1e-10 and |phi|_M = 1
 *   U10 no zero rows, positive Cholesky pivots, every coefficient an unknown or a dropped (zero T column) one
 * Workflow case (U11): WRK_Input_TD_Operator (two-phase XFEM bar, exact state, TD values -0.9 / -3.4), two
 * evaluations with the interface at x = 1.3 and 1.6, for the two-phase route (TDP0 + TDP1) and the cut-ratio route
 * (TDP0 + ( 1 - rho_e ) void value): g_e per cell, sum b = int g, cached direction of the last evaluated ADVs, operator
 * built once.
 */

#include <algorithm>
#include <cmath>
#include <cstdlib>

#include "catch.hpp"

#include "cl_Communication_Tools.hpp"
#include "paths.hpp"
#include "HDF5_Tools.hpp"

#include "cl_HMR.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Mesh_Pair.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Field.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_Module_Parameter_Lists.hpp"
#include "cl_Library_Enums.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_PRM_MORIS_GENERAL_Parameters.hpp"

#include "cl_WRK_TD_Design_Operator.hpp"
#include "cl_WRK_Performer_Manager.hpp"
#include "cl_WRK_Workflow_HMR_XTK.hpp"
#include "cl_Library_Factory.hpp"

using namespace moris;

namespace
{
    //------------------------------------------------------------------------------

    struct TD_Test_Setup
    {
        std::shared_ptr< hmr::HMR >             mHMR;
        std::shared_ptr< mtk::Mesh_Manager >    mMeshManager;
        std::shared_ptr< gen::Geometry_Engine > mGeometryEngine;
        Vector< real >                          mADVs;

        mtk::Interpolation_Mesh*
        get_mesh()
        {
            return mMeshManager->get_mesh_pair( 0 ).get_interpolation_mesh();
        }
    };

    //------------------------------------------------------------------------------

    void
    add_design( Module_Parameter_Lists& aGENParameters, gen::Field_Type aType, const std::string& aName )
    {
        aGENParameters( GEN_Submodule::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( aType ) );
        aGENParameters.set( "name", aName );
        aGENParameters.set( "discretization_mesh_index", 0 );
        aGENParameters.set( "discretization_lower_bound", -1.0e3 );
        aGENParameters.set( "discretization_upper_bound", 1.0e3 );
    }

    //------------------------------------------------------------------------------

    /** HMR on [0,2] x [0,1] with one Lagrange and one B-spline mesh, GEN with the designs Const and LineX */
    TD_Test_Setup
    create_setup( uint aLagrangeOrder, uint aBSplineOrder, uint aNumX, uint aNumY )
    {
        TD_Test_Setup tSetup;

        Module_Parameter_Lists tHMRParameters( Module_Type::HMR );
        tHMRParameters.set( "number_of_elements_per_dimension", aNumX, aNumY );
        tHMRParameters.set( "domain_dimensions", 2.0, 1.0 );
        tHMRParameters.set( "domain_offset", 0.0, 0.0 );
        tHMRParameters.set( "lagrange_output_meshes", "0" );
        tHMRParameters.set( "lagrange_orders", std::to_string( aLagrangeOrder ) );
        tHMRParameters.set( "lagrange_pattern", "0" );
        tHMRParameters.set( "bspline_orders", std::to_string( aBSplineOrder ) );
        tHMRParameters.set( "bspline_pattern", "0" );
        tHMRParameters.set( "lagrange_to_bspline", "0" );
        tHMRParameters.set( "refinement_buffer", 1 );
        tHMRParameters.set( "staircase_buffer", 1 );

        tSetup.mHMR         = std::make_shared< hmr::HMR >( tHMRParameters );
        tSetup.mMeshManager = std::make_shared< mtk::Mesh_Manager >();
        tSetup.mHMR->set_performer( tSetup.mMeshManager );
        tSetup.mHMR->perform_initial_refinement();
        tSetup.mHMR->perform();

        Module_Parameter_Lists tGENParameters( Module_Type::GEN );
        add_design( tGENParameters, gen::Field_Type::CONSTANT, "Const" );
        tGENParameters.set( "constant", -1.0 / std::sqrt( 2.0 ) );
        add_design( tGENParameters, gen::Field_Type::LINE, "LineX" );    // phi = normal . ( x - center ) = x

        tSetup.mGeometryEngine = std::make_shared< gen::Geometry_Engine >( tGENParameters );
        tSetup.mGeometryEngine->distribute_advs( tSetup.mMeshManager->get_mesh_pair( 0 ) );
        tSetup.mADVs = tSetup.mGeometryEngine->get_advs();

        return tSetup;
    }

    //------------------------------------------------------------------------------

    Parameter_List
    create_td_parameters( const std::string& aDesignName, real aGamma )
    {
        Parameter_List tParameters = prm::create_td_operator_parameter_list();
        tParameters.set( "adv_field", aDesignName );
        tParameters.set( "td_field_names", "NotUsedHere" );
        tParameters.set( "filter_gamma", aGamma );
        return tParameters;
    }

    //------------------------------------------------------------------------------

    real
    max_abs_difference( const Vector< real >& aFirst, const Vector< real >& aSecond )
    {
        real tMax = 0.0;
        for ( uint i = 0; i < aFirst.size(); i++ )
        {
            tMax = std::max( tMax, std::abs( aFirst( i ) - aSecond( i ) ) );
        }
        return tMax;
    }

    real
    max_abs( const Vector< real >& aVector )
    {
        real tMax = 0.0;
        for ( real iValue : aVector )
        {
            tMax = std::max( tMax, std::abs( iValue ) );
        }
        return tMax;
    }

    real
    sum( const Vector< real >& aVector )
    {
        real tSum = 0.0;
        for ( real iValue : aVector )
        {
            tSum += iValue;
        }
        return tSum;
    }

    //------------------------------------------------------------------------------

    /** x range of a background cell */
    void
    cell_x_range( const mtk::Cell& aCell, real& aXMin, real& aXMax )
    {
        aXMin = MORIS_REAL_MAX;
        aXMax = -MORIS_REAL_MAX;
        for ( mtk::Vertex* iVertex : aCell.get_vertex_pointers() )
        {
            aXMin = std::min( aXMin, iVertex->get_coords()( 0 ) );
            aXMax = std::max( aXMax, iVertex->get_coords()( 0 ) );
        }
    }

    //------------------------------------------------------------------------------

    Vector< real >
    read_hdf5_vector( hid_t aFile, const std::string& aLabel )
    {
        Vector< real > tVector;
        herr_t         tStatus = 0;
        load_vector_from_hdf5_file( aFile, aLabel, tVector.data(), tStatus );
        REQUIRE( tStatus == 0 );
        return tVector;
    }

    //------------------------------------------------------------------------------
}    // namespace

//------------------------------------------------------------------------------

TEST_CASE( "WRK_TD_Design_Operator_Parameters", "[moris],[WRK_TD_Design_Operator]" )
{
    // the TD_OPERATOR sublist is the fourth MORISGENERAL submodule
    CHECK( get_number_of_sub_parameter_lists_in_module( Module_Type::MORISGENERAL ) == 4 );
    CHECK( static_cast< uint >( MORISGENERAL_Submodule::TD_OPERATOR ) == 3 );
    CHECK( get_submodule_names( Module_Type::MORISGENERAL )( 3 ) == "TD_Operator" );

    Parameter_List tParameters = prm::create_td_operator_parameter_list();
    CHECK( tParameters.get< std::string >( "adv_field" ).empty() );
    CHECK( tParameters.get< std::string >( "td_field_names" ).empty() );
    CHECK( tParameters.get< std::string >( "cut_ratio_field" ).empty() );
    CHECK( tParameters.get< real >( "void_td_value" ) == 0.0 );
    CHECK( tParameters.get< real >( "filter_gamma" ) == 1.0e-4 );
    CHECK( tParameters.get< std::string >( "output_file" ).empty() );
}

//------------------------------------------------------------------------------

TEST_CASE( "WRK_TD_Design_Operator_Matrices", "[moris],[WRK_TD_Design_Operator]" )
{
    if ( par_size() != 1 )
    {
        return;
    }

    // ( Lagrange order, B-spline order )
    const uint tOrders[ 4 ][ 2 ] = { { 1, 1 }, { 2, 2 }, { 3, 3 }, { 2, 1 } };

    for ( uint iCase = 0; iCase < 4; iCase++ )
    {
        uint tLagrangeOrder = tOrders[ iCase ][ 0 ];
        uint tBSplineOrder  = tOrders[ iCase ][ 1 ];

        SECTION( "Lagrange order " + std::to_string( tLagrangeOrder ) + ", B-spline order " + std::to_string( tBSplineOrder ) )
        {
            TD_Test_Setup            tSetup    = create_setup( tLagrangeOrder, tBSplineOrder, 8, 8 );
            mtk::Interpolation_Mesh* tMesh     = tSetup.get_mesh();
            uint                     tNumADVs  = tSetup.mADVs.size();
            uint                     tNumCells = tMesh->get_num_elems();
            Vector< real >           tOnes( tNumADVs, 1.0 );

            // two designs of ( 8 + d )^2 coefficients each
            REQUIRE( tNumADVs == 2 * ( 8 + tBSplineOrder ) * ( 8 + tBSplineOrder ) );

            wrk::TD_Design_Operator tLineX( create_td_parameters( "LineX", 1.0e-2 ) );
            tLineX.build_operators( *tSetup.mGeometryEngine, tNumADVs );

            //------------------------------------------------------------------------------
            // U10: every coefficient of the design is an unknown or dropped; SPD

            Vector< sint > tLinePositions = tSetup.mGeometryEngine->get_design_adv_positions( "LineX" );
            CHECK( tLineX.get_num_unknowns() + tLineX.get_num_dropped_coefficients() == tLinePositions.size() );
            CHECK( tLineX.get_num_unknowns() > 0 );
            CHECK( tLineX.get_min_diagonal() > 0.0 );
            CHECK( tLineX.get_min_pivot() > 0.0 );

            //------------------------------------------------------------------------------
            // U6: partition of unity, row sums, kernel of K, linear reproduction

            CHECK( std::abs( tLineX.get_domain_measure() - 2.0 ) < 1.0e-12 );
            CHECK( std::abs( tLineX.get_mass_sum() - 2.0 ) < 1.0e-12 );

            Vector< real > tRowSums   = tLineX.mass_times( tOnes );
            Vector< real > tIntegrals = tLineX.get_coefficient_integrals();
            CHECK( max_abs_difference( tRowSums, tIntegrals ) < 1.0e-14 );
            CHECK( std::abs( sum( tIntegrals ) - 2.0 ) < 1.0e-12 );

            CHECK( max_abs( tLineX.stiffness_times( tOnes ) ) < 1.0e-11 );

            // c = Greville abscissae in x: c^T M c = int_D x^2 = 8/3, c^T M 1 = int_D x = 2
            CHECK( std::abs( tLineX.mass_product( tSetup.mADVs, tSetup.mADVs ) - 8.0 / 3.0 ) < 1.0e-11 );
            CHECK( std::abs( tLineX.mass_product( tSetup.mADVs, tOnes ) - 2.0 ) < 1.0e-11 );

            // K c_x: int grad B_A . grad x = int d B_A / dx = boundary terms only; sum over A is zero
            CHECK( std::abs( sum( tLineX.stiffness_times( tSetup.mADVs ) ) ) < 1.0e-11 );

            //------------------------------------------------------------------------------
            // U9: constant start projected by GEN, unit M-norm

            wrk::TD_Design_Operator tConst( create_td_parameters( "Const", 1.0e-2 ) );
            tConst.build_operators( *tSetup.mGeometryEngine, tNumADVs );

            Vector< sint > tConstPositions = tSetup.mGeometryEngine->get_design_adv_positions( "Const" );
            real           tMaxDeviation   = 0.0;
            for ( sint iPosition : tConstPositions )
            {
                tMaxDeviation = std::max( tMaxDeviation, std::abs( tSetup.mADVs( iPosition ) + 1.0 / std::sqrt( 2.0 ) ) );
            }
            CHECK( tMaxDeviation < 1.0e-10 );
            CHECK( std::abs( std::sqrt( tConst.mass_product( tSetup.mADVs, tSetup.mADVs ) ) - 1.0 ) < 1.0e-10 );

            //------------------------------------------------------------------------------
            // U8: constants pass the filter; the filter solves its system; gamma = 0 is the L2 projection

            Vector< real > tConstantLoad = tConst.assemble_load( Vector< real >( tNumCells, 3.0 ) );
            Vector< real > tFiltered     = tConst.filter( tConstantLoad );
            real           tMaxConstErr  = 0.0;
            for ( sint iPosition : tConstPositions )
            {
                tMaxConstErr = std::max( tMaxConstErr, std::abs( tFiltered( iPosition ) - 3.0 ) );
            }
            CHECK( tMaxConstErr < 1.0e-10 );

            // a non-trivial cell field
            Vector< real > tCellValues( tNumCells );
            for ( uint iCell = 0; iCell < tNumCells; iCell++ )
            {
                real tXMin, tXMax;
                cell_x_range( tMesh->get_mtk_cell( iCell ), tXMin, tXMax );
                tCellValues( tMesh->get_mtk_cell( iCell ).get_index() ) = std::sin( 3.0 * tXMin ) + 0.1 * iCell;
            }
            Vector< real > tLoad = tConst.assemble_load( tCellValues );
            Vector< real > tG    = tConst.filter( tLoad );

            Vector< real > tResidual = tConst.mass_times( tG );
            Vector< real > tKG       = tConst.stiffness_times( tG );
            for ( uint i = 0; i < tNumADVs; i++ )
            {
                tResidual( i ) += tConst.get_filter_gamma() * tKG( i ) - tLoad( i );
            }
            CHECK( max_abs( tResidual ) < 1.0e-12 * max_abs( tLoad ) );

            wrk::TD_Design_Operator tProjection( create_td_parameters( "Const", 0.0 ) );
            tProjection.build_operators( *tSetup.mGeometryEngine, tNumADVs );
            Vector< real > tProjected = tProjection.filter( tLoad );
            CHECK( max_abs_difference( tProjection.mass_times( tProjected ), tLoad ) < 1.0e-12 * max_abs( tLoad ) );

            // the M-inner product is symmetric and the load is the same for any design on the same mesh
            CHECK( std::abs( tConst.mass_product( tG, tProjected ) - tConst.mass_product( tProjected, tG ) ) < 1.0e-12 );

            //------------------------------------------------------------------------------
            // U7: straight interface at x0 = 0.9 (cut column [0.75, 1.0], cut ratio 0.6)

            const real tX0   = 0.9;
            const real tGIn  = 2.0;
            const real tGOut = -5.0;

            Vector< real > tInterfaceValues( tNumCells, 0.0 );
            Vector< real > tCutIndicator( tNumCells, 0.0 );
            Vector< real > tLeftIndicator( tNumCells, 0.0 );
            Vector< real > tRightIndicator( tNumCells, 0.0 );
            for ( uint iCell = 0; iCell < tNumCells; iCell++ )
            {
                const mtk::Cell& tCell = tMesh->get_mtk_cell( iCell );
                real             tXMin, tXMax;
                cell_x_range( tCell, tXMin, tXMax );
                real tRho = std::min( 1.0, std::max( 0.0, ( tX0 - tXMin ) / ( tXMax - tXMin ) ) );

                tInterfaceValues( tCell.get_index() ) = tRho * tGIn + ( 1.0 - tRho ) * tGOut;
                tCutIndicator( tCell.get_index() )    = ( tRho > 0.0 and tRho < 1.0 ) ? 1.0 : 0.0;
                tLeftIndicator( tCell.get_index() )   = tRho == 1.0 ? 1.0 : 0.0;
                tRightIndicator( tCell.get_index() )  = tRho == 0.0 ? 1.0 : 0.0;
            }
            CHECK( sum( tCutIndicator ) == 8.0 );

            Vector< real > tInterfaceLoad = tLineX.assemble_load( tInterfaceValues );
            CHECK( std::abs( sum( tInterfaceLoad ) - ( tGIn * tX0 + tGOut * ( 2.0 - tX0 ) ) ) < 1.0e-12 );

            Vector< real > tCutLoad   = tLineX.assemble_load( tCutIndicator );
            Vector< real > tLeftLoad  = tLineX.assemble_load( tLeftIndicator );
            Vector< real > tRightLoad = tLineX.assemble_load( tRightIndicator );

            uint tNumLeft  = 0;
            uint tNumRight = 0;
            real tMaxErr   = 0.0;
            for ( sint iPosition : tLinePositions )
            {
                if ( tCutLoad( iPosition ) != 0.0 )
                {
                    continue;    // the cut band: differs from the exact sub-cell load by O( [g] )
                }
                if ( tLeftLoad( iPosition ) > 0.0 and tRightLoad( iPosition ) == 0.0 )
                {
                    tMaxErr = std::max( tMaxErr, std::abs( tInterfaceLoad( iPosition ) - tGIn * tIntegrals( iPosition ) ) );
                    tNumLeft++;
                }
                if ( tRightLoad( iPosition ) > 0.0 and tLeftLoad( iPosition ) == 0.0 )
                {
                    tMaxErr = std::max( tMaxErr, std::abs( tInterfaceLoad( iPosition ) - tGOut * tIntegrals( iPosition ) ) );
                    tNumRight++;
                }
            }
            CHECK( tNumLeft > 0 );
            CHECK( tNumRight > 0 );
            CHECK( tMaxErr < 1.0e-14 );
        }
    }
}

//------------------------------------------------------------------------------

namespace
{
    /** expected g_e of the bar for an interface at aX0 in the column [1, 2] (cut ratio aX0 - 1) */
    real
    expected_cell_value( real aCentroidX, real aX0 )
    {
        if ( aCentroidX < 1.0 )
        {
            return -0.9;
        }
        if ( aCentroidX < 2.0 )
        {
            real tRho = aX0 - 1.0;
            return tRho * -0.9 + ( 1.0 - tRho ) * -3.4;
        }
        return -3.4;
    }

    void
    check_bar_evaluation(
            wrk::Workflow_HMR_XTK&                            aWorkflow,
            const std::shared_ptr< wrk::TD_Design_Operator >& aOperator,
            const Vector< real >&                             aADVs,
            real                                              aX0,
            uint                                              aEvaluation )
    {
        CHECK( aOperator->get_num_evaluations() == aEvaluation );
        CHECK( aOperator->get_num_builds() == 1 );

        // the cache belongs to the ADVs of the last evaluation
        CHECK( max_abs_difference( aOperator->get_evaluated_advs(), aADVs ) == 0.0 );

        // output file of this evaluation: per-cell values by centroid
        hid_t            tFile       = open_hdf5_file( "./WRK_TD_Operator_" + std::to_string( aEvaluation ) + ".hdf5" );
        Vector< real >   tCellValues = read_hdf5_vector( tFile, "CellValues" );
        Vector< real >   tLoad       = read_hdf5_vector( tFile, "Load" );
        Vector< real >   tDirection  = read_hdf5_vector( tFile, "Direction" );
        Vector< real >   tADVs       = read_hdf5_vector( tFile, "ADVs" );
        Matrix< DDRMat > tCentroids;
        herr_t           tStatus = 0;
        load_matrix_from_hdf5_file( tFile, "CellCentroids", tCentroids, tStatus );
        close_hdf5_file( tFile );

        REQUIRE( tCellValues.size() == 8 );
        REQUIRE( tCentroids.n_rows() == 8 );
        for ( uint iCell = 0; iCell < 8; iCell++ )
        {
            CHECK( std::abs( tCellValues( iCell ) - expected_cell_value( tCentroids( iCell, 0 ), aX0 ) ) < 1.0e-8 );
        }

        // sum b = int g ( unit cells ): phase 0 area 2 aX0, phase 1 area 8 - 2 aX0
        CHECK( std::abs( sum( tLoad ) - ( 2.0 * aX0 * -0.9 + ( 8.0 - 2.0 * aX0 ) * -3.4 ) ) < 1.0e-8 );

        // the written and cached results agree; the direction solves the filter system
        CHECK( max_abs_difference( tDirection, aWorkflow.get_design_direction() ) == 0.0 );
        CHECK( max_abs_difference( tLoad, aOperator->get_load() ) == 0.0 );
        CHECK( max_abs_difference( tADVs, aADVs ) == 0.0 );
        CHECK( max_abs_difference( aOperator->filter( tLoad ), tDirection ) == 0.0 );

        // the M-inner product of the workflow: |D| = 8
        Vector< real > tOnes( aADVs.size(), 1.0 );
        CHECK( std::abs( aWorkflow.mass_inner_product( tOnes, tOnes ) - 8.0 ) < 1.0e-12 );
    }

    void
    run_bar( const std::string& aMode )
    {
        if ( aMode.empty() )
        {
            unsetenv( "WRK_TD_OPERATOR_MODE" );
        }
        else
        {
            setenv( "WRK_TD_OPERATOR_MODE", aMode.c_str(), 1 );
        }

        std::string tInputFilePath = moris::get_moris_bin_dir() + "/lib/WRK_Input_TD_Operator.so";

        moris::Library_Factory        tLibraryFactory;
        std::shared_ptr< Library_IO > tLibrary = tLibraryFactory.create_Library( Library_Type::STANDARD );
        tLibrary->load_parameter_list( tInputFilePath, File_Type::SO_FILE );
        tLibrary->finalize();

        unsetenv( "WRK_TD_OPERATOR_MODE" );

        wrk::Performer_Manager tPerformerManager( tLibrary );
        wrk::Workflow_HMR_XTK  tWorkflow( &tPerformerManager );

        std::shared_ptr< wrk::TD_Design_Operator > tOperator = tWorkflow.get_td_design_operator();
        REQUIRE( tOperator != nullptr );
        CHECK( tOperator->get_num_evaluations() == 0 );

        Vector< real >  tADVs;
        Vector< real >  tLowerBounds;
        Vector< real >  tUpperBounds;
        Matrix< IdMat > tIjklIds;
        tWorkflow.initialize( tADVs, tLowerBounds, tUpperBounds, tIjklIds );

        // linear B-splines on the 4 x 2 mesh: 15 coefficients = nodal values of x - 1.3
        REQUIRE( tADVs.size() == 15 );

        // evaluation 1: interface x = 1.3
        Vector< real > tADVs1 = tADVs;
        tWorkflow.perform( tADVs1 );
        check_bar_evaluation( tWorkflow, tOperator, tADVs1, 1.3, 1 );
        Vector< real > tDirection1 = tWorkflow.get_design_direction();

        // evaluation 2: interface x = 1.6, the cache follows
        Vector< real > tADVs2 = tADVs;
        for ( real& iADV : tADVs2 )
        {
            iADV -= 0.3;
        }
        tWorkflow.perform( tADVs2 );
        check_bar_evaluation( tWorkflow, tOperator, tADVs2, 1.6, 2 );
        CHECK( max_abs_difference( tDirection1, tWorkflow.get_design_direction() ) > 1.0e-3 );
    }
}    // namespace

TEST_CASE( "WRK_TD_Design_Operator_Workflow", "[moris],[WRK_TD_Design_Operator]" )
{
    if ( par_size() == 1 )
    {
        SECTION( "two-phase route: TDP0 + TDP1" )
        {
            run_bar( "" );
        }
        SECTION( "cut-ratio route: TDP0 + ( 1 - rho_e ) void value" )
        {
            run_bar( "cut_ratio" );
        }
    }
}
