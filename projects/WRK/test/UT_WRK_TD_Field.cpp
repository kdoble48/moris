/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_WRK_TD_Field.cpp
 *
 * ELEMENTAL FEM fields on a two-phase enriched (XFEM) mesh, with the FEM FIELDS key field_elemental_normalization
 * and IQI TOPOLOGICAL_DERIVATIVE, run through the HMR-XTK workflow on WRK_Input_TD_Field (4 x 2 background cells,
 * h = 1, interface x = 1.3, phase 0 left; exact piecewise-constant state, see the input file).
 *
 * - The fields live on the enriched interpolation mesh: one row per enriched interpolation cell (8 background
 *   cells + 2 for the second phase of the cut column = 10), each written once (the single-write assert in
 *   Interpolation_Element::populate_fields holds in debug builds); rows of the other phase stay unset.
 * - ip_cell divides the cluster integral by the interpolation-cell measure: VOLUME gives the cut ratio rho_e = 0.3
 *   on the cut column (1 elsewhere), and every ip_cell value is rho_e times the default (cluster average) value.
 * - TD values: g = delta dJ, dJ = -3(r-1)/(2r+1) sigma:eps + delta l (chi = 2), sigma = 1, eps_p = 1/E_p:
 *   phase 0 (E0 = 1, r = 0.25, delta = -1, l = 0.6): g = -(1.5 - 0.6) = -0.9;
 *   phase 1 (E1 = 0.25, r = 4, delta = +1): g = -1 * 4 + 0.6 = -3.4.
 *   Sum over the rows of both ip_cell fields times the cell measure = integral of g = 2.6 (-0.9) + 5.4 (-3.4).
 *
 */

#include <algorithm>
#include <cmath>

#include "catch.hpp"

#include "cl_Communication_Tools.hpp"
#include "paths.hpp"
#include "HDF5_Tools.hpp"

#include "cl_WRK_Performer_Manager.hpp"
#include "cl_WRK_Workflow_HMR_XTK.hpp"
#include "cl_Library_Factory.hpp"

using namespace moris;

namespace
{
    Matrix< DDRMat >
    read_td_field( const std::string& aFieldName )
    {
        Matrix< DDRMat > tValues;

        hid_t  tFile   = open_hdf5_file( "./WRK_TD_Field_" + aFieldName + ".hdf5" );
        herr_t tStatus = 0;
        load_matrix_from_hdf5_file( tFile, aFieldName, tValues, tStatus );
        tStatus = close_hdf5_file( tFile );

        REQUIRE( tStatus == 0 );
        return tValues;
    }

    bool
    is_set( real aValue )
    {
        return aValue != MORIS_REAL_MIN;
    }

    uint
    count_close( const Matrix< DDRMat >& aValues, real aTarget, real aTolerance )
    {
        uint tCount = 0;
        for ( uint iRow = 0; iRow < aValues.n_rows(); iRow++ )
        {
            if ( is_set( aValues( iRow ) ) && std::abs( aValues( iRow ) - aTarget ) < aTolerance )
            {
                tCount++;
            }
        }
        return tCount;
    }
}    // namespace

TEST_CASE( "WRK_TD_Field_Elemental_Normalization", "[moris],[WRK_TD_Field],[IQI_Topological_Derivative]" )
{
    if ( par_size() == 1 )
    {
        std::string tInputFilePath = moris::get_moris_bin_dir() + "/lib/WRK_Input_TD_Field.so";

        // load library and run the workflow (HMR, GEN, XTK, FEM solve, field output)
        moris::Library_Factory        tLibraryFactory;
        std::shared_ptr< Library_IO > tLibrary = tLibraryFactory.create_Library( Library_Type::STANDARD );
        tLibrary->load_parameter_list( tInputFilePath, File_Type::SO_FILE );
        tLibrary->finalize();

        {
            wrk::Performer_Manager tPerformerManager( tLibrary );
            wrk::Workflow_HMR_XTK  tWorkflow( &tPerformerManager );

            Vector< real >  tDummy( 1, 0.0 );
            Matrix< IdMat > tDummy1( 1, 1, 0.0 );
            tWorkflow.initialize( tDummy, tDummy, tDummy, tDummy1 );
            tWorkflow.perform( tDummy );
        }

        const real tTol = 1.0e-8;

        Matrix< DDRMat > tCutRatio = read_td_field( "CutRatio" );
        Matrix< DDRMat > tTDP0     = read_td_field( "TDP0" );
        Matrix< DDRMat > tTDP1     = read_td_field( "TDP1" );
        Matrix< DDRMat > tTDP0Avg  = read_td_field( "TDP0Avg" );

        // one row per enriched interpolation cell: 8 background cells + 2 for the second phase of the cut column
        REQUIRE( tCutRatio.n_rows() == 10 );
        REQUIRE( tTDP0.n_rows() == 10 );
        REQUIRE( tTDP1.n_rows() == 10 );
        REQUIRE( tTDP0Avg.n_rows() == 10 );

        // cut ratio: 1 on the two uncut phase-0 cells, 0.3 on the two cut cells, phase-1 rows unset
        CHECK( count_close( tCutRatio, 1.0, tTol ) == 2 );
        CHECK( count_close( tCutRatio, 0.3, tTol ) == 2 );

        // cluster average (default): -0.9 on all four phase-0 cells
        CHECK( count_close( tTDP0Avg, -0.9, tTol ) == 4 );

        // ip_cell: rho_e times the cluster average, row by row
        uint tNumPhase0 = 0;
        uint tNumPhase1 = 0;
        real tIntegral  = 0.0;
        for ( uint iRow = 0; iRow < tTDP0.n_rows(); iRow++ )
        {
            // each enriched cell belongs to exactly one phase
            CHECK( is_set( tTDP0( iRow ) ) != is_set( tTDP1( iRow ) ) );

            if ( is_set( tTDP0( iRow ) ) )
            {
                tNumPhase0++;
                REQUIRE( is_set( tCutRatio( iRow ) ) );
                REQUIRE( is_set( tTDP0Avg( iRow ) ) );
                CHECK( std::abs( tTDP0( iRow ) - tCutRatio( iRow ) * tTDP0Avg( iRow ) ) < tTol );
                tIntegral += tTDP0( iRow );
            }
            else
            {
                tNumPhase1++;
                CHECK( !is_set( tCutRatio( iRow ) ) );
                tIntegral += tTDP1( iRow );
            }
        }
        CHECK( tNumPhase0 == 4 );
        CHECK( tNumPhase1 == 6 );

        // phase-0 values: -0.9 uncut, 0.3 * -0.9 cut; phase-1 values: -3.4 uncut, 0.7 * -3.4 cut
        CHECK( count_close( tTDP0, -0.9, tTol ) == 2 );
        CHECK( count_close( tTDP0, -0.27, tTol ) == 2 );
        CHECK( count_close( tTDP1, -3.4, tTol ) == 4 );
        CHECK( count_close( tTDP1, -2.38, tTol ) == 2 );

        // per cut background cell: TD_p0 + TD_p1 = rho_e g_in + ( 1 - rho_e ) g_out = -2.65; over the domain the
        // ip_cell rows (cell measure 1) sum to the integral of g
        CHECK( std::abs( tIntegral - ( 2.6 * -0.9 + 5.4 * -3.4 ) ) < 1.0e-7 );
    }
}
