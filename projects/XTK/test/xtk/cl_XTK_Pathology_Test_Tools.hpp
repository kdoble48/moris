/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Pathology_Test_Tools.hpp
 *
 * Test-only helpers shared by the XTK pathology sweeps (UT_XTK_Cut_Pathologies.cpp) and the connectivity tests
 * (cl_XTK_Connectivity_Test_Harness.hpp): field regimes, the HMR background box, and FNV-1a hashing with gathering
 * to rank 0 for multi-rank comparisons.
 *
 */

#pragma once

#include <array>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>

#include "cl_Communication_Tools.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_HMR_Background_Mesh_Base.hpp"
#include "cl_HMR_Background_Element_Base.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "cl_XTK_Connectivity_Checker.hpp"

namespace moris::xtk
{
    namespace pathology_test
    {
        using connectivity_check::Hash;
        using connectivity_check::Point;

        enum class Regime
        {
            STORED,
            ANALYTIC,
            MULTILINEAR
        };

        inline std::string
        regime_name( Regime aRegime )
        {
            switch ( aRegime )
            {
                case Regime::STORED:
                    return "stored";
                case Regime::ANALYTIC:
                    return "analytic";
                default:
                    return "multilinear";
            }
        }

        // ----------------------------------------------------------------------------------
        // background meshes: parameter list + perform() + Mesh_Manager (the setter path segfaults on > 1 rank)

        struct Grid_Mesh
        {
            std::shared_ptr< hmr::HMR >          mHMR;
            std::shared_ptr< mtk::Mesh_Manager > mMeshManager;
            mtk::Interpolation_Mesh*             mMesh = nullptr;

            // aRefineAboveX: background cells whose centre lies at x > aRefineAboveX are refined one level
            // aBsplineOrders / aBsplinePatterns / aLagrangeToBspline: HMR strings; aRefineFirstPattern: initial refinement of pattern 0
            // (pattern 1 unrefined). E.g. "1,2", "0,1", "0,1", 1: linear B-splines on the Lagrange mesh plus quadratic ones one level coarser
            Grid_Mesh(
                    const std::array< uint, 3 >& aElements,
                    const Point&                 aDimensions,
                    real                         aRefineAboveX      = MORIS_REAL_MAX,
                    const std::string&           aBsplineOrders     = "1",
                    const std::string&           aBsplinePatterns   = "0",
                    const std::string&           aLagrangeToBspline = "0",
                    uint                         aRefineFirstPattern = 0 )
            {
                this->build( aElements, aDimensions, aBsplineOrders, aBsplinePatterns, aLagrangeToBspline, aRefineFirstPattern, [ aRefineAboveX ]( const Matrix< DDRMat >& aCentre, uint ) { return aCentre( 0 ) > aRefineAboveX; }, aRefineAboveX < MORIS_REAL_MAX ? 1 : 0 );
            }

            // aRefine( centre, pass ): background cells whose centre satisfies it are refined one level, in aNumPasses passes
            Grid_Mesh(
                    const std::array< uint, 3 >&                                  aElements,
                    const Point&                                                  aDimensions,
                    const std::function< bool( const Matrix< DDRMat >&, uint ) >& aRefine,
                    uint                                                          aNumPasses )
            {
                this->build( aElements, aDimensions, "1", "0", "0", 0, aRefine, aNumPasses );
            }

          private:
            void
            build(
                    const std::array< uint, 3 >&                                  aElements,
                    const Point&                                                  aDimensions,
                    const std::string&                                            aBsplineOrders,
                    const std::string&                                            aBsplinePatterns,
                    const std::string&                                            aLagrangeToBspline,
                    uint                                                          aRefineFirstPattern,
                    const std::function< bool( const Matrix< DDRMat >&, uint ) >& aRefine,
                    uint                                                          aNumPasses )
            {
                Module_Parameter_Lists tParameters( Module_Type::HMR );
                tParameters( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );
                tParameters( 0 )( 0 ).set( "number_of_elements_per_dimension", aElements[ 0 ], aElements[ 1 ], aElements[ 2 ] );
                tParameters( 0 )( 0 ).set( "domain_dimensions", aDimensions[ 0 ], aDimensions[ 1 ], aDimensions[ 2 ] );
                tParameters( 0 )( 0 ).set( "domain_offset", 0.0, 0.0, 0.0 );
                tParameters( 0 )( 0 ).set( "lagrange_output_meshes", "0" );
                tParameters( 0 )( 0 ).set( "lagrange_orders", "1" );
                tParameters( 0 )( 0 ).set( "lagrange_pattern", std::string( "0" ) );
                tParameters( 0 )( 0 ).set( "bspline_orders", aBsplineOrders );
                tParameters( 0 )( 0 ).set( "bspline_pattern", aBsplinePatterns );
                tParameters( 0 )( 0 ).set( "lagrange_to_bspline", aLagrangeToBspline );
                tParameters( 0 )( 0 ).set( "refinement_buffer", 1 );
                tParameters( 0 )( 0 ).set( "staircase_buffer", 1 );
                if ( aRefineFirstPattern > 0 )
                {
                    tParameters( 0 )( 0 ).set( "pattern_initial_refinement", aRefineFirstPattern, 0u );    // pattern 0 refined, pattern 1 not
                }

                mHMR         = std::make_shared< hmr::HMR >( tParameters );
                mMeshManager = std::make_shared< mtk::Mesh_Manager >();
                mHMR->set_performer( mMeshManager );
                mHMR->perform_initial_refinement();
                for ( uint iPass = 0; iPass < aNumPasses; iPass++ )
                {
                    // flag as Database::flag_element does, without HMR::flag_element (MB-043: rank-local flag deadlocks)
                    hmr::Background_Mesh_Base* tBackgroundMesh = mHMR->get_database()->get_background_mesh();
                    luint                      tNumElements    = tBackgroundMesh->get_number_of_active_elements_on_proc();
                    Matrix< DDRMat >           tCentre( 3, 1 );
                    uint                       tWorkingPattern = mHMR->get_database()->get_parameters()->get_working_pattern();
                    for ( luint iElement = 0; iElement < tNumElements; iElement++ )
                    {
                        hmr::Background_Element_Base* tElement = tBackgroundMesh->get_element( iElement );
                        tBackgroundMesh->calc_center_of_element( tElement, tCentre );
                        if ( aRefine( tCentre, iPass ) )
                        {
                            tElement->put_on_refinement_queue();
                            tElement->set_refined_flag( tWorkingPattern );
                        }
                    }
                    mHMR->perform_refinement_based_on_working_pattern( 0 );
                }
                mHMR->perform();
                mMesh = mMeshManager->get_interpolation_mesh( 0 );
            }
        };

        // ----------------------------------------------------------------------------------
        // hashing (64-bit FNV-1a, split into two exactly representable reals for gathering as DDRMat)

        // rows gathered on rank 0 (empty elsewhere)
        inline Vector< Matrix< DDRMat > >
        gather_rows( const Vector< Vector< real > >& aRows, uint aNumColumns )
        {
            Matrix< DDRMat > tRows( aRows.size(), aNumColumns );
            for ( uint iRow = 0; iRow < aRows.size(); iRow++ )
            {
                for ( uint iCol = 0; iCol < aNumColumns; iCol++ )
                {
                    tRows( iRow, iCol ) = aRows( iRow )( iCol );
                }
            }
            Vector< Matrix< DDRMat > > tGathered;
            gatherv_mats( tRows, tGathered, 0 );
            return tGathered;
        }
    }    // namespace pathology_test
}    // namespace moris::xtk
