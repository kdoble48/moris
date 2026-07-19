/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_Moment_Fitting.cpp
 *
 * Unit tests for the moment-fitted cut-cell quadrature
 * (fn_FEM_Moment_Fitting: Legendre basis, NNLS, per-cluster rule fit).
 *
 */

#include "catch.hpp"

#include <chrono>
#include <iostream>

#include "fn_FEM_Moment_Fitting.hpp"

#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Integration_Rule.hpp"
#include "cl_MTK_Integrator.hpp"

namespace moris::fem
{
    namespace
    {
        //------------------------------------------------------------------------------
        // exact integral of x^a y^b z^c over the unit simplex ( x,y,z >= 0, x+y+z <= 1 )
        real
        simplex_monomial_integral( uint aA, uint aB, uint aC )
        {
            // a! b! c! / ( a + b + c + 3 )!
            real tValue = 1.0;
            uint tDen   = 0;

            auto tFactorial = []( uint aN ) {
                real tF = 1.0;
                for ( uint iK = 2; iK <= aN; iK++ )
                {
                    tF *= iK;
                }
                return tF;
            };

            tDen = aA + aB + aC + 3;

            tValue = tFactorial( aA ) * tFactorial( aB ) * tFactorial( aC ) / tFactorial( tDen );

            return tValue;
        }

        //------------------------------------------------------------------------------
        // Kuhn (6-tet) tessellation of the box [x0,x1] x [-1,1] x [-1,1];
        // returns tets as ( 4 x 3 ) vertex-coordinate matrices, positively oriented
        // in the MORIS 3-coordinate tet parametrization
        Vector< Matrix< DDRMat > >
        kuhn_tessellation_of_box( real aX0, real aX1 )
        {
            // vertices of the box, index bit pattern ( x, y, z )
            Matrix< DDRMat > tV( 8, 3 );
            for ( uint iV = 0; iV < 8; iV++ )
            {
                tV( iV, 0 ) = ( iV & 1 ) ? aX1 : aX0;
                tV( iV, 1 ) = ( iV & 2 ) ? 1.0 : -1.0;
                tV( iV, 2 ) = ( iV & 4 ) ? 1.0 : -1.0;
            }

            // the 6 Kuhn tets along the vertex 0 -> 7 diagonal (all positively oriented
            // when the path permutations are even, fixed below by a volume check)
            Vector< Vector< uint > > tPaths = {
                { 0, 1, 3, 7 },
                { 0, 3, 2, 7 },
                { 0, 2, 6, 7 },
                { 0, 6, 4, 7 },
                { 0, 4, 5, 7 },
                { 0, 5, 1, 7 }
            };

            Vector< Matrix< DDRMat > > tTets( 6 );

            for ( uint iT = 0; iT < 6; iT++ )
            {
                Matrix< DDRMat > tTet( 4, 3 );
                for ( uint iV = 0; iV < 4; iV++ )
                {
                    for ( uint iD = 0; iD < 3; iD++ )
                    {
                        tTet( iV, iD ) = tV( tPaths( iT )( iV ), iD );
                    }
                }

                // enforce positive orientation in the MORIS tet convention
                Matrix< DDRMat > tMap, tOrigin;
                real             tDetJ = 0.0;
                moment_fitting::simplex_affine_map( tTet, tMap, tOrigin, tDetJ );

                if ( tDetJ < 0.0 )
                {
                    // swap two vertices
                    for ( uint iD = 0; iD < 3; iD++ )
                    {
                        std::swap( tTet( 0, iD ), tTet( 1, iD ) );
                    }
                }

                tTets( iT ) = tTet;
            }

            return tTets;
        }

        //------------------------------------------------------------------------------
        // evaluate a monomial x^a y^b z^c at the columns of aPoints
        real
        quadrature_of_monomial(
                const Matrix< DDRMat >& aPoints,
                const Matrix< DDRMat >& aWeights,
                uint                    aA,
                uint                    aB,
                uint                    aC )
        {
            real tSum = 0.0;
            for ( uint iP = 0; iP < aWeights.numel(); iP++ )
            {
                tSum += aWeights( iP )
                      * std::pow( aPoints( 0, iP ), aA )
                      * std::pow( aPoints( 1, iP ), aB )
                      * std::pow( aPoints( 2, iP ), aC );
            }
            return tSum;
        }

        //------------------------------------------------------------------------------
        // exact integral of x^a y^b z^c over the box [x0,x1] x [-1,1] x [-1,1]
        real
        box_monomial_integral( real aX0, real aX1, uint aA, uint aB, uint aC )
        {
            auto tMom1d = []( real aLo, real aHi, uint aExp ) {
                return ( std::pow( aHi, aExp + 1 ) - std::pow( aLo, aExp + 1 ) ) / ( aExp + 1.0 );
            };
            return tMom1d( aX0, aX1, aA ) * tMom1d( -1.0, 1.0, aB ) * tMom1d( -1.0, 1.0, aC );
        }

        //------------------------------------------------------------------------------
        // get the space part of a MORIS TET Gauss rule
        void
        tet_rule(
                mtk::Integration_Order aOrder,
                Matrix< DDRMat >&      aPoints,
                Matrix< DDRMat >&      aWeights )
        {
            mtk::Integration_Rule tRule(
                    mtk::Geometry_Type::TET,
                    mtk::Integration_Type::GAUSS,
                    aOrder,
                    mtk::Geometry_Type::POINT,
                    mtk::Integration_Type::GAUSS,
                    mtk::Integration_Order::POINT );

            mtk::Integrator tIntegrator( tRule );

            aPoints  = tIntegrator.get_space_points();
            aWeights = tIntegrator.get_space_weights();
        }
    }    // namespace

    //------------------------------------------------------------------------------

    TEST_CASE( "MomentFitting_LegendreBasis", "[moris],[fem],[MomentFitting]" )
    {
        // number of basis functions of the total-degree space
        REQUIRE( moment_fitting::number_of_basis_functions( 3, 0 ) == 1 );
        REQUIRE( moment_fitting::number_of_basis_functions( 3, 2 ) == 10 );
        REQUIRE( moment_fitting::number_of_basis_functions( 3, 4 ) == 35 );

        // spot check values: P2( x ) = ( 3 x^2 - 1 ) / 2
        Matrix< DDRMat > tPoints = { { 0.3 }, { -0.5 }, { 0.7 } };
        Matrix< DDRMat > tBasis;
        moment_fitting::evaluate_legendre_basis( tPoints, 2, tBasis );

        // ordering (matches the Phase-0 harness mflib): per total degree,
        // exponents (i,j,k) with i ascending, then j ascending, k = tot-i-j:
        // (0,0,0) | (0,0,1),(0,1,0),(1,0,0) | (0,0,2),(0,1,1),(0,2,0),(1,0,1),(1,1,0),(2,0,0)
        CHECK( tBasis( 0, 0 ) == Approx( 1.0 ) );
        CHECK( tBasis( 1, 0 ) == Approx( 0.7 ) );
        CHECK( tBasis( 2, 0 ) == Approx( -0.5 ) );
        CHECK( tBasis( 3, 0 ) == Approx( 0.3 ) );
        CHECK( tBasis( 4, 0 ) == Approx( ( 3.0 * 0.49 - 1.0 ) / 2.0 ) );
        CHECK( tBasis( 5, 0 ) == Approx( -0.5 * 0.7 ) );
        CHECK( tBasis( 9, 0 ) == Approx( ( 3.0 * 0.09 - 1.0 ) / 2.0 ) );
    }

    //------------------------------------------------------------------------------

    TEST_CASE( "MomentFitting_TetRuleDegrees", "[moris],[fem],[MomentFitting]" )
    {
        // verify the degree -> TET rule map used by Cluster::compute_moment_fitted_rule:
        // each rule must integrate all monomials of its assigned total degree exactly
        Vector< std::pair< mtk::Integration_Order, uint > > tMap = {
            { mtk::Integration_Order::TET_1, 1 },
            { mtk::Integration_Order::TET_4, 2 },
            { mtk::Integration_Order::TET_5, 3 },
            { mtk::Integration_Order::TET_11, 4 },
            { mtk::Integration_Order::TET_15, 5 },
            { mtk::Integration_Order::TET_35, 6 },
            { mtk::Integration_Order::TET_56, 8 }
        };

        for ( auto& tEntry : tMap )
        {
            Matrix< DDRMat > tPoints, tWeights;
            tet_rule( tEntry.first, tPoints, tWeights );

            const uint tDegree = tEntry.second;

            for ( uint iTot = 0; iTot <= tDegree; iTot++ )
            {
                for ( uint iA = 0; iA <= iTot; iA++ )
                {
                    for ( uint iB = 0; iB <= iTot - iA; iB++ )
                    {
                        const uint tC = iTot - iA - iB;

                        // MORIS tet weights sum to one; the reference volume 1/6
                        // is folded into the det J convention
                        real tQuad  = quadrature_of_monomial( tPoints, tWeights, iA, iB, tC ) / 6.0;
                        real tExact = simplex_monomial_integral( iA, iB, tC );

                        CHECK( tQuad == Approx( tExact ).margin( 1e-12 ) );
                    }
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    TEST_CASE( "MomentFitting_NNLS", "[moris],[fem],[MomentFitting]" )
    {
        // consistent system with a known non-negative solution: recovery to machine precision
        Matrix< DDRMat > tA = {
            { 1.0, 0.0, 1.0 },
            { 0.0, 1.0, 1.0 },
            { 1.0, 1.0, 0.0 }
        };
        Matrix< DDRMat > tXRef = { { 0.5 }, { 1.5 }, { 2.0 } };

        Matrix< DDRMat > tB( 3, 1, 0.0 );
        for ( uint iR = 0; iR < 3; iR++ )
        {
            for ( uint iC = 0; iC < 3; iC++ )
            {
                tB( iR ) += tA( iR, iC ) * tXRef( iC );
            }
        }

        Matrix< DDRMat > tX;
        real             tResidual = moment_fitting::nnls( tA, tB, tX );

        CHECK( tResidual < 1e-12 );
        for ( uint iC = 0; iC < 3; iC++ )
        {
            CHECK( tX( iC ) == Approx( tXRef( iC ) ).margin( 1e-12 ) );
        }

        // inconsistent system whose unconstrained solution is infeasible:
        // NNLS must return non-negative weights and match the known solution
        // min || [1 0; 0 1] x - [ -1; 2 ] ||, x >= 0  ->  x = ( 0, 2 )
        Matrix< DDRMat > tA2 = {
            { 1.0, 0.0 },
            { 0.0, 1.0 }
        };
        Matrix< DDRMat > tB2 = { { -1.0 }, { 2.0 } };

        Matrix< DDRMat > tX2;
        real             tResidual2 = moment_fitting::nnls( tA2, tB2, tX2 );

        CHECK( tX2( 0 ) == Approx( 0.0 ).margin( 1e-14 ) );
        CHECK( tX2( 1 ) == Approx( 2.0 ).margin( 1e-14 ) );
        CHECK( tResidual2 == Approx( 1.0 ).margin( 1e-12 ) );
    }

    //------------------------------------------------------------------------------

    TEST_CASE( "MomentFitting_UncutHex_RigorGate", "[moris],[fem],[MomentFitting]" )
    {
        // rigor gate (Phase 0): moment fitting on an uncut parent cell tessellated
        // into 6 Kuhn tets must reproduce all total-degree-d moments of the cell,
        // with non-negative weights and the NNLS active set pruning the candidates
        Vector< Matrix< DDRMat > > tTets = kuhn_tessellation_of_box( -1.0, 1.0 );

        moment_fitting::Fit_Input tInput;
        tInput.mDegree          = 4;
        tInput.mRelTol          = 1e-10;
        tInput.mCellParamCoords = tTets;

        tet_rule( mtk::Integration_Order::TET_35, tInput.mCandSpacePoints, tInput.mCandSpaceWeights );
        tet_rule( mtk::Integration_Order::TET_11, tInput.mMomSpacePoints, tInput.mMomSpaceWeights );

        moment_fitting::Fit_Result tFit = moment_fitting::fit_cluster_rule( tInput );

        REQUIRE( tFit.mSuccess );

        // material volume = 8 (the full cell)
        CHECK( tFit.mMaterialVolume == Approx( 8.0 ).margin( 1e-12 ) );
        CHECK( tFit.mResidual < 1e-10 );

        // weights are strictly positive and pruned below the candidate count
        real tWeightSum = 0.0;
        for ( uint iP = 0; iP < tFit.mWeights.numel(); iP++ )
        {
            CHECK( tFit.mWeights( iP ) > 0.0 );
            tWeightSum += tFit.mWeights( iP );
        }
        CHECK( tWeightSum == Approx( 8.0 ).margin( 1e-10 ) );
        CHECK( tFit.mWeights.numel() < tFit.mNumCandidates );

        // the fitted rule integrates every total-degree-4 monomial over the cell exactly
        for ( uint iTot = 0; iTot <= 4; iTot++ )
        {
            for ( uint iA = 0; iA <= iTot; iA++ )
            {
                for ( uint iB = 0; iB <= iTot - iA; iB++ )
                {
                    const uint tC     = iTot - iA - iB;
                    real       tQuad  = quadrature_of_monomial( tFit.mPoints, tFit.mWeights, iA, iB, tC );
                    real       tExact = box_monomial_integral( -1.0, 1.0, iA, iB, tC );

                    CHECK( tQuad == Approx( tExact ).margin( 1e-9 ) );
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    TEST_CASE( "MomentFitting_CutSlab", "[moris],[fem],[MomentFitting]" )
    {
        // plane-cut cell: material slab x in [-1, -1 + 2 alpha] of the parent cell;
        // simplicial tessellation is exact for plane cuts, so the fitted rule must
        // reproduce the analytic slab moments to the fitting tolerance
        for ( real tAlpha : { 0.5, 0.05, 1e-3 } )
        {
            const real tX1 = -1.0 + 2.0 * tAlpha;

            Vector< Matrix< DDRMat > > tTets = kuhn_tessellation_of_box( -1.0, tX1 );

            moment_fitting::Fit_Input tInput;
            tInput.mDegree          = 4;
            tInput.mRelTol          = 1e-10;
            tInput.mCellParamCoords = tTets;

            tet_rule( mtk::Integration_Order::TET_35, tInput.mCandSpacePoints, tInput.mCandSpaceWeights );
            tet_rule( mtk::Integration_Order::TET_11, tInput.mMomSpacePoints, tInput.mMomSpaceWeights );

            moment_fitting::Fit_Result tFit = moment_fitting::fit_cluster_rule( tInput );

            CAPTURE( tAlpha, tFit.mResidual, tFit.mMaterialVolume, tFit.mNumCandidates );
            REQUIRE( tFit.mSuccess );

            const real tVolume = 2.0 * tAlpha * 4.0;
            CHECK( tFit.mMaterialVolume == Approx( tVolume ).epsilon( 1e-10 ) );

            // all weights non-negative
            for ( uint iP = 0; iP < tFit.mWeights.numel(); iP++ )
            {
                CHECK( tFit.mWeights( iP ) > 0.0 );
            }

            // fitted rule reproduces the analytic monomial moments of the slab
            // relative to the moment vector scale
            real tScale = tVolume;

            for ( uint iTot = 0; iTot <= 4; iTot++ )
            {
                for ( uint iA = 0; iA <= iTot; iA++ )
                {
                    for ( uint iB = 0; iB <= iTot - iA; iB++ )
                    {
                        const uint tC     = iTot - iA - iB;
                        real       tQuad  = quadrature_of_monomial( tFit.mPoints, tFit.mWeights, iA, iB, tC );
                        real       tExact = box_monomial_integral( -1.0, tX1, iA, iB, tC );

                        CHECK( tQuad == Approx( tExact ).margin( 1e-8 * tScale ) );
                    }
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    TEST_CASE( "MomentFitting_HostTetConversion", "[moris],[fem],[MomentFitting]" )
    {
        // the per-cluster rule is evaluated through a host tet: verify that mapping
        // a parent point to host coordinates and back is the identity, and that the
        // det J convention matches the MORIS tet space interpolator ( det / 6 )
        Matrix< DDRMat > tTet = {
            { 0.2, -0.3, 0.1 },
            { 0.9, 0.4, -0.2 },
            { -0.5, -0.8, 0.3 },
            { 0.1, 0.2, 0.9 }
        };

        Matrix< DDRMat > tMap, tOrigin;
        real             tDetJ = 0.0;
        moment_fitting::simplex_affine_map( tTet, tMap, tOrigin, tDetJ );

        REQUIRE( std::abs( tDetJ ) > 1e-6 );

        // vertices map to the MORIS tet parametric coordinates
        // v0 -> (1,0,0), v1 -> (0,1,0), v2 -> (0,0,0), v3 -> (0,0,1)
        Matrix< DDRMat > tRef = {
            { 1.0, 0.0, 0.0 },
            { 0.0, 1.0, 0.0 },
            { 0.0, 0.0, 0.0 },
            { 0.0, 0.0, 1.0 }
        };

        for ( uint iV = 0; iV < 4; iV++ )
        {
            for ( uint iD = 0; iD < 3; iD++ )
            {
                real tXi = tOrigin( iD )
                         + tMap( iD, 0 ) * tRef( iV, 0 )
                         + tMap( iD, 1 ) * tRef( iV, 1 )
                         + tMap( iD, 2 ) * tRef( iV, 2 );

                CHECK( tXi == Approx( tTet( iV, iD ) ).margin( 1e-14 ) );
            }
        }

        // the tet volume equals det J ( weights of MORIS tet rules sum to one )
        // volume of the tet from the scalar triple product
        real tVolume =
                std::abs(
                        ( tTet( 1, 0 ) - tTet( 0, 0 ) ) * ( ( tTet( 2, 1 ) - tTet( 0, 1 ) ) * ( tTet( 3, 2 ) - tTet( 0, 2 ) ) - ( tTet( 2, 2 ) - tTet( 0, 2 ) ) * ( tTet( 3, 1 ) - tTet( 0, 1 ) ) )
                        - ( tTet( 1, 1 ) - tTet( 0, 1 ) ) * ( ( tTet( 2, 0 ) - tTet( 0, 0 ) ) * ( tTet( 3, 2 ) - tTet( 0, 2 ) ) - ( tTet( 2, 2 ) - tTet( 0, 2 ) ) * ( tTet( 3, 0 ) - tTet( 0, 0 ) ) )
                        + ( tTet( 1, 2 ) - tTet( 0, 2 ) ) * ( ( tTet( 2, 0 ) - tTet( 0, 0 ) ) * ( tTet( 3, 1 ) - tTet( 0, 1 ) ) - ( tTet( 2, 1 ) - tTet( 0, 1 ) ) * ( tTet( 3, 0 ) - tTet( 0, 0 ) ) ) )
                / 6.0;

        CHECK( std::abs( tDetJ ) == Approx( tVolume ).epsilon( 1e-12 ) );
    }

    //------------------------------------------------------------------------------
    // 2D helpers
    //------------------------------------------------------------------------------

    namespace
    {
        //------------------------------------------------------------------------------
        // exact integral of x^a y^b over the unit triangle ( x,y >= 0, x+y <= 1 )
        real
        triangle_monomial_integral( uint aA, uint aB )
        {
            auto tFactorial = []( uint aN ) {
                real tF = 1.0;
                for ( uint iK = 2; iK <= aN; iK++ )
                {
                    tF *= iK;
                }
                return tF;
            };

            // a! b! / ( a + b + 2 )!
            return tFactorial( aA ) * tFactorial( aB ) / tFactorial( aA + aB + 2 );
        }

        //------------------------------------------------------------------------------
        // 2-triangle tessellation of the box [x0,x1] x [-1,1], positively oriented
        // in the MORIS TRI parametrization
        Vector< Matrix< DDRMat > >
        tri_tessellation_of_box( real aX0, real aX1 )
        {
            Vector< Matrix< DDRMat > > tTris( 2 );

            tTris( 0 ) = { { aX1, -1.0 }, { aX1, 1.0 }, { aX0, -1.0 } };
            tTris( 1 ) = { { aX1, 1.0 }, { aX0, 1.0 }, { aX0, -1.0 } };

            for ( uint iT = 0; iT < 2; iT++ )
            {
                Matrix< DDRMat > tMap, tOrigin;
                real             tDetJ = 0.0;
                moment_fitting::simplex_affine_map( tTris( iT ), tMap, tOrigin, tDetJ );

                if ( tDetJ < 0.0 )
                {
                    for ( uint iD = 0; iD < 2; iD++ )
                    {
                        std::swap( tTris( iT )( 0, iD ), tTris( iT )( 1, iD ) );
                    }
                }
            }

            return tTris;
        }

        //------------------------------------------------------------------------------
        // evaluate a monomial x^a y^b at the columns of aPoints
        real
        quadrature_of_monomial_2d(
                const Matrix< DDRMat >& aPoints,
                const Matrix< DDRMat >& aWeights,
                uint                    aA,
                uint                    aB )
        {
            real tSum = 0.0;
            for ( uint iP = 0; iP < aWeights.numel(); iP++ )
            {
                tSum += aWeights( iP )
                      * std::pow( aPoints( 0, iP ), aA )
                      * std::pow( aPoints( 1, iP ), aB );
            }
            return tSum;
        }

        //------------------------------------------------------------------------------
        // exact integral of x^a y^b over the box [x0,x1] x [-1,1]
        real
        box_monomial_integral_2d( real aX0, real aX1, uint aA, uint aB )
        {
            auto tMom1d = []( real aLo, real aHi, uint aExp ) {
                return ( std::pow( aHi, aExp + 1 ) - std::pow( aLo, aExp + 1 ) ) / ( aExp + 1.0 );
            };
            return tMom1d( aX0, aX1, aA ) * tMom1d( -1.0, 1.0, aB );
        }

        //------------------------------------------------------------------------------
        // get the space part of a MORIS TRI Gauss rule
        void
        tri_rule(
                mtk::Integration_Order aOrder,
                Matrix< DDRMat >&      aPoints,
                Matrix< DDRMat >&      aWeights )
        {
            mtk::Integration_Rule tRule(
                    mtk::Geometry_Type::TRI,
                    mtk::Integration_Type::GAUSS,
                    aOrder,
                    mtk::Geometry_Type::POINT,
                    mtk::Integration_Type::GAUSS,
                    mtk::Integration_Order::POINT );

            mtk::Integrator tIntegrator( tRule );

            aPoints  = tIntegrator.get_space_points();
            aWeights = tIntegrator.get_space_weights();
        }
    }    // namespace

    //------------------------------------------------------------------------------

    TEST_CASE( "MomentFitting_LegendreBasis2D", "[moris],[fem],[MomentFitting]" )
    {
        REQUIRE( moment_fitting::number_of_basis_functions( 2, 0 ) == 1 );
        REQUIRE( moment_fitting::number_of_basis_functions( 2, 2 ) == 6 );
        REQUIRE( moment_fitting::number_of_basis_functions( 2, 5 ) == 21 );

        // spot check values at ( 0.3, -0.5 )
        Matrix< DDRMat > tPoints = { { 0.3 }, { -0.5 } };
        Matrix< DDRMat > tBasis;
        moment_fitting::evaluate_legendre_basis( tPoints, 2, tBasis );

        // ordering: (0,0) | (0,1),(1,0) | (0,2),(1,1),(2,0)
        CHECK( tBasis( 0, 0 ) == Approx( 1.0 ) );
        CHECK( tBasis( 1, 0 ) == Approx( -0.5 ) );
        CHECK( tBasis( 2, 0 ) == Approx( 0.3 ) );
        CHECK( tBasis( 3, 0 ) == Approx( ( 3.0 * 0.25 - 1.0 ) / 2.0 ) );
        CHECK( tBasis( 4, 0 ) == Approx( 0.3 * -0.5 ) );
        CHECK( tBasis( 5, 0 ) == Approx( ( 3.0 * 0.09 - 1.0 ) / 2.0 ) );
    }

    //------------------------------------------------------------------------------

    TEST_CASE( "MomentFitting_TriRuleDegrees", "[moris],[fem],[MomentFitting]" )
    {
        // verify the degree -> TRI rule map used by Cluster::compute_moment_fitted_rule
        Vector< std::pair< mtk::Integration_Order, uint > > tMap = {
            { mtk::Integration_Order::TRI_1, 1 },
            { mtk::Integration_Order::TRI_3, 2 },
            { mtk::Integration_Order::TRI_4, 3 },
            { mtk::Integration_Order::TRI_6, 4 },
            { mtk::Integration_Order::TRI_7, 5 },
            { mtk::Integration_Order::TRI_12, 6 },
            { mtk::Integration_Order::TRI_13, 7 },
            { mtk::Integration_Order::TRI_16, 8 }
        };

        for ( auto& tEntry : tMap )
        {
            Matrix< DDRMat > tPoints, tWeights;
            tri_rule( tEntry.first, tPoints, tWeights );

            const uint tDegree = tEntry.second;

            for ( uint iTot = 0; iTot <= tDegree; iTot++ )
            {
                for ( uint iA = 0; iA <= iTot; iA++ )
                {
                    // MORIS tri weights sum to one; the reference area 1/2
                    // is folded into the det J convention
                    real tQuad  = quadrature_of_monomial_2d( tPoints, tWeights, iA, iTot - iA ) / 2.0;
                    real tExact = triangle_monomial_integral( iA, iTot - iA );

                    CHECK( tQuad == Approx( tExact ).margin( 1e-12 ) );
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    TEST_CASE( "MomentFitting_UncutQuad_RigorGate2D", "[moris],[fem],[MomentFitting]" )
    {
        // 2D rigor gate: moment fitting on an uncut parent cell tessellated into
        // 2 triangles must reproduce all total-degree-d moments of the cell,
        // with non-negative weights and the NNLS active set pruning the candidates
        Vector< Matrix< DDRMat > > tTris = tri_tessellation_of_box( -1.0, 1.0 );

        moment_fitting::Fit_Input tInput;
        tInput.mDegree          = 5;
        tInput.mRelTol          = 1e-10;
        tInput.mCellParamCoords = tTris;

        tri_rule( mtk::Integration_Order::TRI_12, tInput.mCandSpacePoints, tInput.mCandSpaceWeights );
        tri_rule( mtk::Integration_Order::TRI_7, tInput.mMomSpacePoints, tInput.mMomSpaceWeights );

        moment_fitting::Fit_Result tFit = moment_fitting::fit_cluster_rule( tInput );

        REQUIRE( tFit.mSuccess );

        // material measure = 4 (the full cell)
        CHECK( tFit.mMaterialVolume == Approx( 4.0 ).margin( 1e-12 ) );
        CHECK( tFit.mResidual < 1e-10 );

        // weights are strictly positive and pruned below the candidate count
        real tWeightSum = 0.0;
        for ( uint iP = 0; iP < tFit.mWeights.numel(); iP++ )
        {
            CHECK( tFit.mWeights( iP ) > 0.0 );
            tWeightSum += tFit.mWeights( iP );
        }
        CHECK( tWeightSum == Approx( 4.0 ).margin( 1e-10 ) );
        CHECK( tFit.mWeights.numel() < tFit.mNumCandidates );

        // the fitted rule integrates every total-degree-5 monomial over the cell exactly
        for ( uint iTot = 0; iTot <= 5; iTot++ )
        {
            for ( uint iA = 0; iA <= iTot; iA++ )
            {
                real tQuad  = quadrature_of_monomial_2d( tFit.mPoints, tFit.mWeights, iA, iTot - iA );
                real tExact = box_monomial_integral_2d( -1.0, 1.0, iA, iTot - iA );

                CHECK( tQuad == Approx( tExact ).margin( 1e-9 ) );
            }
        }
    }

    //------------------------------------------------------------------------------

    TEST_CASE( "MomentFitting_CutSlab2D", "[moris],[fem],[MomentFitting]" )
    {
        // plane-cut 2D cell: material slab x in [-1, -1 + 2 alpha] of the parent cell
        for ( real tAlpha : { 0.5, 0.05, 1e-3 } )
        {
            const real tX1 = -1.0 + 2.0 * tAlpha;

            Vector< Matrix< DDRMat > > tTris = tri_tessellation_of_box( -1.0, tX1 );

            moment_fitting::Fit_Input tInput;
            tInput.mDegree          = 5;
            tInput.mRelTol          = 1e-10;
            tInput.mCellParamCoords = tTris;

            tri_rule( mtk::Integration_Order::TRI_12, tInput.mCandSpacePoints, tInput.mCandSpaceWeights );
            tri_rule( mtk::Integration_Order::TRI_7, tInput.mMomSpacePoints, tInput.mMomSpaceWeights );

            moment_fitting::Fit_Result tFit = moment_fitting::fit_cluster_rule( tInput );

            CAPTURE( tAlpha, tFit.mResidual, tFit.mMaterialVolume, tFit.mNumCandidates );
            REQUIRE( tFit.mSuccess );

            const real tArea = 2.0 * tAlpha * 2.0;
            CHECK( tFit.mMaterialVolume == Approx( tArea ).epsilon( 1e-10 ) );

            // all weights non-negative
            for ( uint iP = 0; iP < tFit.mWeights.numel(); iP++ )
            {
                CHECK( tFit.mWeights( iP ) > 0.0 );
            }

            // fitted rule reproduces the analytic monomial moments of the slab
            real tScale = tArea;

            for ( uint iTot = 0; iTot <= 5; iTot++ )
            {
                for ( uint iA = 0; iA <= iTot; iA++ )
                {
                    real tQuad  = quadrature_of_monomial_2d( tFit.mPoints, tFit.mWeights, iA, iTot - iA );
                    real tExact = box_monomial_integral_2d( -1.0, tX1, iA, iTot - iA );

                    CHECK( tQuad == Approx( tExact ).margin( 1e-8 * tScale ) );
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    TEST_CASE( "MomentFitting_HostTriConversion", "[moris],[fem],[MomentFitting]" )
    {
        // verify the MORIS TRI parametrization of the affine map:
        // v0 -> (1,0), v1 -> (0,1), v2 -> (0,0), det J convention det / 2
        Matrix< DDRMat > tTri = {
            { 0.7, -0.2 },
            { -0.1, 0.8 },
            { -0.6, -0.9 }
        };

        Matrix< DDRMat > tMap, tOrigin;
        real             tDetJ = 0.0;
        moment_fitting::simplex_affine_map( tTri, tMap, tOrigin, tDetJ );

        REQUIRE( std::abs( tDetJ ) > 1e-6 );

        Matrix< DDRMat > tRef = {
            { 1.0, 0.0 },
            { 0.0, 1.0 },
            { 0.0, 0.0 }
        };

        for ( uint iV = 0; iV < 3; iV++ )
        {
            for ( uint iD = 0; iD < 2; iD++ )
            {
                real tXi = tOrigin( iD )
                         + tMap( iD, 0 ) * tRef( iV, 0 )
                         + tMap( iD, 1 ) * tRef( iV, 1 );

                CHECK( tXi == Approx( tTri( iV, iD ) ).margin( 1e-14 ) );
            }
        }

        // the triangle area equals det J ( weights of MORIS tri rules sum to one )
        real tArea = std::abs(
                             ( tTri( 1, 0 ) - tTri( 0, 0 ) ) * ( tTri( 2, 1 ) - tTri( 0, 1 ) )
                             - ( tTri( 1, 1 ) - tTri( 0, 1 ) ) * ( tTri( 2, 0 ) - tTri( 0, 0 ) ) )
                   / 2.0;

        CHECK( std::abs( tDetJ ) == Approx( tArea ).epsilon( 1e-12 ) );
    }

    //------------------------------------------------------------------------------

    TEST_CASE( "MomentFitting_FitBench", "[moris],[fem],[MomentFitting]" )
    {
        // micro-benchmark of the per-cluster fit at production-like sizes
        // (3D: 30 tets x TET_35 = 1050 candidates; 2D: 32 tris x TRI_12 = 384);
        // prints ms/fit per (dim, degree); the loose assertions guard against
        // gross regressions - the before/after speedup gate is evaluated externally

        // 3D cluster: 5 x 1 x 1 boxes tessellated into 6 Kuhn tets each
        Vector< Matrix< DDRMat > > tTets;
        for ( uint iB = 0; iB < 5; iB++ )
        {
            const real tX0 = -1.0 + 2.0 * iB / 5.0;
            const real tX1 = -1.0 + 2.0 * ( iB + 1.0 ) / 5.0;

            Vector< Matrix< DDRMat > > tBox = kuhn_tessellation_of_box( tX0, tX1 );
            for ( uint iT = 0; iT < tBox.size(); iT++ )
            {
                tTets.push_back( tBox( iT ) );
            }
        }

        // 2D cluster: 16 boxes tessellated into 2 tris each
        Vector< Matrix< DDRMat > > tTris;
        for ( uint iB = 0; iB < 16; iB++ )
        {
            const real tX0 = -1.0 + 2.0 * iB / 16.0;
            const real tX1 = -1.0 + 2.0 * ( iB + 1.0 ) / 16.0;

            Vector< Matrix< DDRMat > > tBox = tri_tessellation_of_box( tX0, tX1 );
            for ( uint iT = 0; iT < tBox.size(); iT++ )
            {
                tTris.push_back( tBox( iT ) );
            }
        }

        auto tBench = [ & ]( uint aDim, uint aDegree, uint aReps ) {
            moment_fitting::Fit_Input tInput;
            tInput.mDegree = aDegree;
            tInput.mRelTol = 1e-10;

            if ( aDim == 3 )
            {
                tInput.mCellParamCoords = tTets;
                tet_rule( mtk::Integration_Order::TET_35, tInput.mCandSpacePoints, tInput.mCandSpaceWeights );
                tet_rule( aDegree <= 4 ? mtk::Integration_Order::TET_11 : mtk::Integration_Order::TET_35,
                        tInput.mMomSpacePoints,
                        tInput.mMomSpaceWeights );
            }
            else
            {
                tInput.mCellParamCoords = tTris;
                tri_rule( mtk::Integration_Order::TRI_12, tInput.mCandSpacePoints, tInput.mCandSpaceWeights );
                tri_rule( mtk::Integration_Order::TRI_12, tInput.mMomSpacePoints, tInput.mMomSpaceWeights );
            }

            // warm-up + correctness
            moment_fitting::Fit_Result tFit = moment_fitting::fit_cluster_rule( tInput );
            REQUIRE( tFit.mSuccess );

            const auto tStart = std::chrono::steady_clock::now();
            for ( uint iR = 0; iR < aReps; iR++ )
            {
                tFit = moment_fitting::fit_cluster_rule( tInput );
            }
            const auto tStop = std::chrono::steady_clock::now();

            const real tMs =
                    std::chrono::duration< real, std::milli >( tStop - tStart ).count() / aReps;

            std::cout << "MomentFitting_FitBench: dim " << aDim
                      << " degree " << aDegree
                      << " candidates " << tFit.mNumCandidates
                      << " retained " << tFit.mWeights.numel()
                      << " : " << tMs << " ms/fit" << std::endl;

            return tMs;
        };

        // 3D production-like sizes
        for ( uint tDeg : { 4u, 5u, 6u } )
        {
            const real tMs = tBench( 3, tDeg, 3 );
            CHECK( tMs < 500.0 );    // loose ceiling; the speedup gate is evaluated externally
        }

        // 3D sliver-conditioned cluster (the production-hard case: thin slab,
        // ill-conditioned candidate basis, many Lawson-Hanson iterations)
        tTets.clear();
        for ( uint iB = 0; iB < 5; iB++ )
        {
            const real tX0 = -1.0 + 0.002 * iB / 5.0;
            const real tX1 = -1.0 + 0.002 * ( iB + 1.0 ) / 5.0;

            Vector< Matrix< DDRMat > > tBox = kuhn_tessellation_of_box( tX0, tX1 );
            for ( uint iT = 0; iT < tBox.size(); iT++ )
            {
                tTets.push_back( tBox( iT ) );
            }
        }
        {
            const real tMs = tBench( 3, 5, 3 );
            CHECK( tMs < 500.0 );
        }

        // 2D
        for ( uint tDeg : { 4u, 5u, 6u } )
        {
            const real tMs = tBench( 2, tDeg, 5 );
            CHECK( tMs < 200.0 );
        }
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem
