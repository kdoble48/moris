/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_FEM_Moment_Fitting.cpp
 *
 */

#include "fn_FEM_Moment_Fitting.hpp"

#include "fn_linsolve.hpp"    // LINALG/src - solve_least_squares

#include <cmath>
#include <limits>

namespace moris::fem::moment_fitting
{
    //------------------------------------------------------------------------------

    uint
    number_of_basis_functions( uint aDim, uint aDegree )
    {
        MORIS_ASSERT( aDim == 2 || aDim == 3,
                "moment_fitting::number_of_basis_functions - only 2D and 3D supported." );

        if ( aDim == 2 )
        {
            return ( aDegree + 1 ) * ( aDegree + 2 ) / 2;
        }

        return ( aDegree + 1 ) * ( aDegree + 2 ) * ( aDegree + 3 ) / 6;
    }

    //------------------------------------------------------------------------------

    namespace
    {
        /**
         * 1D Legendre polynomials P_0 ... P_d at the entries of one coordinate row
         * @param[ in ]  aPoints points ( dim x n )
         * @param[ in ]  aRow    coordinate row index
         * @param[ in ]  aDegree maximum degree d
         * @param[ out ] aP      values ( d+1 x n )
         */
        void
        legendre_1d(
                const Matrix< DDRMat >& aPoints,
                uint                    aRow,
                uint                    aDegree,
                Matrix< DDRMat >&       aP )
        {
            const uint tN = aPoints.n_cols();
            aP.set_size( aDegree + 1, tN );

            for ( uint iP = 0; iP < tN; iP++ )
            {
                const real tX = aPoints( aRow, iP );

                aP( 0, iP ) = 1.0;
                if ( aDegree >= 1 )
                {
                    aP( 1, iP ) = tX;
                }
                for ( uint iK = 1; iK < aDegree; iK++ )
                {
                    // ( k + 1 ) P_{k+1} = ( 2 k + 1 ) x P_k - k P_{k-1}
                    aP( iK + 1, iP ) =
                            ( ( 2.0 * iK + 1.0 ) * tX * aP( iK, iP ) - iK * aP( iK - 1, iP ) ) / ( iK + 1.0 );
                }
            }
        }
    }    // namespace

    //------------------------------------------------------------------------------

    void
    evaluate_legendre_basis(
            const Matrix< DDRMat >& aPoints,
            uint                    aDegree,
            Matrix< DDRMat >&       aBasis )
    {
        const uint tDim = aPoints.n_rows();

        MORIS_ASSERT( tDim == 2 || tDim == 3,
                "moment_fitting::evaluate_legendre_basis - points must be 2 x n or 3 x n." );

        const uint tN = aPoints.n_cols();
        const uint tM = number_of_basis_functions( tDim, aDegree );

        // 1D Legendre values per direction
        Matrix< DDRMat > tPx, tPy, tPz;
        legendre_1d( aPoints, 0, aDegree, tPx );
        legendre_1d( aPoints, 1, aDegree, tPy );
        if ( tDim == 3 )
        {
            legendre_1d( aPoints, 2, aDegree, tPz );
        }

        aBasis.set_size( tM, tN );

        // total-degree ordering (see header)
        uint tRow = 0;
        for ( uint iTot = 0; iTot <= aDegree; iTot++ )
        {
            for ( uint iI = 0; iI <= iTot; iI++ )
            {
                if ( tDim == 2 )
                {
                    const uint tJ = iTot - iI;
                    for ( uint iP = 0; iP < tN; iP++ )
                    {
                        aBasis( tRow, iP ) = tPx( iI, iP ) * tPy( tJ, iP );
                    }
                    tRow++;
                }
                else
                {
                    for ( uint iJ = 0; iJ <= iTot - iI; iJ++ )
                    {
                        const uint tK = iTot - iI - iJ;
                        for ( uint iP = 0; iP < tN; iP++ )
                        {
                            aBasis( tRow, iP ) = tPx( iI, iP ) * tPy( iJ, iP ) * tPz( tK, iP );
                        }
                        tRow++;
                    }
                }
            }
        }

        MORIS_ASSERT( tRow == tM,
                "moment_fitting::evaluate_legendre_basis - inconsistent basis count." );
    }

    //------------------------------------------------------------------------------

    namespace
    {
        /**
         * fast passive-set least-squares solve via the backend's compiled
         * LAPACK path (LINALG solve_least_squares): passive columns are
         * gathered into a persistent, geometrically-grown workspace, so the
         * repeated inner solves of Lawson-Hanson allocate nothing on the
         * common path
         * @param[ in ]     aA       full matrix ( m x n )
         * @param[ in ]     aPassive column indices of the passive set ( size k )
         * @param[ in ]     aB       right-hand side ( m x 1 )
         * @param[ in,out ] aApWork  gather workspace ( m x >= k ), grown geometrically
         * @param[ out ]    aZ       subproblem solution ( k x 1 )
         * @return false if the backend flags (numerical) rank deficiency or is
         *         unavailable - the caller then uses the MGS path, which
         *         preserves the original rank-deficiency handling
         */
        bool
        solve_passive_ls_fast(
                const Matrix< DDRMat >& aA,
                const Vector< uint >&   aPassive,
                const Matrix< DDRMat >& aB,
                Matrix< DDRMat >&       aApWork,
                Matrix< DDRMat >&       aZ )
        {
            const uint tM = aA.n_rows();
            const uint tK = aPassive.size();

            // keep the rank-handling MGS path for non-tall subproblems
            if ( tK == 0 || tK > tM )
            {
                return false;
            }

            // grow the gather workspace geometrically, never shrink
            if ( aApWork.n_rows() != tM || aApWork.n_cols() < tK )
            {
                aApWork.set_size( tM, std::max( tK, 2 * (uint)aApWork.n_cols() ) );
            }

            for ( uint iC = 0; iC < tK; iC++ )
            {
                const uint tCol = aPassive( iC );
                for ( uint iR = 0; iR < tM; iR++ )
                {
                    aApWork( iR, iC ) = aA( iR, tCol );
                }
            }

            // compiled LAPACK least squares on the k leading workspace columns
            if ( !solve_least_squares( aApWork, tK, aB, aZ ) )
            {
                return false;
            }

            // reject non-finite solutions (near-rank-deficient dgels artifacts)
            for ( uint iC = 0; iC < tK; iC++ )
            {
                if ( !std::isfinite( aZ( iC ) ) )
                {
                    return false;
                }
            }

            return true;
        }

        /**
         * least-squares solve of the passive-set subproblem min || A_P z - b ||
         * by modified Gram-Schmidt QR with one reorthogonalization pass
         * (rank-deficiency-detecting FALLBACK path of solve_passive_ls_fast)
         * @param[ in ]  aA       full matrix ( m x n )
         * @param[ in ]  aPassive column indices of the passive set ( size k )
         * @param[ in ]  aB       right-hand side ( m x 1 )
         * @param[ out ] aZ       subproblem solution ( k x 1 )
         * @return false if the subproblem is numerically rank deficient
         */
        bool
        solve_passive_ls(
                const Matrix< DDRMat >& aA,
                const Vector< uint >&   aPassive,
                const Matrix< DDRMat >& aB,
                Matrix< DDRMat >&       aZ )
        {
            const uint tM = aA.n_rows();
            const uint tK = aPassive.size();

            // copy passive columns
            Matrix< DDRMat > tQ( tM, tK );
            for ( uint iC = 0; iC < tK; iC++ )
            {
                for ( uint iR = 0; iR < tM; iR++ )
                {
                    tQ( iR, iC ) = aA( iR, aPassive( iC ) );
                }
            }

            // modified Gram-Schmidt with reorthogonalization
            Matrix< DDRMat > tR( tK, tK, 0.0 );
            for ( uint iC = 0; iC < tK; iC++ )
            {
                for ( uint iPass = 0; iPass < 2; iPass++ )
                {
                    for ( uint iPrev = 0; iPrev < iC; iPrev++ )
                    {
                        real tProj = 0.0;
                        for ( uint iR = 0; iR < tM; iR++ )
                        {
                            tProj += tQ( iR, iPrev ) * tQ( iR, iC );
                        }
                        for ( uint iR = 0; iR < tM; iR++ )
                        {
                            tQ( iR, iC ) -= tProj * tQ( iR, iPrev );
                        }
                        tR( iPrev, iC ) += tProj;
                    }
                }

                real tNorm = 0.0;
                for ( uint iR = 0; iR < tM; iR++ )
                {
                    tNorm += tQ( iR, iC ) * tQ( iR, iC );
                }
                tNorm = std::sqrt( tNorm );

                if ( tNorm < 1e-13 )
                {
                    return false;
                }

                tR( iC, iC ) = tNorm;
                for ( uint iR = 0; iR < tM; iR++ )
                {
                    tQ( iR, iC ) /= tNorm;
                }
            }

            // z = R^-1 Q^T b
            Matrix< DDRMat > tQtB( tK, 1, 0.0 );
            for ( uint iC = 0; iC < tK; iC++ )
            {
                for ( uint iR = 0; iR < tM; iR++ )
                {
                    tQtB( iC ) += tQ( iR, iC ) * aB( iR );
                }
            }

            aZ.set_size( tK, 1 );
            for ( sint iC = tK - 1; iC >= 0; iC-- )
            {
                real tSum = tQtB( iC );
                for ( uint iNext = iC + 1; iNext < tK; iNext++ )
                {
                    tSum -= tR( iC, iNext ) * aZ( iNext );
                }
                aZ( iC ) = tSum / tR( iC, iC );
            }

            return true;
        }
    }    // namespace

    //------------------------------------------------------------------------------

    real
    nnls(
            const Matrix< DDRMat >& aA,
            const Matrix< DDRMat >& aB,
            Matrix< DDRMat >&       aX,
            uint                    aMaxIter )
    {
        const uint tM = aA.n_rows();
        const uint tN = aA.n_cols();

        if ( aMaxIter == 0 )
        {
            aMaxIter = 10 * tN;
        }

        aX.set_size( tN, 1, 0.0 );

        // active/passive bookkeeping
        Vector< uint > tInPassive( tN, 0 );
        Vector< uint > tPassive;
        tPassive.reserve( tM );

        // residual r = b - A x  (x = 0 initially)
        Matrix< DDRMat > tRes( tM, 1 );
        for ( uint iR = 0; iR < tM; iR++ )
        {
            tRes( iR ) = aB( iR );
        }

        // gradient tolerance relative to the problem scale
        real tBNorm = 0.0;
        for ( uint iR = 0; iR < tM; iR++ )
        {
            tBNorm += aB( iR ) * aB( iR );
        }
        tBNorm = std::sqrt( tBNorm );

        const real tGradTol = 1e-13 * std::max( tBNorm, 1.0e-300 );

        Matrix< DDRMat > tZ;

        // persistent gather workspace for the fast passive-set solves
        Matrix< DDRMat > tApWork;

        // stagnation guard: terminate when the residual stops decreasing
        // (degenerate candidate sets can otherwise cycle add/remove of a column)
        real tPrevResNorm = tBNorm;

        for ( uint iOuter = 0; iOuter < aMaxIter; iOuter++ )
        {
            // gradient w = A^T r on the active set; pick the most positive entry
            sint tBest     = -1;
            real tBestGrad = tGradTol;
            for ( uint iC = 0; iC < tN; iC++ )
            {
                if ( tInPassive( iC ) )
                {
                    continue;
                }
                real tGrad = 0.0;
                for ( uint iR = 0; iR < tM; iR++ )
                {
                    tGrad += aA( iR, iC ) * tRes( iR );
                }
                if ( tGrad > tBestGrad )
                {
                    tBestGrad = tGrad;
                    tBest     = iC;
                }
            }

            // optimality reached
            if ( tBest < 0 )
            {
                break;
            }

            tInPassive( tBest ) = 1;
            tPassive.push_back( (uint)tBest );

            // inner loop: restore feasibility of the passive set
            for ( uint iInner = 0; iInner < aMaxIter; iInner++ )
            {
                if ( !solve_passive_ls_fast( aA, tPassive, aB, tApWork, tZ ) &&    //
                        !solve_passive_ls( aA, tPassive, aB, tZ ) )
                {
                    // rank-deficient subproblem: drop the newest column and stop growing
                    tInPassive( tPassive( tPassive.size() - 1 ) ) = 0;
                    tPassive.pop_back();
                    break;
                }

                // check feasibility of the subproblem solution
                bool tFeasible = true;
                for ( uint iP = 0; iP < tPassive.size(); iP++ )
                {
                    if ( tZ( iP ) <= 0.0 )
                    {
                        tFeasible = false;
                        break;
                    }
                }

                if ( tFeasible )
                {
                    for ( uint iP = 0; iP < tPassive.size(); iP++ )
                    {
                        aX( tPassive( iP ) ) = tZ( iP );
                    }
                    break;
                }

                // step towards z until the first passive variable hits zero
                real tAlpha = 1.0;
                for ( uint iP = 0; iP < tPassive.size(); iP++ )
                {
                    if ( tZ( iP ) <= 0.0 )
                    {
                        const real tXi = aX( tPassive( iP ) );
                        tAlpha         = std::min( tAlpha, tXi / ( tXi - tZ( iP ) ) );
                    }
                }

                for ( uint iP = 0; iP < tPassive.size(); iP++ )
                {
                    const uint tCol = tPassive( iP );
                    aX( tCol ) += tAlpha * ( tZ( iP ) - aX( tCol ) );
                }

                // move variables at zero back to the active set
                Vector< uint > tNewPassive;
                tNewPassive.reserve( tPassive.size() );
                for ( uint iP = 0; iP < tPassive.size(); iP++ )
                {
                    const uint tCol = tPassive( iP );
                    if ( aX( tCol ) <= 1e-14 * tBNorm || aX( tCol ) <= 0.0 )
                    {
                        aX( tCol )         = 0.0;
                        tInPassive( tCol ) = 0;
                    }
                    else
                    {
                        tNewPassive.push_back( tCol );
                    }
                }
                tPassive = tNewPassive;

                if ( tPassive.size() == 0 )
                {
                    break;
                }
            }

            // update residual r = b - A x (only passive columns contribute)
            for ( uint iR = 0; iR < tM; iR++ )
            {
                tRes( iR ) = aB( iR );
            }
            for ( uint iP = 0; iP < tPassive.size(); iP++ )
            {
                const uint tCol = tPassive( iP );
                const real tW   = aX( tCol );
                for ( uint iR = 0; iR < tM; iR++ )
                {
                    tRes( iR ) -= aA( iR, tCol ) * tW;
                }
            }

            // stagnation check
            real tCurResNorm = 0.0;
            for ( uint iR = 0; iR < tM; iR++ )
            {
                tCurResNorm += tRes( iR ) * tRes( iR );
            }
            tCurResNorm = std::sqrt( tCurResNorm );

            if ( tPrevResNorm - tCurResNorm <= 1e-15 * tBNorm )
            {
                break;
            }
            tPrevResNorm = tCurResNorm;
        }

        // final residual norm
        real tResNorm = 0.0;
        for ( uint iR = 0; iR < tM; iR++ )
        {
            tResNorm += tRes( iR ) * tRes( iR );
        }

        return std::sqrt( tResNorm );
    }

    //------------------------------------------------------------------------------

    void
    simplex_affine_map(
            const Matrix< DDRMat >& aCellCoords,
            Matrix< DDRMat >&       aM,
            Matrix< DDRMat >&       aOrigin,
            real&                   aDetJ )
    {
        const uint tDim = aCellCoords.n_cols();

        MORIS_ASSERT( ( aCellCoords.n_rows() == 4 && tDim == 3 ) ||    //
                              ( aCellCoords.n_rows() == 3 && tDim == 2 ),
                "moment_fitting::simplex_affine_map - cell coordinates must be 4 x 3 (TET) or 3 x 2 (TRI)." );

        // MORIS simplex parametrizations (vertex 2 carries the 1-sum shape function):
        // TET: N = [ z1, z2, 1 - z1 - z2 - z3, z3 ]
        //      xi( zeta ) = v2 + [ v0 - v2 | v1 - v2 | v3 - v2 ] * zeta
        // TRI: N = [ z1, z2, 1 - z1 - z2 ]
        //      xi( zeta ) = v2 + [ v0 - v2 | v1 - v2 ] * zeta
        aOrigin.set_size( tDim, 1 );
        aM.set_size( tDim, tDim );

        for ( uint iD = 0; iD < tDim; iD++ )
        {
            aOrigin( iD ) = aCellCoords( 2, iD );
            aM( iD, 0 )   = aCellCoords( 0, iD ) - aCellCoords( 2, iD );
            aM( iD, 1 )   = aCellCoords( 1, iD ) - aCellCoords( 2, iD );
            if ( tDim == 3 )
            {
                aM( iD, 2 ) = aCellCoords( 3, iD ) - aCellCoords( 2, iD );
            }
        }

        if ( tDim == 2 )
        {
            const real tDet = aM( 0, 0 ) * aM( 1, 1 ) - aM( 0, 1 ) * aM( 1, 0 );

            // MORIS tri space_det_J convention (reference area folded into the det)
            aDetJ = tDet / 2.0;
        }
        else
        {
            const real tDet =
                    +aM( 0, 0 ) * ( aM( 1, 1 ) * aM( 2, 2 ) - aM( 1, 2 ) * aM( 2, 1 ) )
                    - aM( 0, 1 ) * ( aM( 1, 0 ) * aM( 2, 2 ) - aM( 1, 2 ) * aM( 2, 0 ) )
                    + aM( 0, 2 ) * ( aM( 1, 0 ) * aM( 2, 1 ) - aM( 1, 1 ) * aM( 2, 0 ) );

            // MORIS tet space_det_J convention (reference volume folded into the det)
            aDetJ = tDet / 6.0;
        }
    }

    //------------------------------------------------------------------------------

    Fit_Result
    fit_cluster_rule( const Fit_Input& aInput )
    {
        Fit_Result tResult;

        const uint tNumCells = aInput.mCellParamCoords.size();
        const uint tNumCand  = aInput.mCandSpacePoints.n_cols();
        const uint tNumMom   = aInput.mMomSpacePoints.n_cols();

        if ( tNumCells == 0 || tNumCand == 0 || tNumMom == 0 )
        {
            return tResult;
        }

        const uint tDim = aInput.mCandSpacePoints.n_rows();
        const uint tM   = number_of_basis_functions( tDim, aInput.mDegree );

        // all candidate points in parent coordinates
        Matrix< DDRMat > tCandidates( tDim, tNumCells * tNumCand );

        Matrix< DDRMat > tMap, tOrigin, tMomPts( tDim, tNumMom ), tBasisMom;
        real             tDetJ = 0.0;

        for ( uint iCell = 0; iCell < tNumCells; iCell++ )
        {
            simplex_affine_map( aInput.mCellParamCoords( iCell ), tMap, tOrigin, tDetJ );

            // inverted or degenerate simplex in parametric space: refuse to fit
            if ( tDetJ <= 0.0 )
            {
                return tResult;
            }

            // map the candidate points
            for ( uint iP = 0; iP < tNumCand; iP++ )
            {
                for ( uint iD = 0; iD < tDim; iD++ )
                {
                    real tXi = tOrigin( iD );
                    for ( uint iE = 0; iE < tDim; iE++ )
                    {
                        tXi += tMap( iD, iE ) * aInput.mCandSpacePoints( iE, iP );
                    }
                    tCandidates( iD, iCell * tNumCand + iP ) = tXi;
                }
            }
        }

        tResult.mNumCandidates = tNumCells * tNumCand;

        // conditioning: evaluate the fitting basis in coordinates scaled to the
        // material bounding box (a thin sliver otherwise makes the basis columns
        // nearly collinear); the fit equations are equivalent - only the basis
        // of the fitted polynomial space changes
        Matrix< DDRMat > tBoxCenter( tDim, 1 );
        Matrix< DDRMat > tBoxHalf( tDim, 1 );

        for ( uint iD = 0; iD < tDim; iD++ )
        {
            real tLo = std::numeric_limits< real >::max();
            real tHi = std::numeric_limits< real >::lowest();

            for ( uint iCell = 0; iCell < tNumCells; iCell++ )
            {
                const uint tNv = aInput.mCellParamCoords( iCell ).n_rows();
                for ( uint iV = 0; iV < tNv; iV++ )
                {
                    tLo = std::min( tLo, aInput.mCellParamCoords( iCell )( iV, iD ) );
                    tHi = std::max( tHi, aInput.mCellParamCoords( iCell )( iV, iD ) );
                }
            }

            tBoxCenter( iD ) = 0.5 * ( tLo + tHi );
            tBoxHalf( iD )   = std::max( 0.5 * ( tHi - tLo ), 1e-12 );
        }

        auto tScalePoints = [ & ]( const Matrix< DDRMat >& aPoints, Matrix< DDRMat >& aScaled ) {
            aScaled.set_size( tDim, aPoints.n_cols() );
            for ( uint iP = 0; iP < aPoints.n_cols(); iP++ )
            {
                for ( uint iD = 0; iD < tDim; iD++ )
                {
                    aScaled( iD, iP ) = ( aPoints( iD, iP ) - tBoxCenter( iD ) ) / tBoxHalf( iD );
                }
            }
        };

        // moments of the material region in the scaled frame
        Matrix< DDRMat > tMoments( tM, 1, 0.0 );
        Matrix< DDRMat > tScaledMomPts;

        for ( uint iCell = 0; iCell < tNumCells; iCell++ )
        {
            simplex_affine_map( aInput.mCellParamCoords( iCell ), tMap, tOrigin, tDetJ );

            for ( uint iP = 0; iP < tNumMom; iP++ )
            {
                for ( uint iD = 0; iD < tDim; iD++ )
                {
                    real tXi = tOrigin( iD );
                    for ( uint iE = 0; iE < tDim; iE++ )
                    {
                        tXi += tMap( iD, iE ) * aInput.mMomSpacePoints( iE, iP );
                    }
                    tMomPts( iD, iP ) = tXi;
                }
            }

            tScalePoints( tMomPts, tScaledMomPts );
            evaluate_legendre_basis( tScaledMomPts, aInput.mDegree, tBasisMom );

            for ( uint iP = 0; iP < tNumMom; iP++ )
            {
                const real tW = aInput.mMomSpaceWeights( iP ) * tDetJ;
                for ( uint iB = 0; iB < tM; iB++ )
                {
                    tMoments( iB ) += tBasisMom( iB, iP ) * tW;
                }
            }
        }

        tResult.mMaterialVolume = tMoments( 0 );

        // vanishing material measure: keep the tessellated rule
        if ( !( tResult.mMaterialVolume > 0.0 ) )
        {
            return tResult;
        }

        real tMomNorm = 0.0;
        for ( uint iB = 0; iB < tM; iB++ )
        {
            tMomNorm += tMoments( iB ) * tMoments( iB );
        }
        tMomNorm = std::sqrt( tMomNorm );

        // fitting matrix at the candidates (scaled frame)
        Matrix< DDRMat > tScaledCandidates;
        tScalePoints( tCandidates, tScaledCandidates );

        Matrix< DDRMat > tA;
        evaluate_legendre_basis( tScaledCandidates, aInput.mDegree, tA );

        // non-negative fit; the active set is the pruning
        Matrix< DDRMat > tW;
        const real       tResNorm = nnls( tA, tMoments, tW );

        tResult.mResidual = tResNorm / tMomNorm;

        if ( tResult.mResidual > aInput.mRelTol )
        {
            return tResult;
        }

        // collect the retained points
        uint tNumRetained = 0;
        for ( uint iP = 0; iP < tW.numel(); iP++ )
        {
            if ( tW( iP ) > 0.0 )
            {
                tNumRetained++;
            }
        }

        // pruning failed to reduce below the tessellated count: keep the tessellated rule
        if ( tNumRetained == 0 || tNumRetained > tResult.mNumCandidates )
        {
            return tResult;
        }

        tResult.mPoints.set_size( tDim, tNumRetained );
        tResult.mWeights.set_size( 1, tNumRetained );

        uint tCount = 0;
        for ( uint iP = 0; iP < tW.numel(); iP++ )
        {
            if ( tW( iP ) > 0.0 )
            {
                for ( uint iD = 0; iD < tDim; iD++ )
                {
                    tResult.mPoints( iD, tCount ) = tCandidates( iD, iP );
                }
                tResult.mWeights( tCount ) = tW( iP );
                tCount++;
            }
        }

        tResult.mSuccess = true;

        return tResult;
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem::moment_fitting
