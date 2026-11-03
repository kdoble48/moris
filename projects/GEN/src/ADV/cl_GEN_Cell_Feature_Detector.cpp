/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Cell_Feature_Detector.cpp
 *
 */

#include "cl_GEN_Cell_Feature_Detector.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "cl_Vector.hpp"
#include "fn_det.hpp"
#include "fn_inv.hpp"

namespace moris::gen
{
    namespace
    {
        // parametric coordinate of corner aCorner (exodus order), direction aDim
        real
        corner_xi( uint aCorner, uint aDim )
        {
            // QUAD4 / bottom face of HEX8: (-1,-1), (1,-1), (1,1), (-1,1)
            const uint tInFace = aCorner % 4;
            if ( aDim == 0 )
            {
                return ( tInFace == 1 or tInFace == 2 ) ? 1.0 : -1.0;
            }
            if ( aDim == 1 )
            {
                return ( tInFace == 2 or tInFace == 3 ) ? 1.0 : -1.0;
            }
            return aCorner < 4 ? -1.0 : 1.0;
        }

        // multilinear corner shape function
        real
        shape_function( uint aCorner, const real* aXi, uint aNumDims )
        {
            real tValue = 1.0;
            for ( uint iDim = 0; iDim < aNumDims; iDim++ )
            {
                tValue *= 0.5 * ( 1.0 + corner_xi( aCorner, iDim ) * aXi[ iDim ] );
            }
            return tValue;
        }

        // one sample: parametric and global position, value, region (-1, 0, 1)
        struct Sample
        {
            real             mXi[ 3 ] = { 0.0, 0.0, 0.0 };
            Matrix< DDRMat > mX;
            real             mValue  = 0.0;
            signed char      mRegion = 0;
        };

        class Detector
        {
          public:
            Detector(
                    const Point_Field_Function&  aPhi,
                    const Matrix< DDRMat >&      aCorners,
                    real                         aTolerance,
                    const Cell_Feature_Settings& aSettings )
                    : mPhi( aPhi )
                    , mCorners( aCorners )
                    , mTolerance( aTolerance )
                    , mSettings( aSettings )
                    , mNumDims( aCorners.n_cols() )
            {
                MORIS_ERROR( mNumDims == 2 or mNumDims == 3, "detect_cell_features() - 2D or 3D cells only." );
                MORIS_ERROR( aCorners.n_rows() >= ( 1u << mNumDims ), "detect_cell_features() - need 2^d corner coordinates." );

                mM = std::max< uint >( 2, aSettings.mSamples + aSettings.mSamples % 2 );
                mN = mM + 1;
                mMaxCellBisections = aSettings.mMaxCellBisections > 0 ? aSettings.mMaxCellBisections : ( mNumDims == 2 ? 3 : 2 );
            }

            Cell_Feature_Report
            perform()
            {
                this->compute_measure();
                this->sample_lattice();
                this->check_crossings();
                if ( mSettings.mLipschitz > 0.0 )
                {
                    this->search_sub_cells();
                }
                this->check_edges();
                this->check_unseen();
                this->compute_template_values();
                if ( mSettings.mCheckTopology and mReport.mContainsCrossing )
                {
                    this->check_topology();
                }
                this->compute_measures();
                return mReport;
            }

          private:
            const Point_Field_Function&  mPhi;
            const Matrix< DDRMat >&      mCorners;
            real                         mTolerance;
            const Cell_Feature_Settings& mSettings;
            uint                         mNumDims;
            uint                         mM                 = 2;
            uint                         mN                 = 3;
            uint                         mMaxCellBisections = 3;
            std::vector< Sample >        mLattice;
            std::vector< real >          mTemplate;          // template value at each lattice point
            std::vector< Sample >        mFoundPoints;       // crossings found by the sub-cell search
            signed char                  mXTKRegion = 0;     // region of XTK's samples if uniform
            Cell_Feature_Report          mReport;

            signed char
            region( real aValue ) const
            {
                if ( std::abs( aValue ) <= mTolerance )
                {
                    return 0;
                }
                return aValue < 0.0 ? -1 : 1;
            }

            Matrix< DDRMat >
            map( const real* aXi ) const
            {
                Matrix< DDRMat > tX( 1, mNumDims, 0.0 );
                for ( uint iCorner = 0; iCorner < ( 1u << mNumDims ); iCorner++ )
                {
                    real tN = shape_function( iCorner, aXi, mNumDims );
                    for ( uint iDim = 0; iDim < mNumDims; iDim++ )
                    {
                        tX( iDim ) += tN * mCorners( iCorner, iDim );
                    }
                }
                return tX;
            }

            Sample
            evaluate( const real* aXi ) const
            {
                Sample tSample;
                for ( uint iDim = 0; iDim < mNumDims; iDim++ )
                {
                    tSample.mXi[ iDim ] = aXi[ iDim ];
                }
                tSample.mX      = this->map( aXi );
                tSample.mValue  = mPhi( tSample.mX );
                tSample.mRegion = this->region( tSample.mValue );
                return tSample;
            }

            uint
            index( const uint* aIJK ) const
            {
                uint tIndex = 0;
                for ( uint iDim = mNumDims; iDim-- > 0; )
                {
                    tIndex = tIndex * mN + aIJK[ iDim ];
                }
                return tIndex;
            }

            void
            unindex( uint aIndex, uint* aIJK ) const
            {
                for ( uint iDim = 0; iDim < mNumDims; iDim++ )
                {
                    aIJK[ iDim ] = aIndex % mN;
                    aIndex /= mN;
                }
            }

            uint
            num_lattice_points() const
            {
                return mNumDims == 2 ? mN * mN : mN * mN * mN;
            }

            void
            compute_measure()
            {
                // Jacobian of the multilinear map at the cell centre (exact for parallelograms / boxes)
                Matrix< DDRMat > tJacobian( mNumDims, mNumDims, 0.0 );
                real             tXi[ 3 ] = { 0.0, 0.0, 0.0 };
                for ( uint iCorner = 0; iCorner < ( 1u << mNumDims ); iCorner++ )
                {
                    for ( uint iXi = 0; iXi < mNumDims; iXi++ )
                    {
                        real tDerivative = 0.5 * corner_xi( iCorner, iXi );
                        for ( uint iOther = 0; iOther < mNumDims; iOther++ )
                        {
                            if ( iOther != iXi )
                            {
                                tDerivative *= 0.5 * ( 1.0 + corner_xi( iCorner, iOther ) * tXi[ iOther ] );
                            }
                        }
                        for ( uint iDim = 0; iDim < mNumDims; iDim++ )
                        {
                            tJacobian( iXi, iDim ) += tDerivative * mCorners( iCorner, iDim );
                        }
                    }
                }
                mReport.mCellMeasure = std::abs( det( tJacobian ) ) * std::pow( 2.0, mNumDims );
                mReport.mCellSize    = std::pow( mReport.mCellMeasure, 1.0 / mNumDims );
            }

            void
            sample_lattice()
            {
                uint tNumPoints = this->num_lattice_points();
                mLattice.resize( tNumPoints );
                uint tIJK[ 3 ] = { 0, 0, 0 };
                real tXi[ 3 ]  = { 0.0, 0.0, 0.0 };
                for ( uint iPoint = 0; iPoint < tNumPoints; iPoint++ )
                {
                    this->unindex( iPoint, tIJK );
                    for ( uint iDim = 0; iDim < mNumDims; iDim++ )
                    {
                        tXi[ iDim ] = -1.0 + 2.0 * real( tIJK[ iDim ] ) / real( mM );
                    }
                    mLattice[ iPoint ] = this->evaluate( tXi );
                }
            }

            // XTK's samples: corners + regular-subdivision vertices (QUAD4 centre; HEX8 face centres + centre), or all
            // 3^d nodes of a quadratic Lagrange cell
            bool
            is_xtk_sample( const uint* aIJK ) const
            {
                uint tNumOnBoundary = 0;
                for ( uint iDim = 0; iDim < mNumDims; iDim++ )
                {
                    if ( aIJK[ iDim ] == 0 or aIJK[ iDim ] == mM )
                    {
                        tNumOnBoundary++;
                    }
                    else if ( aIJK[ iDim ] != mM / 2 )
                    {
                        return false;
                    }
                }
                if ( mSettings.mLagrangeOrder >= 2 )
                {
                    return true;
                }
                // corner, cell centre, or (3D) face centre
                return tNumOnBoundary == mNumDims or tNumOnBoundary == 0 or ( mNumDims == 3 and tNumOnBoundary == 1 );
            }

            void
            check_crossings()
            {
                signed char tFirst      = mLattice[ 0 ].mRegion;
                bool        tXTKUniform = true;
                uint        tIJK[ 3 ]   = { 0, 0, 0 };
                for ( uint iPoint = 0; iPoint < mLattice.size(); iPoint++ )
                {
                    if ( mLattice[ iPoint ].mRegion != tFirst )
                    {
                        mReport.mContainsCrossing = true;
                        this->unindex( iPoint, tIJK );
                        if ( this->is_xtk_sample( tIJK ) )
                        {
                            tXTKUniform = false;
                        }
                    }
                }
                // corner 0 is an XTK sample, so the XTK samples are uniform iff none differs from it
                mReport.mXTKSeesCut = not tXTKUniform;
                mXTKRegion          = tFirst;
            }

            // geometric depth |phi| / |grad phi| of a sample, relative to the cell size. The slope per direction is the
            // larger one-sided difference, so that a kink of a min / max composition (the axis of a slot or a bridge),
            // where a central difference vanishes, does not inflate the depth.
            real
            relative_depth( const Sample& aSample ) const
            {
                real tStep        = mReport.mCellSize / ( 32.0 * mM );
                real tGradSquared = 0.0;
                for ( uint iDim = 0; iDim < mNumDims; iDim++ )
                {
                    Matrix< DDRMat > tPlus  = aSample.mX;
                    Matrix< DDRMat > tMinus = aSample.mX;
                    tPlus( iDim ) += tStep;
                    tMinus( iDim ) -= tStep;
                    real tForward    = std::abs( mPhi( tPlus ) - aSample.mValue ) / tStep;
                    real tBackward   = std::abs( aSample.mValue - mPhi( tMinus ) ) / tStep;
                    real tDerivative = std::max( tForward, tBackward );
                    tGradSquared += tDerivative * tDerivative;
                }
                if ( tGradSquared <= 0.0 )
                {
                    return std::numeric_limits< real >::max();
                }
                return std::abs( aSample.mValue ) / std::sqrt( tGradSquared ) / mReport.mCellSize;
            }

            // record a hidden feature of the given relative depth; returns whether it is significant
            bool
            record_hidden( real aRelativeDepth )
            {
                if ( aRelativeDepth >= mSettings.mDepthTolerance )
                {
                    mReport.mMaxHiddenDepth = std::max( mReport.mMaxHiddenDepth, std::min( aRelativeDepth, 1.0e3 ) );
                    return true;
                }
                mReport.mNumInsignificant++;
                return false;
            }

            real
            physical_distance( const Sample& aA, const Sample& aB ) const
            {
                real tSquared = 0.0;
                for ( uint iDim = 0; iDim < mNumDims; iDim++ )
                {
                    tSquared += ( aA.mX( iDim ) - aB.mX( iDim ) ) * ( aA.mX( iDim ) - aB.mX( iDim ) );
                }
                return std::sqrt( tSquared );
            }

            // Lipschitz bisection of the segment [aA, aB] (same region, not interface); appends inner samples in order
            void
            bisect_segment( const Sample& aA, const Sample& aB, uint aDepth, std::vector< Sample >& aOut )
            {
                if ( aA.mRegion != aB.mRegion or aA.mRegion == 0 )
                {
                    return;
                }
                if ( std::abs( aA.mValue ) + std::abs( aB.mValue ) > mSettings.mLipschitz * this->physical_distance( aA, aB ) )
                {
                    return;    // certified: no zero between a and b
                }
                if ( aDepth >= mSettings.mMaxEdgeBisections )
                {
                    mReport.mNumUncertified++;
                    return;
                }
                real tXi[ 3 ] = { 0.0, 0.0, 0.0 };
                for ( uint iDim = 0; iDim < mNumDims; iDim++ )
                {
                    tXi[ iDim ] = 0.5 * ( aA.mXi[ iDim ] + aB.mXi[ iDim ] );
                }
                Sample tMid = this->evaluate( tXi );
                this->bisect_segment( aA, tMid, aDepth + 1, aOut );
                aOut.push_back( tMid );
                this->bisect_segment( tMid, aB, aDepth + 1, aOut );
            }

            void
            check_edges()
            {
                const Vector< std::pair< uint, uint > >& tEdges = cell_feature_edges( mNumDims );
                for ( uint iEdge = 0; iEdge < tEdges.size(); iEdge++ )
                {
                    // lattice line from corner a to corner b
                    uint tIJKA[ 3 ] = { 0, 0, 0 };
                    uint tIJKB[ 3 ] = { 0, 0, 0 };
                    for ( uint iDim = 0; iDim < mNumDims; iDim++ )
                    {
                        tIJKA[ iDim ] = corner_xi( tEdges( iEdge ).first, iDim ) > 0.0 ? mM : 0;
                        tIJKB[ iDim ] = corner_xi( tEdges( iEdge ).second, iDim ) > 0.0 ? mM : 0;
                    }
                    std::vector< Sample > tLine;
                    for ( uint iStep = 0; iStep <= mM; iStep++ )
                    {
                        uint tIJK[ 3 ] = { 0, 0, 0 };
                        for ( uint iDim = 0; iDim < mNumDims; iDim++ )
                        {
                            tIJK[ iDim ] = ( tIJKA[ iDim ] * ( mM - iStep ) + tIJKB[ iDim ] * iStep ) / mM;
                        }
                        const Sample& tSample = mLattice[ this->index( tIJK ) ];
                        if ( iStep > 0 and mSettings.mLipschitz > 0.0 )
                        {
                            std::vector< Sample > tInner;
                            this->bisect_segment( tLine.back(), tSample, 0, tInner );
                            for ( const Sample& iInner : tInner )
                            {
                                if ( iInner.mRegion != tLine.back().mRegion )
                                {
                                    mReport.mContainsCrossing = true;
                                }
                                tLine.push_back( iInner );
                            }
                        }
                        tLine.push_back( tSample );
                    }

                    // runs of equal non-interface regions along the edge
                    std::vector< signed char >          tRunRegion;
                    std::vector< std::vector< uint > >  tRunSamples;
                    for ( uint iSample = 0; iSample < tLine.size(); iSample++ )
                    {
                        signed char tRegion = tLine[ iSample ].mRegion;
                        if ( tRegion == 0 )
                        {
                            continue;
                        }
                        if ( tRunRegion.empty() or tRunRegion.back() != tRegion )
                        {
                            tRunRegion.push_back( tRegion );
                            tRunSamples.emplace_back();
                        }
                        tRunSamples.back().push_back( iSample );
                    }

                    // runs the cut misses: the cut puts at most one crossing on the edge, between its end regions
                    signed char tEndA = tLine.front().mRegion;
                    signed char tEndB = tLine.back().mRegion;
                    std::vector< uint > tHiddenRuns;
                    for ( uint iRun = 0; iRun < tRunRegion.size(); iRun++ )
                    {
                        bool tHidden = false;
                        if ( tEndA != 0 and tEndB != 0 )
                        {
                            tHidden = iRun > 0 and iRun + 1 < tRunRegion.size();
                        }
                        else if ( tEndA == 0 and tEndB != 0 )
                        {
                            tHidden = tRunRegion[ iRun ] != tEndB;
                        }
                        else if ( tEndB == 0 and tEndA != 0 )
                        {
                            tHidden = tRunRegion[ iRun ] != tEndA;
                        }
                        else
                        {
                            tHidden = tRunRegion.size() > 1 and tRunRegion[ iRun ] != tRunRegion[ 0 ];
                        }
                        if ( tHidden )
                        {
                            tHiddenRuns.push_back( iRun );
                        }
                    }
                    if ( tHiddenRuns.empty() )
                    {
                        continue;
                    }

                    // depth of the deepest hidden run
                    real tDepth = 0.0;
                    for ( uint iRun : tHiddenRuns )
                    {
                        for ( uint iSample : tRunSamples[ iRun ] )
                        {
                            tDepth = std::max( tDepth, this->relative_depth( tLine[ iSample ] ) );
                        }
                    }
                    if ( this->record_hidden( tDepth ) )
                    {
                        mReport.mHiddenEdgeMask |= ( 1u << iEdge );
                    }
                }
            }

            // Lipschitz search of the lattice sub-cells with uniform corner regions
            void
            search_box( const real* aLower, const real* aUpper, signed char aRegion, uint aDepth )
            {
                real tCentre[ 3 ] = { 0.0, 0.0, 0.0 };
                real tHalfDiagonal = 0.0;
                for ( uint iDim = 0; iDim < mNumDims; iDim++ )
                {
                    tCentre[ iDim ] = 0.5 * ( aLower[ iDim ] + aUpper[ iDim ] );
                    real tHalfExtent = 0.25 * ( aUpper[ iDim ] - aLower[ iDim ] ) * mReport.mCellSize;
                    tHalfDiagonal += tHalfExtent * tHalfExtent;
                }
                tHalfDiagonal = std::sqrt( tHalfDiagonal );
                Sample tSample = this->evaluate( tCentre );
                if ( tSample.mRegion != aRegion )
                {
                    mFoundPoints.push_back( tSample );
                    return;
                }
                if ( std::abs( tSample.mValue ) > mSettings.mLipschitz * tHalfDiagonal )
                {
                    return;    // certified
                }
                if ( aDepth >= mMaxCellBisections )
                {
                    mReport.mNumUncertified++;
                    return;
                }
                for ( uint iChild = 0; iChild < ( 1u << mNumDims ); iChild++ )
                {
                    real tLower[ 3 ] = { 0.0, 0.0, 0.0 };
                    real tUpper[ 3 ] = { 0.0, 0.0, 0.0 };
                    for ( uint iDim = 0; iDim < mNumDims; iDim++ )
                    {
                        bool tHigh      = ( iChild >> iDim ) & 1;
                        tLower[ iDim ] = tHigh ? tCentre[ iDim ] : aLower[ iDim ];
                        tUpper[ iDim ] = tHigh ? aUpper[ iDim ] : tCentre[ iDim ];
                    }
                    this->search_box( tLower, tUpper, aRegion, aDepth + 1 );
                }
            }

            void
            search_sub_cells()
            {
                uint tNumSub = mNumDims == 2 ? mM * mM : mM * mM * mM;
                for ( uint iSub = 0; iSub < tNumSub; iSub++ )
                {
                    uint tIJK[ 3 ] = { 0, 0, 0 };
                    uint tRest     = iSub;
                    for ( uint iDim = 0; iDim < mNumDims; iDim++ )
                    {
                        tIJK[ iDim ] = tRest % mM;
                        tRest /= mM;
                    }
                    // corner regions of the sub-cell
                    signed char tRegion  = 0;
                    bool        tUniform = true;
                    for ( uint iCorner = 0; iCorner < ( 1u << mNumDims ) and tUniform; iCorner++ )
                    {
                        uint tCornerIJK[ 3 ] = { 0, 0, 0 };
                        for ( uint iDim = 0; iDim < mNumDims; iDim++ )
                        {
                            tCornerIJK[ iDim ] = tIJK[ iDim ] + ( ( iCorner >> iDim ) & 1 );
                        }
                        signed char tCornerRegion = mLattice[ this->index( tCornerIJK ) ].mRegion;
                        if ( iCorner == 0 )
                        {
                            tRegion = tCornerRegion;
                        }
                        tUniform = tUniform and tCornerRegion == tRegion;
                    }
                    if ( not tUniform or tRegion == 0 )
                    {
                        continue;
                    }
                    real tLower[ 3 ] = { 0.0, 0.0, 0.0 };
                    real tUpper[ 3 ] = { 0.0, 0.0, 0.0 };
                    for ( uint iDim = 0; iDim < mNumDims; iDim++ )
                    {
                        tLower[ iDim ] = -1.0 + 2.0 * real( tIJK[ iDim ] ) / real( mM );
                        tUpper[ iDim ] = -1.0 + 2.0 * real( tIJK[ iDim ] + 1 ) / real( mM );
                    }
                    this->search_box( tLower, tUpper, tRegion, 0 );
                }
                if ( not mFoundPoints.empty() )
                {
                    mReport.mContainsCrossing = true;
                }
            }

            void
            check_unseen()
            {
                if ( not mReport.mContainsCrossing or mReport.mXTKSeesCut )
                {
                    return;
                }
                // the cell is cut but XTK's samples all share mXTKRegion
                real tDepth = 0.0;
                bool tAny   = false;
                for ( const Sample& iSample : mLattice )
                {
                    if ( iSample.mRegion != 0 and iSample.mRegion != mXTKRegion )
                    {
                        tDepth = std::max( tDepth, this->relative_depth( iSample ) );
                        tAny   = true;
                    }
                }
                for ( const Sample& iSample : mFoundPoints )
                {
                    if ( iSample.mRegion != 0 and iSample.mRegion != mXTKRegion )
                    {
                        tDepth = std::max( tDepth, this->relative_depth( iSample ) );
                        tAny   = true;
                    }
                }
                if ( tAny )
                {
                    mReport.mUnseenCrossing = this->record_hidden( tDepth );
                }
            }

            // value of the piecewise-linear template at a parametric point: the P1 interpolant of the XTK samples on
            // the regular subdivision (4 triangles around the centre in 2D, 24 tets around face and cell centres in 3D)
            real
            template_value( const real* aXi ) const
            {
                uint tCentreIJK[ 3 ] = { mM / 2, mM / 2, mM / 2 };
                real tCentre         = mLattice[ this->index( tCentreIJK ) ].mValue;

                auto tValueAt = [ & ]( const real* aPointXi ) -> real {
                    uint tIJK[ 3 ] = { 0, 0, 0 };
                    for ( uint iDim = 0; iDim < mNumDims; iDim++ )
                    {
                        tIJK[ iDim ] = uint( std::lround( 0.5 * ( aPointXi[ iDim ] + 1.0 ) * mM ) );
                    }
                    return mLattice[ this->index( tIJK ) ].mValue;
                };

                if ( mNumDims == 2 )
                {
                    // triangle (A, B, centre): the side the point faces
                    real tA[ 2 ];
                    real tB[ 2 ];
                    if ( std::abs( aXi[ 0 ] ) >= std::abs( aXi[ 1 ] ) )
                    {
                        real tS = aXi[ 0 ] >= 0.0 ? 1.0 : -1.0;
                        tA[ 0 ] = tS, tA[ 1 ] = -1.0, tB[ 0 ] = tS, tB[ 1 ] = 1.0;
                    }
                    else
                    {
                        real tS = aXi[ 1 ] >= 0.0 ? 1.0 : -1.0;
                        tA[ 0 ] = -1.0, tA[ 1 ] = tS, tB[ 0 ] = 1.0, tB[ 1 ] = tS;
                    }
                    // aXi = a A + b B (centre at the origin)
                    real tDet = tA[ 0 ] * tB[ 1 ] - tA[ 1 ] * tB[ 0 ];
                    real tWa  = ( aXi[ 0 ] * tB[ 1 ] - aXi[ 1 ] * tB[ 0 ] ) / tDet;
                    real tWb  = ( tA[ 0 ] * aXi[ 1 ] - tA[ 1 ] * aXi[ 0 ] ) / tDet;
                    return tCentre + tWa * ( tValueAt( tA ) - tCentre ) + tWb * ( tValueAt( tB ) - tCentre );
                }

                // 3D: face of the dominant direction, then the edge of that face the point faces
                uint tK = 0;
                for ( uint iDim = 1; iDim < 3; iDim++ )
                {
                    if ( std::abs( aXi[ iDim ] ) > std::abs( aXi[ tK ] ) )
                    {
                        tK = iDim;
                    }
                }
                real tS  = aXi[ tK ] >= 0.0 ? 1.0 : -1.0;
                uint tI  = ( tK + 1 ) % 3;
                uint tJ  = ( tK + 2 ) % 3;
                real tF[ 3 ] = { 0.0, 0.0, 0.0 };
                tF[ tK ]     = tS;
                real tA[ 3 ] = { 0.0, 0.0, 0.0 };
                real tB[ 3 ] = { 0.0, 0.0, 0.0 };
                tA[ tK ] = tS, tB[ tK ] = tS;
                if ( std::abs( aXi[ tI ] ) >= std::abs( aXi[ tJ ] ) )
                {
                    real tU = aXi[ tI ] >= 0.0 ? 1.0 : -1.0;
                    tA[ tI ] = tU, tA[ tJ ] = -1.0, tB[ tI ] = tU, tB[ tJ ] = 1.0;
                }
                else
                {
                    real tV = aXi[ tJ ] >= 0.0 ? 1.0 : -1.0;
                    tA[ tJ ] = tV, tA[ tI ] = -1.0, tB[ tJ ] = tV, tB[ tI ] = 1.0;
                }
                // aXi = f F + a A + b B (cell centre at the origin)
                Matrix< DDRMat > tMatrix( 3, 3 );
                for ( uint iDim = 0; iDim < 3; iDim++ )
                {
                    tMatrix( iDim, 0 ) = tF[ iDim ];
                    tMatrix( iDim, 1 ) = tA[ iDim ];
                    tMatrix( iDim, 2 ) = tB[ iDim ];
                }
                Matrix< DDRMat > tRhs( 3, 1 );
                for ( uint iDim = 0; iDim < 3; iDim++ )
                {
                    tRhs( iDim ) = aXi[ iDim ];
                }
                Matrix< DDRMat > tWeights = inv( tMatrix ) * tRhs;
                return tCentre + tWeights( 0 ) * ( tValueAt( tF ) - tCentre ) + tWeights( 1 ) * ( tValueAt( tA ) - tCentre )
                     + tWeights( 2 ) * ( tValueAt( tB ) - tCentre );
            }

            void
            compute_template_values()
            {
                mTemplate.resize( mLattice.size() );
                for ( uint iPoint = 0; iPoint < mLattice.size(); iPoint++ )
                {
                    mTemplate[ iPoint ] = this->template_value( mLattice[ iPoint ].mXi );
                }
            }

            // number of connected components of the lattice points with aMask == aValue; face connectivity for the
            // positive side, full (face, edge, corner) connectivity for the negative side, as in digital topology
            uint
            count_components( const std::vector< char >& aMask, char aValue, bool aFullConnectivity ) const
            {
                std::vector< char > tVisited( aMask.size(), 0 );
                std::vector< uint > tStack;
                uint                tCount = 0;
                for ( uint iSeed = 0; iSeed < aMask.size(); iSeed++ )
                {
                    if ( aMask[ iSeed ] != aValue or tVisited[ iSeed ] )
                    {
                        continue;
                    }
                    tCount++;
                    tVisited[ iSeed ] = 1;
                    tStack.push_back( iSeed );
                    while ( not tStack.empty() )
                    {
                        uint tPoint = tStack.back();
                        tStack.pop_back();
                        uint tIJK[ 3 ] = { 0, 0, 0 };
                        this->unindex( tPoint, tIJK );
                        int tNumNeighbors = mNumDims == 2 ? 9 : 27;
                        for ( int iNeighbor = 0; iNeighbor < tNumNeighbors; iNeighbor++ )
                        {
                            int  tOffset[ 3 ] = { iNeighbor % 3 - 1, ( iNeighbor / 3 ) % 3 - 1, iNeighbor / 9 - 1 };
                            int  tNonZero     = 0;
                            bool tInside      = true;
                            uint tOther[ 3 ]  = { 0, 0, 0 };
                            for ( uint iDim = 0; iDim < mNumDims; iDim++ )
                            {
                                tNonZero += tOffset[ iDim ] != 0;
                                int tCoordinate = int( tIJK[ iDim ] ) + tOffset[ iDim ];
                                tInside         = tInside and tCoordinate >= 0 and tCoordinate < int( mN );
                                tOther[ iDim ]  = uint( std::max( tCoordinate, 0 ) );
                            }
                            if ( tNonZero == 0 or not tInside or ( not aFullConnectivity and tNonZero > 1 ) )
                            {
                                continue;
                            }
                            uint tOtherIndex = this->index( tOther );
                            if ( aMask[ tOtherIndex ] == aValue and not tVisited[ tOtherIndex ] )
                            {
                                tVisited[ tOtherIndex ] = 1;
                                tStack.push_back( tOtherIndex );
                            }
                        }
                    }
                }
                return tCount;
            }

            void
            check_topology()
            {
                std::vector< char > tTrue( mLattice.size() );
                std::vector< char > tModel( mLattice.size() );
                for ( uint iPoint = 0; iPoint < mLattice.size(); iPoint++ )
                {
                    tTrue[ iPoint ]  = mLattice[ iPoint ].mValue >= -mTolerance;
                    tModel[ iPoint ] = mTemplate[ iPoint ] >= -mTolerance;
                }
                bool tMismatch = this->count_components( tTrue, 1, false ) != this->count_components( tModel, 1, false )
                              or this->count_components( tTrue, 0, true ) != this->count_components( tModel, 0, true );

                // depth of the samples whose side differs from the template's
                real tDepth = 0.0;
                if ( tMismatch )
                {
                    for ( uint iPoint = 0; iPoint < mLattice.size(); iPoint++ )
                    {
                        if ( tTrue[ iPoint ] != tModel[ iPoint ] and mLattice[ iPoint ].mRegion != 0 )
                        {
                            tDepth = std::max( tDepth, this->relative_depth( mLattice[ iPoint ] ) );
                        }
                    }
                }

                // crossings found inside sub-cells by the certificate change the components wherever the template has
                // the other side there
                for ( const Sample& iSample : mFoundPoints )
                {
                    if ( iSample.mRegion != 0 and ( iSample.mValue >= -mTolerance ) != ( this->template_value( iSample.mXi ) >= -mTolerance ) )
                    {
                        tMismatch = true;
                        tDepth    = std::max( tDepth, this->relative_depth( iSample ) );
                    }
                }

                if ( tMismatch )
                {
                    mReport.mTopologyMismatch = this->record_hidden( tDepth );
                }
            }

            void
            compute_measures()
            {
                // trapezoidal weights on the lattice
                real tTotal    = 0.0;
                real tSampled  = 0.0;
                real tTemplate = 0.0;
                uint tIJK[ 3 ] = { 0, 0, 0 };
                for ( uint iPoint = 0; iPoint < mLattice.size(); iPoint++ )
                {
                    this->unindex( iPoint, tIJK );
                    real tWeight = 1.0;
                    for ( uint iDim = 0; iDim < mNumDims; iDim++ )
                    {
                        tWeight *= ( tIJK[ iDim ] == 0 or tIJK[ iDim ] == mM ) ? 0.5 : 1.0;
                    }
                    tTotal += tWeight;
                    tSampled += mLattice[ iPoint ].mValue > 0.0 ? tWeight : 0.0;
                    tTemplate += mTemplate[ iPoint ] > 0.0 ? tWeight : 0.0;
                }
                mReport.mSampledPositive  = mReport.mCellMeasure * tSampled / tTotal;
                mReport.mTemplatePositive = mReport.mCellMeasure * tTemplate / tTotal;
            }
        };
    }    // namespace

    //--------------------------------------------------------------------------------------------------------------

    uint
    Cell_Feature_Report::num_hidden_edges() const
    {
        uint tCount = 0;
        for ( uint tMask = mHiddenEdgeMask; tMask != 0; tMask >>= 1 )
        {
            tCount += tMask & 1;
        }
        return tCount;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Vector< std::pair< uint, uint > >&
    cell_feature_edges( uint aNumDimensions )
    {
        static const Vector< std::pair< uint, uint > > tQuadEdges = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 } };
        static const Vector< std::pair< uint, uint > > tHexEdges  = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 }, { 4, 5 }, { 5, 6 },
             { 6, 7 }, { 7, 4 }, { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 } };
        return aNumDimensions == 2 ? tQuadEdges : tHexEdges;
    }

    //--------------------------------------------------------------------------------------------------------------

    Cell_Feature_Report
    detect_cell_features(
            const Point_Field_Function&  aPhi,
            const Matrix< DDRMat >&      aCorners,
            real                         aInterfaceTolerance,
            const Cell_Feature_Settings& aSettings )
    {
        Detector tDetector( aPhi, aCorners, aInterfaceTolerance, aSettings );
        return tDetector.perform();
    }

    //--------------------------------------------------------------------------------------------------------------

    uint
    cell_feature_flags( const Cell_Feature_Report& aReport )
    {
        uint tFlags = 0;
        tFlags |= aReport.mContainsCrossing ? CELL_FEATURE_CROSSING : 0;
        tFlags |= aReport.hides_feature() ? CELL_FEATURE_HIDDEN : 0;
        tFlags |= aReport.mXTKSeesCut ? CELL_FEATURE_XTK_CUT : 0;
        tFlags |= aReport.mNumInsignificant > 0 ? CELL_FEATURE_INSIGNIFICANT : 0;
        return tFlags;
    }

    //--------------------------------------------------------------------------------------------------------------
    // Geometry-aware query points
    //--------------------------------------------------------------------------------------------------------------

    namespace
    {
        std::string
        trim_blanks( const std::string& aString )
        {
            std::string tOut;
            for ( char iChar : aString )
            {
                if ( not std::isspace( static_cast< unsigned char >( iChar ) ) )
                {
                    tOut.push_back( iChar );
                }
            }
            return tOut;
        }

        std::vector< std::string >
        split( const std::string& aString, char aSeparator )
        {
            std::vector< std::string > tParts;
            std::string                tPart;
            for ( char iChar : aString )
            {
                if ( iChar == aSeparator )
                {
                    tParts.push_back( tPart );
                    tPart.clear();
                }
                else
                {
                    tPart.push_back( iChar );
                }
            }
            tParts.push_back( tPart );
            return tParts;
        }

        real
        to_real( const std::string& aString, const std::string& aItem )
        {
            std::size_t tUsed  = 0;
            real        tValue = 0.0;
            try
            {
                tValue = std::stod( aString, &tUsed );
            } catch ( ... )
            {
                tUsed = 0;
            }
            MORIS_ERROR( tUsed > 0 and tUsed == aString.size(), "feature_query_points: cannot read '%s' in '%s'.", aString.c_str(), aItem.c_str() );
            return tValue;
        }

        void
        check_query_point( const Feature_Query_Point& aPoint, const std::string& aItem )
        {
            MORIS_ERROR( aPoint.mCoordinates.numel() == 2 or aPoint.mCoordinates.numel() == 3,
                    "feature_query_points: point '%s' needs 2 or 3 coordinates.",
                    aItem.c_str() );
            MORIS_ERROR( aPoint.mPhase >= -1 and aPoint.mPhase <= 1,
                    "feature_query_points: phase of point '%s' must be +1, -1 or 0.",
                    aItem.c_str() );
            MORIS_ERROR( aPoint.mRadius >= 0.0, "feature_query_points: radius of point '%s' must not be negative.", aItem.c_str() );
        }

        // determinant of the leading d x d block of aRows (d = 2, 3)
        real
        small_det( const real aRows[ 3 ][ 3 ], uint aDim )
        {
            if ( aDim == 2 )
            {
                return aRows[ 0 ][ 0 ] * aRows[ 1 ][ 1 ] - aRows[ 0 ][ 1 ] * aRows[ 1 ][ 0 ];
            }
            return aRows[ 0 ][ 0 ] * ( aRows[ 1 ][ 1 ] * aRows[ 2 ][ 2 ] - aRows[ 1 ][ 2 ] * aRows[ 2 ][ 1 ] )
                 - aRows[ 0 ][ 1 ] * ( aRows[ 1 ][ 0 ] * aRows[ 2 ][ 2 ] - aRows[ 1 ][ 2 ] * aRows[ 2 ][ 0 ] )
                 + aRows[ 0 ][ 2 ] * ( aRows[ 1 ][ 0 ] * aRows[ 2 ][ 1 ] - aRows[ 1 ][ 1 ] * aRows[ 2 ][ 0 ] );
        }

        using Point3 = std::array< real, 3 >;

        // whether aPoint lies in the simplex aVertices (d + 1 vertices), barycentric coordinates by Cramer's rule
        bool
        simplex_contains( const std::vector< Point3 >& aVertices, const Point3& aPoint, uint aDim )
        {
            real tEdges[ 3 ][ 3 ] = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } };
            for ( uint iEdge = 0; iEdge < aDim; iEdge++ )
            {
                for ( uint iDim = 0; iDim < aDim; iDim++ )
                {
                    tEdges[ iDim ][ iEdge ] = aVertices[ iEdge + 1 ][ iDim ] - aVertices[ 0 ][ iDim ];
                }
            }
            real tDet = small_det( tEdges, aDim );
            if ( std::abs( tDet ) <= std::numeric_limits< real >::min() )
            {
                return false;    // collapsed piece
            }
            real tSum = 0.0;
            for ( uint iEdge = 0; iEdge < aDim; iEdge++ )
            {
                real tReplaced[ 3 ][ 3 ];
                for ( uint iRow = 0; iRow < 3; iRow++ )
                {
                    for ( uint iCol = 0; iCol < 3; iCol++ )
                    {
                        tReplaced[ iRow ][ iCol ] = tEdges[ iRow ][ iCol ];
                    }
                }
                for ( uint iDim = 0; iDim < aDim; iDim++ )
                {
                    tReplaced[ iDim ][ iEdge ] = aPoint[ iDim ] - aVertices[ 0 ][ iDim ];
                }
                real tLambda = small_det( tReplaced, aDim ) / tDet;
                if ( tLambda < -1.0e-10 )
                {
                    return false;
                }
                tSum += tLambda;
            }
            return tSum <= 1.0 + 1.0e-10;
        }

        struct Template_Vertex
        {
            Point3      mX      = { 0.0, 0.0, 0.0 };
            real        mValue  = 0.0;
            signed char mRegion = 0;
        };

        // side of the level set that the cut of one template simplex represents at a point
        class Represented_Region
        {
          public:
            Represented_Region(
                    const Point_Field_Function& aPhi,
                    uint                        aDim,
                    real                        aTolerance,
                    bool                        aExactRoots )
                    : mPhi( aPhi )
                    , mDim( aDim )
                    , mTolerance( aTolerance )
                    , mExactRoots( aExactRoots )
            {
            }

            Template_Vertex
            vertex( const Point3& aX ) const
            {
                Template_Vertex tVertex;
                tVertex.mX      = aX;
                tVertex.mValue  = this->value( aX );
                tVertex.mRegion = std::abs( tVertex.mValue ) <= mTolerance ? 0 : ( tVertex.mValue < 0.0 ? -1 : 1 );
                return tVertex;
            }

            sint
            region_in_simplex( const std::vector< Template_Vertex >& aSimplex, const Point3& aPoint ) const
            {
                std::vector< uint > tPositive;
                std::vector< uint > tNegative;
                std::vector< uint > tZero;
                for ( uint iVertex = 0; iVertex < aSimplex.size(); iVertex++ )
                {
                    signed char tRegion = aSimplex[ iVertex ].mRegion;
                    if ( tRegion > 0 )
                    {
                        tPositive.push_back( iVertex );
                    }
                    else if ( tRegion < 0 )
                    {
                        tNegative.push_back( iVertex );
                    }
                    else
                    {
                        tZero.push_back( iVertex );
                    }
                }
                if ( tPositive.empty() and tNegative.empty() )
                {
                    return 0;
                }
                if ( tNegative.empty() )
                {
                    return 1;
                }
                if ( tPositive.empty() )
                {
                    return -1;
                }

                // the side with fewer vertices
                bool                       tLonePositive = tPositive.size() <= tNegative.size();
                const std::vector< uint >& tLone         = tLonePositive ? tPositive : tNegative;
                const std::vector< uint >& tOther        = tLonePositive ? tNegative : tPositive;
                sint                       tLoneRegion   = tLonePositive ? 1 : -1;

                if ( tLone.size() == 1 )
                {
                    // the lone vertex, the intersections of its edges to the other side, and the interface vertices
                    std::vector< Point3 > tPiece = { aSimplex[ tLone[ 0 ] ].mX };
                    for ( uint iOther : tOther )
                    {
                        tPiece.push_back( this->root( aSimplex[ tLone[ 0 ] ], aSimplex[ iOther ] ) );
                    }
                    for ( uint iZero : tZero )
                    {
                        tPiece.push_back( aSimplex[ iZero ].mX );
                    }
                    return simplex_contains( tPiece, aPoint, mDim ) ? tLoneRegion : -tLoneRegion;
                }

                // 3D, two vertices on each side: the prism (a0, a1, a2 | b0, b1, b2) of the positive side in three tets
                const Template_Vertex&               tA    = aSimplex[ tPositive[ 0 ] ];
                const Template_Vertex&               tB    = aSimplex[ tPositive[ 1 ] ];
                const Template_Vertex&               tN0   = aSimplex[ tNegative[ 0 ] ];
                const Template_Vertex&               tN1   = aSimplex[ tNegative[ 1 ] ];
                Point3                               tA0   = tA.mX;
                Point3                               tA1   = this->root( tA, tN0 );
                Point3                               tA2   = this->root( tA, tN1 );
                Point3                               tB0   = tB.mX;
                Point3                               tB1   = this->root( tB, tN0 );
                Point3                               tB2   = this->root( tB, tN1 );
                std::vector< std::vector< Point3 > > tTets = { { tA0, tA1, tA2, tB0 }, { tA1, tA2, tB0, tB1 }, { tA2, tB0, tB1, tB2 } };
                for ( const std::vector< Point3 >& iTet : tTets )
                {
                    if ( simplex_contains( iTet, aPoint, mDim ) )
                    {
                        return 1;
                    }
                }
                return -1;
            }

          private:
            const Point_Field_Function& mPhi;
            uint                        mDim;
            real                        mTolerance;
            bool                        mExactRoots;

            real
            value( const Point3& aX ) const
            {
                Matrix< DDRMat > tX( 1, mDim );
                for ( uint iDim = 0; iDim < mDim; iDim++ )
                {
                    tX( iDim ) = aX[ iDim ];
                }
                return mPhi( tX );
            }

            Point3
            along( const Point3& aA, const Point3& aB, real aT ) const
            {
                Point3 tX = { 0.0, 0.0, 0.0 };
                for ( uint iDim = 0; iDim < mDim; iDim++ )
                {
                    tX[ iDim ] = aA[ iDim ] + aT * ( aB[ iDim ] - aA[ iDim ] );
                }
                return tX;
            }

            // intersection on the edge between two vertices on opposite sides: exact root (bisection) or linear crossing
            Point3
            root( const Template_Vertex& aA, const Template_Vertex& aB ) const
            {
                real tT = aA.mValue / ( aA.mValue - aB.mValue );
                if ( mExactRoots )
                {
                    real tLow           = 0.0;
                    real tHigh          = 1.0;
                    bool tLowIsNegative = aA.mValue < 0.0;
                    for ( uint iIteration = 0; iIteration < 100 and tHigh - tLow > 1.0e-15; iIteration++ )
                    {
                        real tMid      = 0.5 * ( tLow + tHigh );
                        real tMidValue = this->value( this->along( aA.mX, aB.mX, tMid ) );
                        if ( ( tMidValue < 0.0 ) == tLowIsNegative )
                        {
                            tLow = tMid;
                        }
                        else
                        {
                            tHigh = tMid;
                        }
                    }
                    tT = 0.5 * ( tLow + tHigh );
                }
                return this->along( aA.mX, aB.mX, tT );
            }
        };

        // lower / upper bound of the corners in one direction
        void
        corner_bounds( const Matrix< DDRMat >& aCorners, uint aDim, real& aLower, real& aUpper )
        {
            aLower = std::numeric_limits< real >::max();
            aUpper = -std::numeric_limits< real >::max();
            for ( uint iCorner = 0; iCorner < aCorners.n_rows(); iCorner++ )
            {
                aLower = std::min( aLower, aCorners( iCorner, aDim ) );
                aUpper = std::max( aUpper, aCorners( iCorner, aDim ) );
            }
        }
    }    // namespace

    //--------------------------------------------------------------------------------------------------------------

    Vector< Feature_Query_Point >
    parse_feature_query_points( const std::string& aString )
    {
        Vector< Feature_Query_Point > tPoints;
        for ( const std::string& iItem : split( trim_blanks( aString ), ';' ) )
        {
            if ( iItem.empty() )
            {
                continue;
            }
            std::vector< std::string > tFields = split( iItem, ':' );
            MORIS_ERROR( tFields.size() <= 3, "feature_query_points: point '%s' has more fields than coordinates, phase and radius.", iItem.c_str() );

            std::vector< std::string > tCoordinates = split( tFields[ 0 ], ',' );
            Feature_Query_Point        tPoint;
            tPoint.mCoordinates.set_size( 1, tCoordinates.size() );
            for ( uint iDim = 0; iDim < tCoordinates.size(); iDim++ )
            {
                tPoint.mCoordinates( iDim ) = to_real( tCoordinates[ iDim ], iItem );
            }
            if ( tFields.size() > 1 and not tFields[ 1 ].empty() )
            {
                real tPhase   = to_real( tFields[ 1 ], iItem );
                tPoint.mPhase = sint( std::lround( tPhase ) );
                MORIS_ERROR( real( tPoint.mPhase ) == tPhase, "feature_query_points: phase of point '%s' must be +1, -1 or 0.", iItem.c_str() );
            }
            if ( tFields.size() > 2 and not tFields[ 2 ].empty() )
            {
                tPoint.mRadius = to_real( tFields[ 2 ], iItem );
            }
            check_query_point( tPoint, iItem );
            tPoints.push_back( tPoint );
        }
        return tPoints;
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< Feature_Query_Point >
    feature_query_points_from_rows( const Matrix< DDRMat >& aRows )
    {
        Vector< Feature_Query_Point > tPoints;
        if ( aRows.n_rows() == 0 )
        {
            return tPoints;
        }
        MORIS_ERROR( aRows.n_cols() == 4 or aRows.n_cols() == 5,
                "feature_query_function: rows need x, y[, z], phase, radius (4 or 5 columns), got %u.",
                uint( aRows.n_cols() ) );
        uint tNumDims = aRows.n_cols() - 2;
        for ( uint iRow = 0; iRow < aRows.n_rows(); iRow++ )
        {
            Feature_Query_Point tPoint;
            tPoint.mCoordinates.set_size( 1, tNumDims );
            for ( uint iDim = 0; iDim < tNumDims; iDim++ )
            {
                tPoint.mCoordinates( iDim ) = aRows( iRow, iDim );
            }
            tPoint.mPhase  = sint( std::lround( aRows( iRow, tNumDims ) ) );
            tPoint.mRadius = aRows( iRow, tNumDims + 1 );
            check_query_point( tPoint, "row " + std::to_string( iRow ) + " of feature_query_function" );
            tPoints.push_back( tPoint );
        }
        return tPoints;
    }

    //--------------------------------------------------------------------------------------------------------------

    sint
    represented_region_at_point(
            const Point_Field_Function& aPhi,
            const Matrix< DDRMat >&     aCorners,
            real                        aInterfaceTolerance,
            const Matrix< DDRMat >&     aPoint,
            bool                        aExactRoots )
    {
        uint tNumDims = aCorners.n_cols();
        MORIS_ERROR( tNumDims == 2 or tNumDims == 3, "represented_region_at_point() - 2D or 3D cells only." );
        MORIS_ERROR( aPoint.numel() == tNumDims, "represented_region_at_point() - point and cell dimensions differ." );

        Represented_Region tRegion( aPhi, tNumDims, aInterfaceTolerance, aExactRoots );

        // the point, clamped into the cell's box (a point on the cell boundary rounds into it)
        Point3 tPoint = { 0.0, 0.0, 0.0 };
        for ( uint iDim = 0; iDim < tNumDims; iDim++ )
        {
            real tLower = 0.0;
            real tUpper = 0.0;
            corner_bounds( aCorners, iDim, tLower, tUpper );
            tPoint[ iDim ] = std::min( std::max( aPoint( iDim ), tLower ), tUpper );
        }

        // template vertices: corners, cell centre (and face centres in 3D)
        uint                           tNumCorners = 1u << tNumDims;
        std::vector< Template_Vertex > tCorners( tNumCorners );
        Point3                         tCentre = { 0.0, 0.0, 0.0 };
        for ( uint iCorner = 0; iCorner < tNumCorners; iCorner++ )
        {
            Point3 tX = { 0.0, 0.0, 0.0 };
            for ( uint iDim = 0; iDim < tNumDims; iDim++ )
            {
                tX[ iDim ] = aCorners( iCorner, iDim );
                tCentre[ iDim ] += aCorners( iCorner, iDim ) / real( tNumCorners );
            }
            tCorners[ iCorner ] = tRegion.vertex( tX );
        }
        Template_Vertex tCentreVertex = tRegion.vertex( tCentre );

        if ( tNumDims == 2 )
        {
            for ( uint iEdge = 0; iEdge < 4; iEdge++ )
            {
                std::vector< Template_Vertex > tTriangle = { tCorners[ iEdge ], tCorners[ ( iEdge + 1 ) % 4 ], tCentreVertex };
                if ( simplex_contains( { tTriangle[ 0 ].mX, tTriangle[ 1 ].mX, tTriangle[ 2 ].mX }, tPoint, 2 ) )
                {
                    return tRegion.region_in_simplex( tTriangle, tPoint );
                }
            }
        }
        else
        {
            // HEX8 faces (exodus), corners in cyclic order
            static const uint tFaces[ 6 ][ 4 ] = { { 0, 1, 5, 4 }, { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 0, 4, 7, 3 }, { 0, 3, 2, 1 }, { 4, 5, 6, 7 } };
            for ( uint iFaceIndex = 0; iFaceIndex < 6; iFaceIndex++ )
            {
                const uint* iFace       = tFaces[ iFaceIndex ];
                Point3      tFaceCentre = { 0.0, 0.0, 0.0 };
                for ( uint iCorner = 0; iCorner < 4; iCorner++ )
                {
                    for ( uint iDim = 0; iDim < 3; iDim++ )
                    {
                        tFaceCentre[ iDim ] += 0.25 * tCorners[ iFace[ iCorner ] ].mX[ iDim ];
                    }
                }
                Template_Vertex tFaceVertex = tRegion.vertex( tFaceCentre );
                for ( uint iEdge = 0; iEdge < 4; iEdge++ )
                {
                    std::vector< Template_Vertex > tTet = { tCentreVertex, tFaceVertex, tCorners[ iFace[ iEdge ] ], tCorners[ iFace[ ( iEdge + 1 ) % 4 ] ] };
                    if ( simplex_contains( { tTet[ 0 ].mX, tTet[ 1 ].mX, tTet[ 2 ].mX, tTet[ 3 ].mX }, tPoint, 3 ) )
                    {
                        return tRegion.region_in_simplex( tTet, tPoint );
                    }
                }
            }
        }

        MORIS_ERROR( false, "represented_region_at_point() - the point lies in no piece of the cell's subdivision." );
        return 0;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    simplex_contains_point( const Matrix< DDRMat >& aVertexCoordinates, const Matrix< DDRMat >& aPoint )
    {
        uint tNumDims = aVertexCoordinates.n_cols();
        MORIS_ERROR( aVertexCoordinates.n_rows() > tNumDims, "simplex_contains_point() - need d + 1 vertices." );
        std::vector< Point3 > tVertices( tNumDims + 1, Point3{ 0.0, 0.0, 0.0 } );
        Point3                tPoint = { 0.0, 0.0, 0.0 };
        for ( uint iDim = 0; iDim < tNumDims; iDim++ )
        {
            for ( uint iVertex = 0; iVertex <= tNumDims; iVertex++ )
            {
                tVertices[ iVertex ][ iDim ] = aVertexCoordinates( iVertex, iDim );
            }
            tPoint[ iDim ] = aPoint( iDim );
        }
        return simplex_contains( tVertices, tPoint, tNumDims );
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    cell_contains_point( const Matrix< DDRMat >& aCorners, const Matrix< DDRMat >& aPoint )
    {
        for ( uint iDim = 0; iDim < aCorners.n_cols(); iDim++ )
        {
            real tLower = 0.0;
            real tUpper = 0.0;
            corner_bounds( aCorners, iDim, tLower, tUpper );
            real tTolerance = 1.0e-12 * ( tUpper - tLower );
            if ( aPoint( iDim ) < tLower - tTolerance or aPoint( iDim ) > tUpper + tTolerance )
            {
                return false;
            }
        }
        return true;
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    cell_distance_to_point( const Matrix< DDRMat >& aCorners, const Matrix< DDRMat >& aPoint )
    {
        real tSquared = 0.0;
        for ( uint iDim = 0; iDim < aCorners.n_cols(); iDim++ )
        {
            real tLower = 0.0;
            real tUpper = 0.0;
            corner_bounds( aCorners, iDim, tLower, tUpper );
            real tGap = std::max( { tLower - aPoint( iDim ), 0.0, aPoint( iDim ) - tUpper } );
            tSquared += tGap * tGap;
        }
        return std::sqrt( tSquared );
    }

}    // namespace moris::gen
