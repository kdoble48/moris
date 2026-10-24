/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Cut_Pathology_Reference.hpp
 *
 * Reference geometry for the conformal-cut pathology sweeps: the trilinear corner-value and plane level-set fields,
 * GEN geometries in the three field regimes, and the exact moments of the piecewise-linear interpolant on the 24-tet
 * subdivision of [0,1]^3 that the decomposition is checked against.
 *
 */

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "cl_Communication_Tools.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR_Parameters.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_HMR_Background_Mesh_Base.hpp"
#include "cl_HMR_Background_Element_Base.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_GEN_Level_Set_Geometry.hpp"
#include "cl_GEN_User_Defined_Field.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_MTK_Set.hpp"
#include "cl_MTK_Cluster.hpp"
#include "cl_MTK_Vertex_Interpolation.hpp"
#include "cl_XTK_Pathology_Test_Tools.hpp"

namespace moris::xtk
{
    namespace cut_pathology_test
    {
        using namespace pathology_test;

        using Corner_Values = std::array< real, 8 >;

        // corner i of [0,1]^3 sits at ( i & 1, ( i >> 1 ) & 1, ( i >> 2 ) & 1 )
        inline Corner_Values gCornersA = {};
        inline Corner_Values gCornersB = {};

        // planes a x + b y + c z + d for the multi-geometry sweeps (one field function per geometry slot)
        using Plane                  = std::array< real, 4 >;
        const uint             gMaxPlanes = 6;
        inline std::array< Plane, 6 > gPlanes    = {};

        template< uint N >
        inline real
        plane_field( const Matrix< DDRMat >& aPoint, const Vector< real >& aParameters )
        {
            return gPlanes[ N ][ 0 ] * aPoint( 0 ) + gPlanes[ N ][ 1 ] * aPoint( 1 ) + gPlanes[ N ][ 2 ] * aPoint( 2 ) + gPlanes[ N ][ 3 ];
        }

        const gen::Field_Function gPlaneFields[ 6 ] = {
            &plane_field< 0 >, &plane_field< 1 >, &plane_field< 2 >, &plane_field< 3 >, &plane_field< 4 >, &plane_field< 5 >
        };

        // volumes below this count as slivers (relative to the unit background cell)
        const real gSliverVolume = 1.0e-12;

        // tolerance on volumes compared with the reference
        const real gVolumeTolerance = 1.0e-10;

        inline Point
        corner_point( uint aCorner )
        {
            return { real( aCorner & 1 ), real( ( aCorner >> 1 ) & 1 ), real( ( aCorner >> 2 ) & 1 ) };
        }

        inline real
        trilinear( const Matrix< DDRMat >& aPoint, const Corner_Values& aCorners )
        {
            real tValue = 0.0;
            for ( uint iCorner = 0; iCorner < 8; iCorner++ )
            {
                real tWeight = 1.0;
                for ( uint iDim = 0; iDim < 3; iDim++ )
                {
                    tWeight *= ( ( iCorner >> iDim ) & 1 ) ? aPoint( iDim ) : 1.0 - aPoint( iDim );
                }
                tValue += aCorners[ iCorner ] * tWeight;
            }
            return tValue;
        }

        inline real
        field_a( const Matrix< DDRMat >& aPoint, const Vector< real >& aParameters )
        {
            return trilinear( aPoint, gCornersA );
        }

        inline real
        field_b( const Matrix< DDRMat >& aPoint, const Vector< real >& aParameters )
        {
            return trilinear( aPoint, gCornersB );
        }

        inline std::shared_ptr< gen::Geometry >
        create_geometry( gen::Field_Function aFunction, Regime aRegime )
        {
            gen::Level_Set_Parameters tParameters;
            tParameters.mDiscretizationIndex          = aRegime == Regime::STORED ? -1 : -2;
            tParameters.mUseMultilinearInterpolation = aRegime == Regime::MULTILINEAR;

            std::shared_ptr< gen::Field > tField = std::make_shared< gen::User_Defined_Field >( aFunction );
            return std::make_shared< gen::Level_Set_Geometry >( tField, tParameters );
        }

        // ----------------------------------------------------------------------------------
        // independent reference: negative volume of the piecewise-linear interpolant on the 24-tet subdivision

        inline real
        tet_volume( const Point& aP0, const Point& aP1, const Point& aP2, const Point& aP3 )
        {
            real tA[ 3 ][ 3 ];
            for ( uint iDim = 0; iDim < 3; iDim++ )
            {
                tA[ 0 ][ iDim ] = aP1[ iDim ] - aP0[ iDim ];
                tA[ 1 ][ iDim ] = aP2[ iDim ] - aP0[ iDim ];
                tA[ 2 ][ iDim ] = aP3[ iDim ] - aP0[ iDim ];
            }
            real tDet = tA[ 0 ][ 0 ] * ( tA[ 1 ][ 1 ] * tA[ 2 ][ 2 ] - tA[ 1 ][ 2 ] * tA[ 2 ][ 1 ] )
                      - tA[ 0 ][ 1 ] * ( tA[ 1 ][ 0 ] * tA[ 2 ][ 2 ] - tA[ 1 ][ 2 ] * tA[ 2 ][ 0 ] )
                      + tA[ 0 ][ 2 ] * ( tA[ 1 ][ 0 ] * tA[ 2 ][ 1 ] - tA[ 1 ][ 1 ] * tA[ 2 ][ 0 ] );
            return std::abs( tDet ) / 6.0;
        }

        // point where the linear function crosses zero on the edge from a negative vertex u to a non-negative vertex w
        inline Point
        zero_on_edge( const Point& aU, real aFU, const Point& aW, real aFW )
        {
            if ( aFW == 0.0 )
            {
                return aW;
            }
            real  tT = aFU / ( aFU - aFW );
            Point tP;
            for ( uint iDim = 0; iDim < 3; iDim++ )
            {
                tP[ iDim ] = aU[ iDim ] + tT * ( aW[ iDim ] - aU[ iDim ] );
            }
            return tP;
        }

        // integrals of ( 1, x, y, z, x^2 ) over a domain: the quadrature moments checked by the sweeps
        using Moments = std::array< real, 5 >;

        inline Moments
        operator+( const Moments& aA, const Moments& aB )
        {
            Moments tSum;
            for ( uint iM = 0; iM < 5; iM++ )
            {
                tSum[ iM ] = aA[ iM ] + aB[ iM ];
            }
            return tSum;
        }

        inline Moments
        operator-( const Moments& aA, const Moments& aB )
        {
            Moments tDifference;
            for ( uint iM = 0; iM < 5; iM++ )
            {
                tDifference[ iM ] = aA[ iM ] - aB[ iM ];
            }
            return tDifference;
        }

        inline real
        max_abs_difference( const Moments& aA, const Moments& aB )
        {
            real tMax = 0.0;
            for ( uint iM = 0; iM < 5; iM++ )
            {
                tMax = std::max( tMax, std::abs( aA[ iM ] - aB[ iM ] ) );
            }
            return tMax;
        }

        // exact moments of a tet for a given (signed or absolute) volume: int x = V x_c, int x^2 = V/20 ( sum x_i^2 + (sum x_i)^2 )
        inline Moments
        tet_moments( const Point& aP0, const Point& aP1, const Point& aP2, const Point& aP3, real aVolume )
        {
            Moments tMoments;
            tMoments[ 0 ] = aVolume;
            for ( uint iDim = 0; iDim < 3; iDim++ )
            {
                tMoments[ 1 + iDim ] = aVolume * ( aP0[ iDim ] + aP1[ iDim ] + aP2[ iDim ] + aP3[ iDim ] ) / 4.0;
            }
            real tSum        = aP0[ 0 ] + aP1[ 0 ] + aP2[ 0 ] + aP3[ 0 ];
            real tSumSquares = aP0[ 0 ] * aP0[ 0 ] + aP1[ 0 ] * aP1[ 0 ] + aP2[ 0 ] * aP2[ 0 ] + aP3[ 0 ] * aP3[ 0 ];
            tMoments[ 4 ]    = aVolume / 20.0 * ( tSumSquares + tSum * tSum );
            return tMoments;
        }

        inline Moments
        tet_moments( const Point& aP0, const Point& aP1, const Point& aP2, const Point& aP3 )
        {
            return tet_moments( aP0, aP1, aP2, aP3, tet_volume( aP0, aP1, aP2, aP3 ) );
        }

        // moments of { f < 0 } in a tet with a linear f given by its vertex values
        inline Moments
        tet_negative_moments( const std::array< Point, 4 >& aP, const std::array< real, 4 >& aF )
        {
            Vector< uint > tNeg;
            Vector< uint > tNonNeg;
            bool           tHasPositive = false;
            for ( uint iV = 0; iV < 4; iV++ )
            {
                if ( aF[ iV ] < 0.0 )
                {
                    tNeg.push_back( iV );
                }
                else
                {
                    tNonNeg.push_back( iV );
                    tHasPositive = tHasPositive or aF[ iV ] > 0.0;
                }
            }

            Moments tWhole = tet_moments( aP[ 0 ], aP[ 1 ], aP[ 2 ], aP[ 3 ] );
            if ( tNeg.size() == 0 )
            {
                return Moments{};
            }
            if ( not tHasPositive )
            {
                return tWhole;
            }
            if ( tNeg.size() == 1 )
            {
                uint tU = tNeg( 0 );
                return tet_moments( aP[ tU ],
                        zero_on_edge( aP[ tU ], aF[ tU ], aP[ tNonNeg( 0 ) ], aF[ tNonNeg( 0 ) ] ),
                        zero_on_edge( aP[ tU ], aF[ tU ], aP[ tNonNeg( 1 ) ], aF[ tNonNeg( 1 ) ] ),
                        zero_on_edge( aP[ tU ], aF[ tU ], aP[ tNonNeg( 2 ) ], aF[ tNonNeg( 2 ) ] ) );
            }
            if ( tNeg.size() == 3 )
            {
                // the single non-negative vertex is positive (otherwise handled above): cut off its corner tet
                uint tW = tNonNeg( 0 );
                return tWhole
                     - tet_moments( aP[ tW ],
                             zero_on_edge( aP[ tNeg( 0 ) ], aF[ tNeg( 0 ) ], aP[ tW ], aF[ tW ] ),
                             zero_on_edge( aP[ tNeg( 1 ) ], aF[ tNeg( 1 ) ], aP[ tW ], aF[ tW ] ),
                             zero_on_edge( aP[ tNeg( 2 ) ], aF[ tNeg( 2 ) ], aP[ tW ], aF[ tW ] ) );
            }

            // two negative vertices A, B and two non-negative C, D: the negative part is the prism
            // ( A, AC, AD ) - ( B, BC, BD ), split into three tets
            uint  tA   = tNeg( 0 );
            uint  tB   = tNeg( 1 );
            uint  tC   = tNonNeg( 0 );
            uint  tD   = tNonNeg( 1 );
            Point tAC  = zero_on_edge( aP[ tA ], aF[ tA ], aP[ tC ], aF[ tC ] );
            Point tAD  = zero_on_edge( aP[ tA ], aF[ tA ], aP[ tD ], aF[ tD ] );
            Point tBC  = zero_on_edge( aP[ tB ], aF[ tB ], aP[ tC ], aF[ tC ] );
            Point tBD  = zero_on_edge( aP[ tB ], aF[ tB ], aP[ tD ], aF[ tD ] );
            return tet_moments( aP[ tA ], tAC, tAD, aP[ tB ] ) + tet_moments( tAC, tAD, aP[ tB ], tBC ) + tet_moments( tAD, aP[ tB ], tBC, tBD );
        }

        // moments of the negative part of the piecewise-linear interpolant on the 24-tet subdivision of [0,1]^3
        inline Moments
        reference_negative_moments( const Corner_Values& aCorners )
        {
            Point tCellCentre = { 0.5, 0.5, 0.5 };
            real  tCellValue  = 0.0;
            for ( real iValue : aCorners )
            {
                tCellValue += iValue / 8.0;
            }

            Moments tMoments = {};
            for ( uint iAxis = 0; iAxis < 3; iAxis++ )
            {
                uint tB = ( iAxis + 1 ) % 3;
                uint tC = ( iAxis + 2 ) % 3;
                for ( uint iSide = 0; iSide < 2; iSide++ )
                {
                    // face corners in cyclic order
                    std::array< uint, 4 > tFaceCorners;
                    const uint            tCycle[ 4 ][ 2 ] = { { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } };
                    for ( uint iK = 0; iK < 4; iK++ )
                    {
                        tFaceCorners[ iK ] = ( iSide << iAxis ) | ( tCycle[ iK ][ 0 ] << tB ) | ( tCycle[ iK ][ 1 ] << tC );
                    }

                    Point tFaceCentre;
                    tFaceCentre[ iAxis ] = iSide;
                    tFaceCentre[ tB ]    = 0.5;
                    tFaceCentre[ tC ]    = 0.5;
                    real tFaceValue      = 0.0;
                    for ( uint iCorner : tFaceCorners )
                    {
                        tFaceValue += aCorners[ iCorner ] / 4.0;
                    }

                    for ( uint iK = 0; iK < 4; iK++ )
                    {
                        uint tFirst  = tFaceCorners[ iK ];
                        uint tSecond = tFaceCorners[ ( iK + 1 ) % 4 ];
                        tMoments = tMoments
                                 + tet_negative_moments(
                                         { corner_point( tFirst ), corner_point( tSecond ), tFaceCentre, tCellCentre },
                                         { aCorners[ tFirst ], aCorners[ tSecond ], tFaceValue, tCellValue } );
                    }
                }
            }
            return tMoments;
        }

        // negative volume of the piecewise-linear interpolant on the 24-tet subdivision of [0,1]^3
        inline real
        reference_negative_volume( const Corner_Values& aCorners )
        {
            return reference_negative_moments( aCorners )[ 0 ];
        }
    }    // namespace cut_pathology_test
}    // namespace moris::xtk
