/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_XTK_Cut_Pathologies.cpp
 *
 * Characterization sweep of the conformal decomposition (regular subdivision + node hierarchy) on a single HEX8
 * background cell [0,1]^3. The level set is the trilinear interpolant of 8 corner values; every combination of corner
 * values from a small value set is decomposed and checked against invariants that must hold for any input:
 *
 *   - the decomposition does not throw,
 *   - the integration cells tile the cell (total volume 1),
 *   - every integration cell has a bulk phase,
 *   - no integration cell has both a POSITIVE and a NEGATIVE vertex for any geometry,
 *   - every facet has one or two cells, and facets with one cell lie on the background cell boundary,
 *   - (interpolated and analytic fields) the negative volume equals the volume of the piecewise-linear interpolant on
 *     the 24-tet subdivision (face and cell centre values are corner averages), computed independently below.
 *
 * Three GEN field regimes are swept: STORED (discretization_mesh_index = -1: derived values interpolated from their
 * parents, the regime of discretized level sets), ANALYTIC (derived values evaluated exactly) and MULTILINEAR
 * (use_multilinear_interpolation). With two geometries the STORED regime must not depend on the geometry order: cutting
 * by one geometry only subdivides the piecewise-linear representation of the other.
 *
 * Three plane sweeps report how cuts accumulate (1 to 6 planes, both orders), two parallel planes approaching, and two
 * planes crossing at a shrinking angle. They print the quantities that matter for an immersed discretization, where
 * the integration cells are quadrature domains only: cell count, smallest cell and subphase volume, collapsed cells,
 * and the number and smallest area of interface facets.
 *
 * A multi-rank sweep ([XTK_PATHOLOGY_MPI]) decomposes a 4 x 2 x 2 box [0,2] x [0,1] x [0,1] with planes on, near
 * and across the processor boundary (x = 1 on 2 ranks). Besides the global invariants (total volume, negative volume
 * against the exact value, the same result on every number of ranks) it gathers every integration vertex on rank 0
 * and checks that the ranks agree on it: one id per location, and per id the same coordinates, owner and geometric
 * region. Run it on 1, 2 and 4 ranks, e.g. mpirun -n 2 XTK-test.exe "[XTK_PATHOLOGY_MPI]", and compare the CSV files
 * (MORIS_XTK_PATHOLOGY_MPI_CSV=<file>).
 *
 * The sweeps are hidden ([.]) because they characterize current behaviour rather than gate it. Run one with
 *   XTK-test.exe "[XTK_PATHOLOGY_HEX_SIGNS]"
 * Set MORIS_XTK_PATHOLOGY_CSV=<file> to write one row per decomposed case.
 *
 */

#include "catch.hpp"

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

#include "cl_XTK_Pathology_Test_Tools.hpp"
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

namespace moris::xtk
{
    using namespace pathology_test;

    namespace
    {
        using Corner_Values = std::array< real, 8 >;

        // corner i of [0,1]^3 sits at ( i & 1, ( i >> 1 ) & 1, ( i >> 2 ) & 1 )
        Corner_Values gCornersA = {};
        Corner_Values gCornersB = {};

        // planes a x + b y + c z + d for the multi-geometry sweeps (one field function per geometry slot)
        using Plane                  = std::array< real, 4 >;
        const uint             gMaxPlanes = 6;
        std::array< Plane, 6 > gPlanes    = {};

        template< uint N >
        real
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

        Point
        corner_point( uint aCorner )
        {
            return { real( aCorner & 1 ), real( ( aCorner >> 1 ) & 1 ), real( ( aCorner >> 2 ) & 1 ) };
        }

        real
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

        real
        field_a( const Matrix< DDRMat >& aPoint, const Vector< real >& aParameters )
        {
            return trilinear( aPoint, gCornersA );
        }

        real
        field_b( const Matrix< DDRMat >& aPoint, const Vector< real >& aParameters )
        {
            return trilinear( aPoint, gCornersB );
        }

        std::shared_ptr< gen::Geometry >
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

        real
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
        Point
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

        Moments
        operator+( const Moments& aA, const Moments& aB )
        {
            Moments tSum;
            for ( uint iM = 0; iM < 5; iM++ )
            {
                tSum[ iM ] = aA[ iM ] + aB[ iM ];
            }
            return tSum;
        }

        Moments
        operator-( const Moments& aA, const Moments& aB )
        {
            Moments tDifference;
            for ( uint iM = 0; iM < 5; iM++ )
            {
                tDifference[ iM ] = aA[ iM ] - aB[ iM ];
            }
            return tDifference;
        }

        real
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
        Moments
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

        Moments
        tet_moments( const Point& aP0, const Point& aP1, const Point& aP2, const Point& aP3 )
        {
            return tet_moments( aP0, aP1, aP2, aP3, tet_volume( aP0, aP1, aP2, aP3 ) );
        }

        // moments of { f < 0 } in a tet with a linear f given by its vertex values
        Moments
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
        Moments
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
        real
        reference_negative_volume( const Corner_Values& aCorners )
        {
            return reference_negative_moments( aCorners )[ 0 ];
        }

        // ----------------------------------------------------------------------------------
        // decomposition and invariant checks

        struct Case_Result
        {
            bool        mThrew   = false;
            std::string mMessage = "";

            uint mNumCells        = 0;
            uint mNumNoPhase      = 0;    // cells without a valid bulk phase
            uint mNumMixedSign    = 0;    // cells with a POSITIVE and a NEGATIVE vertex for some geometry
            uint mNumBadFacets    = 0;    // facets with 0 or > 2 cells, or one cell off the background boundary
            uint mNumSlivers      = 0;    // cells with volume < gSliverVolume
            uint mNumUndetermined = 0;    // cells with every vertex on the interface of some geometry

            real mTotalVolume   = 0.0;
            real mMinCellVolume = MORIS_REAL_MAX;

            // quantities that matter for an immersed discretization (the IG cells are integration domains only):
            // small subphases condition the enriched basis, small interface facets carry the interface normal
            uint mNumNonPositive         = 0;    // cells with volume <= 0 (collapsed or inverted)
            uint mNumSubphases           = 0;
            real mMinSubphaseVolume      = MORIS_REAL_MAX;
            uint mNumInterfaceFacets     = 0;    // facets between cells of different bulk phase
            real mMinInterfaceFacetArea  = MORIS_REAL_MAX;

            // negative volume per geometry, from the vertex labels of each cell (cells with undetermined sign excluded)
            Vector< real > mNegativeVolume;

            // moments ( 1, x, y, z, x^2 ) of the negative cells per geometry, with each cell's signed volume (the
            // measure MORIS returns) and with its absolute volume (an integrator using |det J|)
            Vector< Moments > mNegativeMomentsSigned;
            Vector< Moments > mNegativeMomentsAbsolute;

            // volume per sign combination ( bit g set = positive for geometry g ), from the vertex labels
            Vector< real > mSignVolume;
        };

        bool
        on_background_boundary( const Vector< mtk::Vertex* >& aVertices )
        {
            for ( uint iDim = 0; iDim < 3; iDim++ )
            {
                for ( real tSide : { 0.0, 1.0 } )
                {
                    bool tAllOnSide = true;
                    for ( mtk::Vertex* iVertex : aVertices )
                    {
                        tAllOnSide = tAllOnSide and std::abs( iVertex->get_coords()( iDim ) - tSide ) < 1.0e-12;
                    }
                    if ( tAllOnSide )
                    {
                        return true;
                    }
                }
            }
            return false;
        }

        // moments of one integration cell with a given volume: a TET4, or the uncut background HEX8 [0,1]^3
        Moments
        cell_moments( const Vector< mtk::Vertex* >& aVertices, real aVolume )
        {
            if ( aVertices.size() == 8 )
            {
                return { aVolume, 0.5 * aVolume, 0.5 * aVolume, 0.5 * aVolume, aVolume / 3.0 };
            }
            MORIS_ERROR( aVertices.size() == 4, "cell_moments() - expected a TET4 or the unit HEX8" );
            std::array< Point, 4 > tPoints;
            for ( uint iV = 0; iV < 4; iV++ )
            {
                Matrix< DDRMat > tCoords = aVertices( iV )->get_coords();
                tPoints[ iV ]            = { tCoords( 0 ), tCoords( 1 ), tCoords( 2 ) };
            }
            return tet_moments( tPoints[ 0 ], tPoints[ 1 ], tPoints[ 2 ], tPoints[ 3 ], aVolume );
        }

        // area of a planar facet given by its vertices in order (triangle fan)
        real
        facet_area( const Vector< mtk::Vertex* >& aVertices )
        {
            Matrix< DDRMat > tP0  = aVertices( 0 )->get_coords();
            real             tSum[ 3 ] = { 0.0, 0.0, 0.0 };
            for ( uint iV = 1; iV + 1 < aVertices.size(); iV++ )
            {
                Matrix< DDRMat > tA = aVertices( iV )->get_coords() - tP0;
                Matrix< DDRMat > tB = aVertices( iV + 1 )->get_coords() - tP0;
                tSum[ 0 ] += tA( 1 ) * tB( 2 ) - tA( 2 ) * tB( 1 );
                tSum[ 1 ] += tA( 2 ) * tB( 0 ) - tA( 0 ) * tB( 2 );
                tSum[ 2 ] += tA( 0 ) * tB( 1 ) - tA( 1 ) * tB( 0 );
            }
            return 0.5 * std::sqrt( tSum[ 0 ] * tSum[ 0 ] + tSum[ 1 ] * tSum[ 1 ] + tSum[ 2 ] * tSum[ 2 ] );
        }

        Case_Result
        decompose_case( hmr::Interpolation_Mesh_HMR* aMesh, const Vector< std::shared_ptr< gen::Geometry > >& aGeometries )
        {
            Case_Result tResult;
            uint        tNumGeometries = aGeometries.size();
            tResult.mNegativeVolume.resize( tNumGeometries, 0.0 );
            tResult.mNegativeMomentsSigned.resize( tNumGeometries, Moments{} );
            tResult.mNegativeMomentsAbsolute.resize( tNumGeometries, Moments{} );
            tResult.mSignVolume.resize( 1u << tNumGeometries, 0.0 );

            try
            {
                gen::Geometry_Engine_Parameters tParameters;
                tParameters.mGeometries = aGeometries;
                gen::Geometry_Engine tGeometryEngine( aMesh, tParameters );

                Model tModel( 3, aMesh, &tGeometryEngine );
                tModel.mVerbose = false;
                tModel.decompose( { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 } );

                Cut_Integration_Mesh* tCutMesh     = tModel.get_cut_integration_mesh();
                moris_index           tNumPhases   = (moris_index)tGeometryEngine.get_num_bulk_phase();
                uint                  tNumSubphase = tCutMesh->get_num_subphases();

                tResult.mNumSubphases = tNumSubphase;
                for ( uint iSP = 0; iSP < tNumSubphase; iSP++ )
                {
                    real tSubphaseVolume = 0.0;
                    for ( mtk::Cell* iCell : tCutMesh->get_subphase_ig_cells( iSP )->mIgCellGroup )
                    {
                        real tVolume = iCell->compute_cell_measure();
                        tSubphaseVolume += tVolume;
                        if ( tVolume <= 0.0 )
                        {
                            tResult.mNumNonPositive++;
                        }
                        tResult.mNumCells++;
                        tResult.mTotalVolume += tVolume;
                        tResult.mMinCellVolume = std::min( tResult.mMinCellVolume, tVolume );
                        if ( tVolume < gSliverVolume )
                        {
                            tResult.mNumSlivers++;
                        }

                        moris_index tBulkPhase = tCutMesh->get_cell_bulk_phase( iCell->get_index() );
                        if ( tBulkPhase < 0 or tBulkPhase >= tNumPhases )
                        {
                            tResult.mNumNoPhase++;
                        }

                        // vertex labels per geometry
                        Vector< mtk::Vertex* > tVertices    = iCell->get_vertex_pointers();
                        bool                   tMixed       = false;
                        bool                   tUndetermined = false;
                        uint                   tSignBits    = 0;
                        for ( uint iGeom = 0; iGeom < tNumGeometries; iGeom++ )
                        {
                            bool tHasPositive = false;
                            bool tHasNegative = false;
                            for ( mtk::Vertex* iVertex : tVertices )
                            {
                                gen::Geometric_Region tRegion = tGeometryEngine.get_geometric_region( iGeom, iVertex->get_index(), iVertex->get_coords() );
                                tHasPositive = tHasPositive or tRegion == gen::Geometric_Region::POSITIVE;
                                tHasNegative = tHasNegative or tRegion == gen::Geometric_Region::NEGATIVE;
                            }
                            tMixed        = tMixed or ( tHasPositive and tHasNegative );
                            tUndetermined = tUndetermined or ( not tHasPositive and not tHasNegative );
                            if ( tHasNegative and not tHasPositive )
                            {
                                tResult.mNegativeVolume( iGeom ) += tVolume;
                                tResult.mNegativeMomentsSigned( iGeom )   = tResult.mNegativeMomentsSigned( iGeom ) + cell_moments( tVertices, tVolume );
                                tResult.mNegativeMomentsAbsolute( iGeom ) = tResult.mNegativeMomentsAbsolute( iGeom ) + cell_moments( tVertices, std::abs( tVolume ) );
                            }
                            if ( tHasPositive )
                            {
                                tSignBits |= 1u << iGeom;
                            }
                        }
                        tResult.mNumMixedSign += tMixed;
                        tResult.mNumUndetermined += tUndetermined;
                        if ( not tMixed and not tUndetermined )
                        {
                            tResult.mSignVolume( tSignBits ) += tVolume;
                        }
                    }
                    tResult.mMinSubphaseVolume = std::min( tResult.mMinSubphaseVolume, tSubphaseVolume );
                }

                std::shared_ptr< Facet_Based_Connectivity > tFacets = tCutMesh->get_face_connectivity();
                if ( tFacets != nullptr )
                {
                    for ( uint iFacet = 0; iFacet < tFacets->mFacetToCell.size(); iFacet++ )
                    {
                        uint tNumAttached = tFacets->mFacetToCell( iFacet ).size();
                        bool tBad         = tNumAttached == 0 or tNumAttached > 2
                                 or ( tNumAttached == 1 and not on_background_boundary( tFacets->mFacetVertices( iFacet ) ) );
                        tResult.mNumBadFacets += tBad;

                        // interface facet: two cells of different bulk phase (cf. Integration_Mesh_Generator::deduce_interfaces)
                        if ( tNumAttached == 2
                                and tCutMesh->get_cell_bulk_phase( tFacets->mFacetToCell( iFacet )( 0 )->get_index() )
                                            != tCutMesh->get_cell_bulk_phase( tFacets->mFacetToCell( iFacet )( 1 )->get_index() ) )
                        {
                            tResult.mNumInterfaceFacets++;
                            tResult.mMinInterfaceFacetArea = std::min( tResult.mMinInterfaceFacetArea, facet_area( tFacets->mFacetVertices( iFacet ) ) );
                        }
                    }
                }
            }
            catch ( const std::exception& aException )
            {
                tResult.mThrew   = true;
                tResult.mMessage = aException.what();
            }
            return tResult;
        }

        // ----------------------------------------------------------------------------------
        // sweep bookkeeping

        struct Sweep_Summary
        {
            uint mNumCases       = 0;
            uint mNumThrew       = 0;
            uint mNumBadVolume   = 0;
            uint mNumNoPhase     = 0;
            uint mNumMixedSign   = 0;
            uint mNumBadFacets   = 0;
            uint mNumReference   = 0;    // negative volume differs from the reference
            uint mNumOrder       = 0;    // two geometries: result depends on the geometry order
            uint mNumWithSlivers = 0;
            uint mNumUndetermined = 0;
            uint mNumNonPositive  = 0;    // cases with a cell of volume <= 0
            uint mNumMoment       = 0;    // negative-phase moments (signed cell volumes) differ from the reference
            uint mNumMomentAbs    = 0;    // the same with absolute cell volumes
            real mMaxMomentError    = 0.0;
            real mMaxMomentErrorAbs = 0.0;
            uint mNumGuard        = 0;    // cases stopped by the MB-001 guard (every vertex on the interface): expected
            real mMaxReferenceError = 0.0;
            real mMaxOrderDifference = 0.0;
            real mMinCellVolume  = MORIS_REAL_MAX;

            std::map< std::string, uint >          mMessages;
            std::map< std::string, Corner_Values > mFirstCase;
            std::map< std::string, Corner_Values > mFirstViolation;

            void
            record_violation( const std::string& aName, const Corner_Values& aCorners )
            {
                if ( mFirstViolation.find( aName ) == mFirstViolation.end() )
                {
                    mFirstViolation[ aName ] = aCorners;
                }
            }
        };

        std::string
        corners_to_string( const Corner_Values& aCorners )
        {
            std::ostringstream tStream;
            for ( uint iCorner = 0; iCorner < 8; iCorner++ )
            {
                tStream << ( iCorner ? ";" : "" ) << aCorners[ iCorner ];
            }
            return tStream.str();
        }

        // CSV output, open only if MORIS_XTK_PATHOLOGY_CSV is set (a static stream, so it is flushed at exit)
        std::ofstream*
        csv_stream()
        {
            static std::ofstream sStream;
            static bool          sOpened = false;
            if ( not sOpened )
            {
                sOpened           = true;
                const char* tPath = std::getenv( "MORIS_XTK_PATHOLOGY_CSV" );
                if ( tPath != nullptr )
                {
                    sStream.open( tPath );
                    sStream << "sweep,regime,second,order,corners_a,threw,message,num_cells,total_volume,neg_volume_a,ref_a,neg_volume_b,ref_b,"
                               "no_phase,mixed_sign,bad_facets,undetermined,slivers,min_cell_volume\n";
                }
            }
            return sStream.is_open() ? &sStream : nullptr;
        }

        void
        write_row(
                const std::string&   aSweep,
                Regime               aRegime,
                const std::string&   aSecond,
                const std::string&   aOrder,
                const Corner_Values& aCorners,
                const Case_Result&   aResult,
                real                 aReferenceA,
                real                 aReferenceB )
        {
            std::ofstream* tStream = csv_stream();
            if ( tStream == nullptr )
            {
                return;
            }
            std::string tMessage = aResult.mMessage.substr( 0, 160 );
            for ( char& iChar : tMessage )
            {
                if ( iChar == ',' or iChar == '\n' or iChar == '"' )
                {
                    iChar = ' ';
                }
            }
            *tStream << aSweep << "," << regime_name( aRegime ) << "," << aSecond << "," << aOrder << "," << corners_to_string( aCorners ) << ","
                     << aResult.mThrew << "," << tMessage << "," << aResult.mNumCells << "," << aResult.mTotalVolume << ","
                     << ( aResult.mNegativeVolume.size() > 0 ? aResult.mNegativeVolume( 0 ) : 0.0 ) << "," << aReferenceA << ","
                     << ( aResult.mNegativeVolume.size() > 1 ? aResult.mNegativeVolume( 1 ) : 0.0 ) << "," << aReferenceB << ","
                     << aResult.mNumNoPhase << "," << aResult.mNumMixedSign << "," << aResult.mNumBadFacets << "," << aResult.mNumUndetermined << ","
                     << aResult.mNumSlivers << "," << aResult.mMinCellVolume << "\n"
                     << std::flush;    // flushed per row, so the last row locates a crash
        }

        // checks shared by all sweeps; returns false if the case threw
        bool
        accumulate_invariants( Sweep_Summary& aSummary, const Case_Result& aResult, const Corner_Values& aCorners )
        {
            aSummary.mNumCases++;
            if ( aResult.mThrew )
            {
                // the deliberate MB-001 stop (deduce_ig_cell_bulk_phase_from_facets) is expected, not a failure
                if ( aResult.mMessage.find( "lies on a geometry interface, so its bulk phase cannot be determined" ) != std::string::npos )
                {
                    aSummary.mNumGuard++;
                }
                else
                {
                    aSummary.mNumThrew++;
                }
                std::string tKey = aResult.mMessage.substr( 0, 200 );
                if ( aSummary.mMessages[ tKey ]++ == 0 )
                {
                    aSummary.mFirstCase[ tKey ] = aCorners;
                }
                return false;
            }
            if ( std::abs( aResult.mTotalVolume - 1.0 ) > gVolumeTolerance )
            {
                aSummary.mNumBadVolume++;
                aSummary.record_violation( "total volume", aCorners );
            }
            if ( aResult.mNumNoPhase > 0 )
            {
                aSummary.mNumNoPhase++;
                aSummary.record_violation( "no bulk phase", aCorners );
            }
            if ( aResult.mNumMixedSign > 0 )
            {
                aSummary.mNumMixedSign++;
                aSummary.record_violation( "mixed-sign cell", aCorners );
            }
            if ( aResult.mNumBadFacets > 0 )
            {
                aSummary.mNumBadFacets++;
                aSummary.record_violation( "bad facet", aCorners );
            }
            if ( aResult.mNumNonPositive > 0 )
            {
                aSummary.mNumNonPositive++;
                aSummary.record_violation( "cell volume <= 0", aCorners );
            }
            if ( aResult.mNumSlivers > 0 )
            {
                aSummary.mNumWithSlivers++;
            }
            if ( aResult.mNumUndetermined > 0 )
            {
                aSummary.mNumUndetermined++;
            }
            aSummary.mMinCellVolume = std::min( aSummary.mMinCellVolume, aResult.mMinCellVolume );
            return true;
        }

        void
        print_summary( const std::string& aTitle, const Sweep_Summary& aSummary )
        {
            std::cout << "\n=== " << aTitle << " ===\n"
                      << "  cases                         " << aSummary.mNumCases << "\n"
                      << "  threw                         " << aSummary.mNumThrew << "\n"
                      << "  stopped by MB-001 guard       " << aSummary.mNumGuard << "  (expected)\n"
                      << "  cell volume <= 0              " << aSummary.mNumNonPositive << "\n"
                      << "  total volume != 1             " << aSummary.mNumBadVolume << "\n"
                      << "  cell without bulk phase       " << aSummary.mNumNoPhase << "\n"
                      << "  mixed-sign cell               " << aSummary.mNumMixedSign << "\n"
                      << "  bad facet                     " << aSummary.mNumBadFacets << "\n"
                      << "  reference volume mismatch     " << aSummary.mNumReference << "  (max error " << aSummary.mMaxReferenceError << ")\n"
                      << "  moments, signed volumes       " << aSummary.mNumMoment << "  (max error " << aSummary.mMaxMomentError << ")\n"
                      << "  moments, absolute volumes     " << aSummary.mNumMomentAbs << "  (max error " << aSummary.mMaxMomentErrorAbs << ")\n"
                      << "  geometry-order dependence     " << aSummary.mNumOrder << "  (max difference " << aSummary.mMaxOrderDifference << ")\n"
                      << "  cases with slivers (<1e-12)   " << aSummary.mNumWithSlivers << "\n"
                      << "  cases with all-interface cell " << aSummary.mNumUndetermined << "\n"
                      << "  smallest cell volume          " << aSummary.mMinCellVolume << "\n";
            for ( const std::pair< const std::string, uint >& iMessage : aSummary.mMessages )
            {
                std::cout << "  throw x" << iMessage.second << " first at corners [" << corners_to_string( aSummary.mFirstCase.at( iMessage.first ) )
                          << "]: " << iMessage.first << "\n";
            }
            for ( const std::pair< const std::string, Corner_Values >& iViolation : aSummary.mFirstViolation )
            {
                std::cout << "  first " << iViolation.first << " at corners [" << corners_to_string( iViolation.second ) << "]\n";
            }
        }

        // the unit background cell; HMR keeps a reference to its parameters, so they live here as long as the mesh
        struct Unit_Hex
        {
            hmr::Parameters              mParameters;
            hmr::HMR*                    mHMR  = nullptr;
            hmr::Interpolation_Mesh_HMR* mMesh = nullptr;

            Unit_Hex()
            {
                mParameters.set_number_of_elements_per_dimension( { 1, 1, 1 } );
                mParameters.set_domain_dimensions( { 1.0, 1.0, 1.0 } );
                mParameters.set_domain_offset( { 0.0, 0.0, 0.0 } );
                mParameters.set_bspline_truncation( true );
                mParameters.set_lagrange_orders( { 1 } );
                mParameters.set_lagrange_patterns( { 0 } );
                mParameters.set_bspline_orders( { 1 } );
                mParameters.set_bspline_patterns( { 0 } );
                mParameters.set_create_side_sets( true );

                mHMR = new hmr::HMR( mParameters );
                mHMR->finalize();
                mMesh = mHMR->create_interpolation_mesh( 0 );
            }

            ~Unit_Hex()
            {
                delete mMesh;
                delete mHMR;
            }
        };

        // all aValues.size()^8 corner combinations
        Vector< Corner_Values >
        all_corner_combinations( const Vector< real >& aValues )
        {
            Vector< Corner_Values > tCombinations;
            uint                    tNumValues = aValues.size();
            uint                    tNumCases  = 1;
            for ( uint iCorner = 0; iCorner < 8; iCorner++ )
            {
                tNumCases *= tNumValues;
            }
            tCombinations.reserve( tNumCases );
            for ( uint iCase = 0; iCase < tNumCases; iCase++ )
            {
                Corner_Values tCorners;
                uint          tCode = iCase;
                for ( uint iCorner = 0; iCorner < 8; iCorner++ )
                {
                    tCorners[ iCorner ] = aValues( tCode % tNumValues );
                    tCode /= tNumValues;
                }
                tCombinations.push_back( tCorners );
            }
            return tCombinations;
        }

        // single-geometry sweep over all corner combinations of aValues in one regime
        Sweep_Summary
        single_geometry_sweep( const std::string& aSweep, const Vector< real >& aValues, Regime aRegime )
        {
            Unit_Hex                     tHex;
            hmr::Interpolation_Mesh_HMR* tMesh = tHex.mMesh;

            Sweep_Summary tSummary;
            for ( const Corner_Values& iCorners : all_corner_combinations( aValues ) )
            {
                gCornersA = iCorners;

                Case_Result tResult    = decompose_case( tMesh, { create_geometry( &field_a, aRegime ) } );
                real        tReference = reference_negative_volume( iCorners );
                write_row( aSweep, aRegime, "-", "-", iCorners, tResult, tReference, 0.0 );

                if ( accumulate_invariants( tSummary, tResult, iCorners ) and aRegime != Regime::MULTILINEAR )
                {
                    real tError = std::abs( tResult.mNegativeVolume( 0 ) - tReference );
                    tSummary.mMaxReferenceError = std::max( tSummary.mMaxReferenceError, tError );
                    if ( tError > gVolumeTolerance )
                    {
                        tSummary.mNumReference++;
                        tSummary.record_violation( "reference volume", iCorners );
                    }

                    // quadrature moments of the negative phase, with signed and with absolute cell volumes
                    Moments tReferenceMoments = reference_negative_moments( iCorners );
                    real    tMomentError      = max_abs_difference( tResult.mNegativeMomentsSigned( 0 ), tReferenceMoments );
                    real    tMomentErrorAbs   = max_abs_difference( tResult.mNegativeMomentsAbsolute( 0 ), tReferenceMoments );
                    tSummary.mMaxMomentError    = std::max( tSummary.mMaxMomentError, tMomentError );
                    tSummary.mMaxMomentErrorAbs = std::max( tSummary.mMaxMomentErrorAbs, tMomentErrorAbs );
                    if ( tMomentError > gVolumeTolerance )
                    {
                        tSummary.mNumMoment++;
                        tSummary.record_violation( "moments (signed volumes)", iCorners );
                    }
                    if ( tMomentErrorAbs > gVolumeTolerance )
                    {
                        tSummary.mNumMomentAbs++;
                        tSummary.record_violation( "moments (absolute volumes)", iCorners );
                    }
                }
            }
            return tSummary;
        }

        // fixed second geometries for the two-geometry sweep, as corner values of their trilinear interpolant
        Vector< std::pair< std::string, Corner_Values > >
        second_geometries()
        {
            Vector< std::pair< std::string, Corner_Values > > tGeometries;

            // planes, evaluated at the corners (trilinear interpolation reproduces a linear function)
            std::map< std::string, std::array< real, 4 > > tPlanes = {
                { "x=0.5", { 1.0, 0.0, 0.0, -0.5 } },           // through face centres, off the corners
                { "x+y=1", { 1.0, 1.0, 0.0, -1.0 } },           // through two vertical edges (4 corners on the interface)
                { "x+y+z=1", { 1.0, 1.0, 1.0, -1.0 } },         // through three corners
                { "x+y+z=1.5", { 1.0, 1.0, 1.0, -1.5 } },       // through the cell centre, off every corner
                { "z=0", { 0.0, 0.0, 1.0, 0.0 } },              // on the bottom face of the cell
                { "x=0.5+1e-11", { 1.0, 0.0, 0.0, -0.5 - 1.0e-11 } }    // 1e-11 off the face and cell centre vertices
            };
            for ( const std::pair< const std::string, std::array< real, 4 > >& iPlane : tPlanes )
            {
                Corner_Values tCorners;
                for ( uint iCorner = 0; iCorner < 8; iCorner++ )
                {
                    Point tP            = corner_point( iCorner );
                    tCorners[ iCorner ] = iPlane.second[ 0 ] * tP[ 0 ] + iPlane.second[ 1 ] * tP[ 1 ] + iPlane.second[ 2 ] * tP[ 2 ] + iPlane.second[ 3 ];
                }
                tGeometries.push_back( { iPlane.first, tCorners } );
            }

            // one negative corner, and the trilinear saddle (sign alternates over the corners)
            tGeometries.push_back( { "corner", { -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 } } );
            tGeometries.push_back( { "saddle", { 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0 } } );
            return tGeometries;
        }

        // two-geometry sweep: geometry A over all corner combinations of aValues, geometry B from second_geometries(),
        // in both orders
        Sweep_Summary
        two_geometry_sweep( const std::string& aSweep, const Vector< real >& aValues, Regime aRegime )
        {
            Unit_Hex                     tHex;
            hmr::Interpolation_Mesh_HMR* tMesh = tHex.mMesh;

            Sweep_Summary tSummary;
            for ( const std::pair< std::string, Corner_Values >& iSecond : second_geometries() )
            {
                gCornersB         = iSecond.second;
                real tReferenceB  = reference_negative_volume( iSecond.second );
                for ( const Corner_Values& iCorners : all_corner_combinations( aValues ) )
                {
                    gCornersA        = iCorners;
                    real tReferenceA = reference_negative_volume( iCorners );

                    Case_Result tAB = decompose_case( tMesh, { create_geometry( &field_a, aRegime ), create_geometry( &field_b, aRegime ) } );
                    Case_Result tBA = decompose_case( tMesh, { create_geometry( &field_b, aRegime ), create_geometry( &field_a, aRegime ) } );
                    write_row( aSweep, aRegime, iSecond.first, "AB", iCorners, tAB, tReferenceA, tReferenceB );
                    write_row( aSweep, aRegime, iSecond.first, "BA", iCorners, tBA, tReferenceB, tReferenceA );

                    bool tAB_OK = accumulate_invariants( tSummary, tAB, iCorners );
                    bool tBA_OK = accumulate_invariants( tSummary, tBA, iCorners );
                    if ( not( tAB_OK and tBA_OK ) )
                    {
                        continue;
                    }

                    // the volumes of the four sign combinations must not depend on the order (STORED); in BA the bits
                    // of A and B are swapped
                    real tOrderDifference = 0.0;
                    for ( uint iBits = 0; iBits < 4; iBits++ )
                    {
                        uint tSwapped    = ( ( iBits & 1 ) << 1 ) | ( ( iBits >> 1 ) & 1 );
                        tOrderDifference = std::max( tOrderDifference, std::abs( tAB.mSignVolume( iBits ) - tBA.mSignVolume( tSwapped ) ) );
                    }
                    tSummary.mMaxOrderDifference = std::max( tSummary.mMaxOrderDifference, tOrderDifference );
                    if ( tOrderDifference > gVolumeTolerance )
                    {
                        tSummary.mNumOrder++;
                        tSummary.record_violation( "order dependence (second " + iSecond.first + ")", iCorners );
                    }

                    // STORED: the other geometry only subdivides, so each geometry keeps its own reference volume
                    if ( aRegime == Regime::STORED )
                    {
                        real tError = std::max( {
                                std::abs( tAB.mNegativeVolume( 0 ) - tReferenceA ),
                                std::abs( tAB.mNegativeVolume( 1 ) - tReferenceB ),
                                std::abs( tBA.mNegativeVolume( 0 ) - tReferenceB ),
                                std::abs( tBA.mNegativeVolume( 1 ) - tReferenceA ) } );
                        tSummary.mMaxReferenceError = std::max( tSummary.mMaxReferenceError, tError );
                        if ( tError > gVolumeTolerance )
                        {
                            tSummary.mNumReference++;
                            tSummary.record_violation( "reference volume (second " + iSecond.first + ")", iCorners );
                        }
                    }
                }
            }
            return tSummary;
        }

        void
        check_invariants( const Sweep_Summary& aSummary, bool aCheckReference, bool aCheckOrder )
        {
            CHECK( aSummary.mNumThrew == 0 );
            CHECK( aSummary.mNumNonPositive == 0 );
            CHECK( aSummary.mNumBadVolume == 0 );
            CHECK( aSummary.mNumNoPhase == 0 );
            CHECK( aSummary.mNumMixedSign == 0 );
            CHECK( aSummary.mNumBadFacets == 0 );
            if ( aCheckReference )
            {
                CHECK( aSummary.mNumReference == 0 );
                CHECK( aSummary.mNumMoment == 0 );
            }
            if ( aCheckOrder )
            {
                CHECK( aSummary.mNumOrder == 0 );
            }
        }

        // ----------------------------------------------------------------------------------
        // multi-plane sweeps: how cell count, subphases and interface facets evolve as cuts accumulate

        // signed distance to the plane through aPoint with normal aNormal, shifted by aShift along the normal
        Plane
        plane_through( const Point& aNormal, const Point& aPoint, real aShift )
        {
            real tNorm = std::sqrt( aNormal[ 0 ] * aNormal[ 0 ] + aNormal[ 1 ] * aNormal[ 1 ] + aNormal[ 2 ] * aNormal[ 2 ] );
            Plane tPlane;
            tPlane[ 3 ] = -aShift;
            for ( uint iDim = 0; iDim < 3; iDim++ )
            {
                tPlane[ iDim ] = aNormal[ iDim ] / tNorm;
                tPlane[ 3 ] -= tPlane[ iDim ] * aPoint[ iDim ];
            }
            return tPlane;
        }

        Corner_Values
        plane_corner_values( const Plane& aPlane )
        {
            Corner_Values tCorners;
            for ( uint iCorner = 0; iCorner < 8; iCorner++ )
            {
                Point tP            = corner_point( iCorner );
                tCorners[ iCorner ] = aPlane[ 0 ] * tP[ 0 ] + aPlane[ 1 ] * tP[ 1 ] + aPlane[ 2 ] * tP[ 2 ] + aPlane[ 3 ];
            }
            return tCorners;
        }

        Case_Result
        decompose_planes( hmr::Interpolation_Mesh_HMR* aMesh, const Vector< Plane >& aPlanes )
        {
            MORIS_ERROR( aPlanes.size() <= gMaxPlanes, "decompose_planes() - at most 6 planes" );
            Vector< std::shared_ptr< gen::Geometry > > tGeometries;
            for ( uint iPlane = 0; iPlane < aPlanes.size(); iPlane++ )
            {
                gPlanes[ iPlane ] = aPlanes( iPlane );
                tGeometries.push_back( create_geometry( gPlaneFields[ iPlane ], Regime::STORED ) );
            }
            return decompose_case( aMesh, tGeometries );
        }

        void
        print_metrics_header( const std::string& aTitle, const std::string& aExtraColumns = "" )
        {
            std::cout << "\n=== " << aTitle << " ===\n"
                      << std::setw( 26 ) << "case" << std::setw( 7 ) << "cells" << std::setw( 6 ) << "subph"
                      << std::setw( 13 ) << "min cell V" << std::setw( 6 ) << "V<=0" << std::setw( 13 ) << "min subph V"
                      << std::setw( 7 ) << "ifacet" << std::setw( 13 ) << "min ifacet A" << std::setw( 13 ) << "total V - 1"
                      << aExtraColumns << "\n";
        }

        void
        print_metrics_row( const std::string& aLabel, const Case_Result& aResult, const std::string& aExtra = "" )
        {
            std::cout << std::setw( 26 ) << aLabel;
            if ( aResult.mThrew )
            {
                std::cout << "  THREW: " << aResult.mMessage.substr( 0, 120 ) << "\n";
                return;
            }
            std::cout << std::setw( 7 ) << aResult.mNumCells << std::setw( 6 ) << aResult.mNumSubphases
                      << std::setw( 13 ) << std::setprecision( 3 ) << aResult.mMinCellVolume << std::setw( 6 ) << aResult.mNumNonPositive
                      << std::setw( 13 ) << aResult.mMinSubphaseVolume << std::setw( 7 ) << aResult.mNumInterfaceFacets
                      << std::setw( 13 ) << aResult.mMinInterfaceFacetArea << std::setw( 13 ) << aResult.mTotalVolume - 1.0
                      << aExtra << std::setprecision( 6 ) << "\n";
        }
    }    // namespace

    TEST_CASE( "XTK pathology sweep: single hex, corner signs", "[.][XTK][XTK_PATHOLOGY][XTK_PATHOLOGY_HEX_SIGNS]" )
    {
        if ( par_size() == 1 )
        {
            for ( Regime iRegime : { Regime::STORED, Regime::ANALYTIC, Regime::MULTILINEAR } )
            {
                Sweep_Summary tSummary = single_geometry_sweep( "signs", { -1.0, 0.0, 1.0 }, iRegime );
                print_summary( "single hex, corner values {-1,0,1}, " + regime_name( iRegime ), tSummary );
                check_invariants( tSummary, iRegime != Regime::MULTILINEAR, false );
            }
        }
    }

    TEST_CASE( "XTK pathology sweep: single hex, near-zero corner values", "[.][XTK][XTK_PATHOLOGY][XTK_PATHOLOGY_HEX_NEAR_ZERO]" )
    {
        if ( par_size() == 1 )
        {
            // 1e-13 lies inside the isocontour tolerance (1e-12, values label as INTERFACE), 1e-11 just outside it
            for ( Regime iRegime : { Regime::STORED, Regime::ANALYTIC, Regime::MULTILINEAR } )
            {
                for ( real tDelta : { 1.0e-13, 1.0e-11 } )
                {
                    Sweep_Summary      tSummary = single_geometry_sweep( "near_zero", { -1.0, -tDelta, tDelta, 1.0 }, iRegime );
                    std::ostringstream tTitle;
                    tTitle << "single hex, corner values {-1,-d,+d,1}, d = " << tDelta << ", " << regime_name( iRegime );
                    print_summary( tTitle.str(), tSummary );
                    check_invariants( tSummary, iRegime != Regime::MULTILINEAR, false );
                }
            }
        }
    }

    TEST_CASE( "XTK pathology sweep: single hex, two geometries in both orders", "[.][XTK][XTK_PATHOLOGY][XTK_PATHOLOGY_HEX_TWO_GEOMETRIES]" )
    {
        if ( par_size() == 1 )
        {
            for ( Regime iRegime : { Regime::STORED, Regime::ANALYTIC, Regime::MULTILINEAR } )
            {
                Sweep_Summary tSummary = two_geometry_sweep( "two_geometries", { -1.0, 0.0, 1.0 }, iRegime );
                print_summary( "single hex, geometry A corner values {-1,0,1} x 8 second geometries, both orders, " + regime_name( iRegime ), tSummary );
                check_invariants( tSummary, iRegime == Regime::STORED, iRegime == Regime::STORED );
            }
        }
    }

    // one case of the two-geometry sweep in a fresh process, for triage:
    //   MORIS_XTK_PATHOLOGY_ONE="<stored|analytic|multilinear>|<second geometry name>|<AB|BA>|c0;c1;...;c7"
    TEST_CASE( "XTK pathology: one two-geometry case from the environment", "[.][XTK][XTK_PATHOLOGY][XTK_PATHOLOGY_ONE]" )
    {
        if ( par_size() == 1 )
        {
            const char* tSpec = std::getenv( "MORIS_XTK_PATHOLOGY_ONE" );
            REQUIRE( tSpec != nullptr );

            std::stringstream          tStream( tSpec );
            std::string                tRegimeName, tSecondName, tOrder, tCornerList;
            std::getline( tStream, tRegimeName, '|' );
            std::getline( tStream, tSecondName, '|' );
            std::getline( tStream, tOrder, '|' );
            std::getline( tStream, tCornerList, '|' );

            Regime tRegime = tRegimeName == "stored" ? Regime::STORED : tRegimeName == "analytic" ? Regime::ANALYTIC : Regime::MULTILINEAR;

            std::stringstream tCornerStream( tCornerList );
            std::string       tValue;
            for ( uint iCorner = 0; iCorner < 8; iCorner++ )
            {
                REQUIRE( std::getline( tCornerStream, tValue, ';' ) );
                gCornersA[ iCorner ] = std::stod( tValue );
            }

            bool tFound = false;
            for ( const std::pair< std::string, Corner_Values >& iSecond : second_geometries() )
            {
                if ( iSecond.first == tSecondName )
                {
                    gCornersB = iSecond.second;
                    tFound    = true;
                }
            }
            REQUIRE( tFound );

            Unit_Hex tHex;
            Vector< std::shared_ptr< gen::Geometry > > tGeometries =
                    tOrder == "AB" ? Vector< std::shared_ptr< gen::Geometry > >{ create_geometry( &field_a, tRegime ), create_geometry( &field_b, tRegime ) }
                                   : Vector< std::shared_ptr< gen::Geometry > >{ create_geometry( &field_b, tRegime ), create_geometry( &field_a, tRegime ) };
            Case_Result tResult = decompose_case( tHex.mMesh, tGeometries );

            std::cout << "case " << tSpec << ": threw " << tResult.mThrew << " (" << tResult.mMessage.substr( 0, 160 ) << "), cells "
                      << tResult.mNumCells << ", total volume " << tResult.mTotalVolume << ", volume<=0 " << tResult.mNumNonPositive
                      << ", mixed " << tResult.mNumMixedSign << ", bad facets " << tResult.mNumBadFacets << "\n";
            CHECK( not tResult.mThrew );
        }
    }

    TEST_CASE( "XTK pathology sweep: single hex, 1 to 6 planes accumulating", "[.][XTK][XTK_PATHOLOGY][XTK_PATHOLOGY_HEX_LADDER]" )
    {
        if ( par_size() == 1 )
        {
            Unit_Hex                     tHex;
            hmr::Interpolation_Mesh_HMR* tMesh = tHex.mMesh;

            // six generic planes through points near the cell centre (off the grid and off the subdivision vertices)
            const Point tNormals[ 6 ] = { { 1.0, 0.31, 0.17 }, { 0.23, 1.0, -0.41 }, { -0.37, 0.19, 1.0 },
                { 1.0, -0.71, 0.53 }, { 0.61, 0.83, -0.29 }, { -0.47, 1.0, 0.67 } };
            const Point tPoints[ 6 ] = { { 0.47, 0.52, 0.49 }, { 0.53, 0.46, 0.51 }, { 0.50, 0.55, 0.45 },
                { 0.44, 0.50, 0.56 }, { 0.56, 0.48, 0.47 }, { 0.49, 0.44, 0.53 } };

            Sweep_Summary tSummary;
            print_metrics_header( "planes accumulating in one hex (first k of 6, then the same k in reverse order)" );
            for ( uint tK = 1; tK <= 6; tK++ )
            {
                Vector< Plane > tForward;
                for ( uint iPlane = 0; iPlane < tK; iPlane++ )
                {
                    tForward.push_back( plane_through( tNormals[ iPlane ], tPoints[ iPlane ], 0.0 ) );
                }
                Vector< Plane > tReverse;
                for ( uint iPlane = tK; iPlane > 0; iPlane-- )
                {
                    tReverse.push_back( tForward( iPlane - 1 ) );
                }

                Case_Result tForwardResult = decompose_planes( tMesh, tForward );
                Case_Result tReverseResult = decompose_planes( tMesh, tReverse );
                accumulate_invariants( tSummary, tForwardResult, {} );
                accumulate_invariants( tSummary, tReverseResult, {} );
                print_metrics_row( "k = " + std::to_string( tK ) + " forward", tForwardResult );
                print_metrics_row( "k = " + std::to_string( tK ) + " reverse", tReverseResult );
            }
            print_summary( "planes accumulating: invariants", tSummary );
            check_invariants( tSummary, false, false );
        }
    }

    TEST_CASE( "XTK pathology sweep: single hex, two parallel planes approaching", "[.][XTK][XTK_PATHOLOGY][XTK_PATHOLOGY_HEX_PARALLEL]" )
    {
        if ( par_size() == 1 )
        {
            Unit_Hex                     tHex;
            hmr::Interpolation_Mesh_HMR* tMesh = tHex.mMesh;

            // A = n.(x - p), B = A - delta: the layer { A >= 0, B < 0 } has thickness delta (h = 1)
            const Point tNormal = { 1.0, 0.3, 0.2 };
            const Point tPoint  = { 0.45, 0.5, 0.5 };

            Sweep_Summary tSummary;
            print_metrics_header( "two parallel planes a distance delta apart", "     layer V  reference V" );
            for ( real tDelta : { 1.0e-1, 1.0e-3, 1.0e-6, 1.0e-9, 1.0e-11, 1.0e-12, 1.0e-13, 1.0e-14 } )
            {
                Plane tA = plane_through( tNormal, tPoint, 0.0 );
                Plane tB = plane_through( tNormal, tPoint, tDelta );

                Case_Result tResult = decompose_planes( tMesh, { tA, tB } );
                accumulate_invariants( tSummary, tResult, {} );

                // layer: A positive (bit 0 set), B negative (bit 1 clear)
                real tReference = reference_negative_volume( plane_corner_values( tB ) ) - reference_negative_volume( plane_corner_values( tA ) );
                std::ostringstream tExtra;
                tExtra << std::setw( 13 ) << std::setprecision( 3 ) << ( tResult.mThrew ? 0.0 : tResult.mSignVolume( 1 ) )
                       << std::setw( 13 ) << tReference;
                std::ostringstream tLabel;
                tLabel << "delta = " << tDelta;
                print_metrics_row( tLabel.str(), tResult, tExtra.str() );
            }
            print_summary( "parallel planes: invariants", tSummary );
            check_invariants( tSummary, false, false );
        }
    }

    TEST_CASE( "XTK pathology sweep: single hex, two planes crossing at a shrinking angle", "[.][XTK][XTK_PATHOLOGY][XTK_PATHOLOGY_HEX_ANGLE]" )
    {
        if ( par_size() == 1 )
        {
            Unit_Hex                     tHex;
            hmr::Interpolation_Mesh_HMR* tMesh = tHex.mMesh;

            // two vertical planes through a point off the subdivision vertices, crossing at angle theta
            const Point tPoint = { 0.513, 0.507, 0.5 };
            const real  tAlpha = 0.3;

            Sweep_Summary tSummary;
            print_metrics_header( "two planes crossing at angle theta" );
            for ( real tThetaDegrees : { 90.0, 30.0, 5.0, 1.0, 0.1, 0.01 } )
            {
                real  tTheta = tThetaDegrees * M_PI / 180.0;
                Plane tA     = plane_through( { std::cos( tAlpha ), std::sin( tAlpha ), 0.0 }, tPoint, 0.0 );
                Plane tB     = plane_through( { std::cos( tAlpha + tTheta ), std::sin( tAlpha + tTheta ), 0.0 }, tPoint, 0.0 );

                Case_Result tResult = decompose_planes( tMesh, { tA, tB } );
                accumulate_invariants( tSummary, tResult, {} );

                std::ostringstream tLabel;
                tLabel << "theta = " << tThetaDegrees << " deg";
                print_metrics_row( tLabel.str(), tResult );
            }
            print_summary( "crossing planes: invariants", tSummary );
            check_invariants( tSummary, false, false );
        }
    }

    namespace
    {
        // ----------------------------------------------------------------------------------
        // multi-rank sweep: planes on, near and across the processor boundary of a 4 x 2 x 2 box [0,2] x [0,1] x [0,1]

        const real gBoxCellSize = 0.5;
        const real gBoxVolume   = 2.0;

        // integration vertices closer than this count as one location
        const real gLocationTolerance = 1.0e-12;

        // exact negative volume of a plane in the box: the 24-tet reference is exact for a linear field, per background cell
        real
        box_reference_negative_volume( const Plane& aPlane )
        {
            real tVolume = 0.0;
            for ( uint iX = 0; iX < 4; iX++ )
            {
                for ( uint iY = 0; iY < 2; iY++ )
                {
                    for ( uint iZ = 0; iZ < 2; iZ++ )
                    {
                        // the plane in the unit coordinates of the cell
                        Plane tLocal = { aPlane[ 0 ] * gBoxCellSize,
                            aPlane[ 1 ] * gBoxCellSize,
                            aPlane[ 2 ] * gBoxCellSize,
                            aPlane[ 0 ] * iX * gBoxCellSize + aPlane[ 1 ] * iY * gBoxCellSize + aPlane[ 2 ] * iZ * gBoxCellSize + aPlane[ 3 ] };
                        tVolume += reference_negative_volume( plane_corner_values( tLocal ) ) * gBoxCellSize * gBoxCellSize * gBoxCellSize;
                    }
                }
            }
            return tVolume;
        }

        // the owned part of the background mesh, as the boxes of the owned background cells; a point on a face shared by
        // two cells belongs to the cell above it, so every point of the domain belongs to exactly one rank
        struct Owned_Region
        {
            Vector< std::array< real, 6 > > mBoxes;

            explicit Owned_Region( mtk::Interpolation_Mesh* aMesh )
            {
                for ( uint iCell = 0; iCell < aMesh->get_num_elems(); iCell++ )
                {
                    mtk::Cell& tCell = aMesh->get_mtk_cell( iCell );
                    if ( tCell.get_owner() != (moris_id)par_rank() )
                    {
                        continue;
                    }
                    std::array< real, 6 > tBox = { MORIS_REAL_MAX, MORIS_REAL_MAX, MORIS_REAL_MAX, -MORIS_REAL_MAX, -MORIS_REAL_MAX, -MORIS_REAL_MAX };
                    for ( mtk::Vertex* iVertex : tCell.get_vertex_pointers() )
                    {
                        Matrix< DDRMat > tCoords = iVertex->get_coords();
                        for ( uint iDim = 0; iDim < 3; iDim++ )
                        {
                            tBox[ iDim ]     = std::min( tBox[ iDim ], tCoords( iDim ) );
                            tBox[ 3 + iDim ] = std::max( tBox[ 3 + iDim ], tCoords( iDim ) );
                        }
                    }
                    mBoxes.push_back( tBox );
                }
            }

            bool
            contains( const Point& aPoint ) const
            {
                const Point tDomainMax = { 2.0, 1.0, 1.0 };
                for ( const std::array< real, 6 >& iBox : mBoxes )
                {
                    bool tInside = true;
                    for ( uint iDim = 0; iDim < 3; iDim++ )
                    {
                        bool tBelowTop = aPoint[ iDim ] < iBox[ 3 + iDim ] or ( iBox[ 3 + iDim ] == tDomainMax[ iDim ] and aPoint[ iDim ] <= iBox[ 3 + iDim ] );
                        tInside        = tInside and aPoint[ iDim ] >= iBox[ iDim ] and tBelowTop;
                    }
                    if ( tInside )
                    {
                        return true;
                    }
                }
                return false;
            }
        };

        Point
        centroid( const Vector< mtk::Vertex* >& aVertices )
        {
            Point tCentroid = { 0.0, 0.0, 0.0 };
            for ( mtk::Vertex* iVertex : aVertices )
            {
                Matrix< DDRMat > tCoords = iVertex->get_coords();
                for ( uint iDim = 0; iDim < 3; iDim++ )
                {
                    tCentroid[ iDim ] += tCoords( iDim ) / aVertices.size();
                }
            }
            return tCentroid;
        }

        struct Parallel_Result
        {
            bool        mThrew   = false;    // on any rank
            std::string mMessage = "";       // first message on this rank

            // global quantities over the owned subphases and owned interface facets
            uint mNumCells           = 0;
            uint mNumSubphases       = 0;
            uint mNumNonPositive     = 0;
            uint mNumNoPhase         = 0;
            uint mNumMixedSign       = 0;
            uint mNumInterfaceFacets = 0;
            real mInterfaceArea      = 0.0;
            real mTotalVolume        = 0.0;
            real mMinCellVolume      = MORIS_REAL_MAX;
            real mMinSubphaseVolume  = MORIS_REAL_MAX;

            Vector< real > mNegativeVolume;

            // agreement of the ranks on the integration vertices (evaluated on rank 0)
            uint mNumVertexIds          = 0;    // distinct ids
            uint mNumSharedIds          = 0;    // ids held by more than one rank
            uint mNumUnassignedIds      = 0;    // vertices without a valid id
            uint mNumCoordMismatch      = 0;    // ids whose coordinates differ between ranks
            real mMaxCoordMismatch      = 0.0;
            uint mNumOwnerMismatch      = 0;    // ids whose owner differs between ranks
            uint mNumLabelMismatch      = 0;    // ids whose geometric region differs between ranks for some geometry
            uint mNumDuplicateLocations = 0;    // locations carrying more than one id
        };

        // agreement of the ranks on the integration vertices; aRows: one matrix per rank, rows ( id, x, y, z, owner, region per geometry )
        void
        check_vertex_agreement( const Vector< Matrix< DDRMat > >& aRows, uint aNumGeometries, Parallel_Result& aResult )
        {
            std::map< moris_id, Matrix< DDRMat > >                              tFirstRow;
            std::map< moris_id, std::set< uint > >                              tRanksOfId;
            std::set< moris_id >                                                tCoordMismatch, tOwnerMismatch, tLabelMismatch;
            std::map< std::tuple< long long, long long, long long >, std::set< moris_id > > tIdsAtLocation;

            for ( uint iRank = 0; iRank < aRows.size(); iRank++ )
            {
                const Matrix< DDRMat >& tRows = aRows( iRank );
                for ( uint iRow = 0; iRow < tRows.n_rows(); iRow++ )
                {
                    moris_id tId = (moris_id)tRows( iRow, 0 );
                    if ( tId < 0 or tId == MORIS_ID_MAX )
                    {
                        aResult.mNumUnassignedIds++;
                        continue;
                    }
                    tRanksOfId[ tId ].insert( iRank );
                    tIdsAtLocation[ { std::llround( tRows( iRow, 1 ) / gLocationTolerance ),
                                            std::llround( tRows( iRow, 2 ) / gLocationTolerance ),
                                            std::llround( tRows( iRow, 3 ) / gLocationTolerance ) } ]
                            .insert( tId );

                    std::map< moris_id, Matrix< DDRMat > >::iterator tFirst = tFirstRow.find( tId );
                    if ( tFirst == tFirstRow.end() )
                    {
                        tFirstRow[ tId ] = tRows.get_row( iRow );
                        continue;
                    }
                    const Matrix< DDRMat >& tReference = tFirst->second;
                    for ( uint iDim = 1; iDim <= 3; iDim++ )
                    {
                        real tDifference = std::abs( tRows( iRow, iDim ) - tReference( iDim ) );
                        if ( tDifference > 0.0 )
                        {
                            tCoordMismatch.insert( tId );
                            aResult.mMaxCoordMismatch = std::max( aResult.mMaxCoordMismatch, tDifference );
                        }
                    }
                    if ( tRows( iRow, 4 ) != tReference( 4 ) )
                    {
                        tOwnerMismatch.insert( tId );
                    }
                    for ( uint iGeom = 0; iGeom < aNumGeometries; iGeom++ )
                    {
                        if ( tRows( iRow, 5 + iGeom ) != tReference( 5 + iGeom ) )
                        {
                            tLabelMismatch.insert( tId );
                        }
                    }
                }
            }

            aResult.mNumVertexIds     = tRanksOfId.size();
            aResult.mNumCoordMismatch = tCoordMismatch.size();
            aResult.mNumOwnerMismatch = tOwnerMismatch.size();
            aResult.mNumLabelMismatch = tLabelMismatch.size();
            for ( const std::pair< const moris_id, std::set< uint > >& iId : tRanksOfId )
            {
                aResult.mNumSharedIds += iId.second.size() > 1;
            }
            for ( const std::pair< const std::tuple< long long, long long, long long >, std::set< moris_id > >& iLocation : tIdsAtLocation )
            {
                aResult.mNumDuplicateLocations += iLocation.second.size() > 1;
            }
        }

        // decompose on every rank and reduce; the geometries are planes in the given regime
        Parallel_Result
        decompose_parallel_case( mtk::Interpolation_Mesh* aMesh, const Vector< Plane >& aPlanes, Regime aRegime )
        {
            MORIS_ERROR( aPlanes.size() <= gMaxPlanes, "decompose_parallel_case() - at most 6 planes" );
            uint tNumGeometries = aPlanes.size();

            Parallel_Result tResult;
            tResult.mNegativeVolume.resize( tNumGeometries, 0.0 );

            Owned_Region     tOwned( aMesh );
            Vector< real >   tNegativeVolume( tNumGeometries, 0.0 );
            Matrix< DDRMat > tVertexRows( 0, 5 + tNumGeometries );
            bool             tThrew = false;

            // a rank that throws inside a collective call leaves the others waiting; the sweep then hangs (ctest timeout)
            try
            {
                Vector< std::shared_ptr< gen::Geometry > > tGeometries;
                for ( uint iPlane = 0; iPlane < tNumGeometries; iPlane++ )
                {
                    gPlanes[ iPlane ] = aPlanes( iPlane );
                    tGeometries.push_back( create_geometry( gPlaneFields[ iPlane ], aRegime ) );
                }

                gen::Geometry_Engine_Parameters tParameters;
                tParameters.mGeometries = tGeometries;
                gen::Geometry_Engine tGeometryEngine( aMesh, tParameters );

                Model tModel( 3, aMesh, &tGeometryEngine );
                tModel.mVerbose = false;
                tModel.decompose( { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 } );

                Cut_Integration_Mesh* tCutMesh   = tModel.get_cut_integration_mesh();
                moris_index           tNumPhases = (moris_index)tGeometryEngine.get_num_bulk_phase();

                // owned subphases: their background cell is owned by this rank
                for ( uint iSP = 0; iSP < tCutMesh->get_num_subphases(); iSP++ )
                {
                    if ( tCutMesh->get_subphase_parent_cell( iSP )->get_owner() != (moris_id)par_rank() )
                    {
                        continue;
                    }
                    tResult.mNumSubphases++;
                    real tSubphaseVolume = 0.0;
                    for ( mtk::Cell* iCell : tCutMesh->get_subphase_ig_cells( iSP )->mIgCellGroup )
                    {
                        real tVolume = iCell->compute_cell_measure();
                        tSubphaseVolume += tVolume;
                        tResult.mNumCells++;
                        tResult.mTotalVolume += tVolume;
                        tResult.mMinCellVolume = std::min( tResult.mMinCellVolume, tVolume );
                        tResult.mNumNonPositive += tVolume <= 0.0;

                        moris_index tBulkPhase = tCutMesh->get_cell_bulk_phase( iCell->get_index() );
                        tResult.mNumNoPhase += tBulkPhase < 0 or tBulkPhase >= tNumPhases;

                        Vector< mtk::Vertex* > tVertices = iCell->get_vertex_pointers();
                        bool                   tMixed    = false;
                        for ( uint iGeom = 0; iGeom < tNumGeometries; iGeom++ )
                        {
                            bool tHasPositive = false;
                            bool tHasNegative = false;
                            for ( mtk::Vertex* iVertex : tVertices )
                            {
                                gen::Geometric_Region tRegion = tGeometryEngine.get_geometric_region( iGeom, iVertex->get_index(), iVertex->get_coords() );
                                tHasPositive = tHasPositive or tRegion == gen::Geometric_Region::POSITIVE;
                                tHasNegative = tHasNegative or tRegion == gen::Geometric_Region::NEGATIVE;
                            }
                            tMixed = tMixed or ( tHasPositive and tHasNegative );
                            if ( tHasNegative and not tHasPositive )
                            {
                                tNegativeVolume( iGeom ) += tVolume;
                            }
                        }
                        tResult.mNumMixedSign += tMixed;
                    }
                    tResult.mMinSubphaseVolume = std::min( tResult.mMinSubphaseVolume, tSubphaseVolume );
                }

                // owned interface facets: their centroid lies in the owned region
                std::shared_ptr< Facet_Based_Connectivity > tFacets = tCutMesh->get_face_connectivity();
                if ( tFacets != nullptr )
                {
                    for ( uint iFacet = 0; iFacet < tFacets->mFacetToCell.size(); iFacet++ )
                    {
                        if ( tFacets->mFacetToCell( iFacet ).size() == 2
                                and tCutMesh->get_cell_bulk_phase( tFacets->mFacetToCell( iFacet )( 0 )->get_index() )
                                            != tCutMesh->get_cell_bulk_phase( tFacets->mFacetToCell( iFacet )( 1 )->get_index() )
                                and tOwned.contains( centroid( tFacets->mFacetVertices( iFacet ) ) ) )
                        {
                            tResult.mNumInterfaceFacets++;
                            tResult.mInterfaceArea += facet_area( tFacets->mFacetVertices( iFacet ) );
                        }
                    }
                }

                // every integration vertex this rank holds, owned or not
                uint tNumVertices = tCutMesh->get_num_entities( mtk::EntityRank::NODE, 0 );
                tVertexRows.set_size( tNumVertices, 5 + tNumGeometries );
                for ( uint iVertex = 0; iVertex < tNumVertices; iVertex++ )
                {
                    mtk::Vertex&     tVertex = tCutMesh->get_mtk_vertex( iVertex );
                    Matrix< DDRMat > tCoords = tVertex.get_coords();
                    tVertexRows( iVertex, 0 ) = tVertex.get_id();
                    tVertexRows( iVertex, 1 ) = tCoords( 0 );
                    tVertexRows( iVertex, 2 ) = tCoords( 1 );
                    tVertexRows( iVertex, 3 ) = tCoords( 2 );
                    tVertexRows( iVertex, 4 ) = tVertex.get_owner();
                    for ( uint iGeom = 0; iGeom < tNumGeometries; iGeom++ )
                    {
                        tVertexRows( iVertex, 5 + iGeom ) = (real)tGeometryEngine.get_geometric_region( iGeom, tVertex.get_index(), tCoords );
                    }
                }
            }
            catch ( const std::exception& aException )
            {
                tThrew           = true;
                tResult.mMessage = aException.what();
                tVertexRows.set_size( 0, 5 + tNumGeometries );
            }

            // reduce over the ranks
            tResult.mThrew              = sum_all( (uint)tThrew ) > 0;
            tResult.mNumCells           = sum_all( tResult.mNumCells );
            tResult.mNumSubphases       = sum_all( tResult.mNumSubphases );
            tResult.mNumNonPositive     = sum_all( tResult.mNumNonPositive );
            tResult.mNumNoPhase         = sum_all( tResult.mNumNoPhase );
            tResult.mNumMixedSign       = sum_all( tResult.mNumMixedSign );
            tResult.mNumInterfaceFacets = sum_all( tResult.mNumInterfaceFacets );
            tResult.mInterfaceArea      = sum_all( tResult.mInterfaceArea );
            tResult.mTotalVolume        = sum_all( tResult.mTotalVolume );
            tResult.mMinCellVolume      = min_all( tResult.mMinCellVolume );
            tResult.mMinSubphaseVolume  = min_all( tResult.mMinSubphaseVolume );
            for ( uint iGeom = 0; iGeom < tNumGeometries; iGeom++ )
            {
                tResult.mNegativeVolume( iGeom ) = sum_all( tNegativeVolume( iGeom ) );
            }

            Vector< Matrix< DDRMat > > tGatheredRows;
            gatherv_mats( tVertexRows, tGatheredRows, 0 );
            if ( par_rank() == 0 )
            {
                check_vertex_agreement( tGatheredRows, tNumGeometries, tResult );
            }
            return tResult;
        }

        struct Parallel_Case
        {
            std::string     mName;
            Vector< Plane > mPlanes;
        };

        // the processor boundary is x = 1 on 2 ranks; on 4 ranks HMR adds a second one
        Vector< Parallel_Case >
        parallel_cases()
        {
            const Point tFace     = { 1.0, 0.5, 0.5 };    // a vertex on x = 1 (an edge of four background cells)
            const Point tXAxis    = { 1.0, 0.0, 0.0 };
            const Point tSlanted  = { 1.0, 0.3, 0.2 };
            const Point tDiagonal = { 1.0, 1.0, 0.0 };

            Vector< Parallel_Case > tCases;
            tCases.push_back( { "x=3 (no cut, control)", { plane_through( tXAxis, { 3.0, 0.5, 0.5 }, 0.0 ) } } );
            tCases.push_back( { "x=0.3 (away from the face)", { plane_through( tXAxis, { 0.3, 0.5, 0.5 }, 0.0 ) } } );
            tCases.push_back( { "x=1 (processor face)", { plane_through( tXAxis, tFace, 0.0 ) } } );
            for ( real tOffset : { 1.0e-11, 1.0e-12, 1.0e-13, -1.0e-11, -1.0e-12 } )
            {
                std::ostringstream tName;
                tName << "x=1" << ( tOffset > 0.0 ? "+" : "" ) << tOffset;
                tCases.push_back( { tName.str(), { plane_through( tXAxis, tFace, tOffset ) } } );
            }
            tCases.push_back( { "x=0.75 (cell centres)", { plane_through( tXAxis, { 0.75, 0.5, 0.5 }, 0.0 ) } } );
            tCases.push_back( { "x=1.25 (cell centres)", { plane_through( tXAxis, { 1.25, 0.5, 0.5 }, 0.0 ) } } );
            tCases.push_back( { "slanted across the face", { plane_through( tSlanted, tFace, 0.0 ) } } );
            tCases.push_back( { "slanted +1e-11", { plane_through( tSlanted, tFace, 1.0e-11 ) } } );
            tCases.push_back( { "x+y=1.5 (face edge)", { plane_through( tDiagonal, tFace, 0.0 ) } } );
            tCases.push_back( { "x+y=1.5+1e-12", { plane_through( tDiagonal, tFace, 1.0e-12 ) } } );
            tCases.push_back( { "y=0.5+1e-12 (face nodes)", { plane_through( { 0.0, 1.0, 0.0 }, tFace, 1.0e-12 ) } } );

            // two parallel planes straddling the processor face
            for ( real tDelta : { 1.0e-3, 1.0e-6, 1.0e-9, 1.0e-11, 1.0e-12 } )
            {
                std::ostringstream tName;
                tName << "layer at x=1, d=" << tDelta;
                tCases.push_back( { tName.str(), { plane_through( tXAxis, tFace, -0.5 * tDelta ), plane_through( tXAxis, tFace, 0.5 * tDelta ) } } );
            }
            for ( real tDelta : { 1.0e-6, 1.0e-11 } )
            {
                std::ostringstream tName;
                tName << "slanted layer, d=" << tDelta;
                tCases.push_back( { tName.str(), { plane_through( tSlanted, tFace, -0.5 * tDelta ), plane_through( tSlanted, tFace, 0.5 * tDelta ) } } );
            }

            // two geometries crossing on the processor face
            tCases.push_back( { "x=1 and y=0.5", { plane_through( tXAxis, tFace, 0.0 ), plane_through( { 0.0, 1.0, 0.0 }, tFace, 0.0 ) } } );
            tCases.push_back( { "x=1+1e-11 and slanted", { plane_through( tXAxis, tFace, 1.0e-11 ), plane_through( tSlanted, tFace, 0.0 ) } } );
            return tCases;
        }

        // CSV output on rank 0, open only if MORIS_XTK_PATHOLOGY_MPI_CSV is set
        std::ofstream*
        parallel_csv_stream()
        {
            static std::ofstream sStream;
            static bool          sOpened = false;
            if ( not sOpened )
            {
                sOpened           = true;
                const char* tPath = std::getenv( "MORIS_XTK_PATHOLOGY_MPI_CSV" );
                if ( tPath != nullptr and par_rank() == 0 )
                {
                    sStream.open( tPath );
                    sStream << "num_ranks,regime,case,threw,num_cells,num_subphases,total_volume,neg_volume_a,ref_a,neg_volume_b,ref_b,"
                               "non_positive,no_phase,mixed_sign,interface_facets,interface_area,min_cell_volume,min_subphase_volume,"
                               "vertex_ids,shared_ids,unassigned_ids,coord_mismatch,max_coord_mismatch,owner_mismatch,label_mismatch,duplicate_locations\n";
                }
            }
            return sStream.is_open() ? &sStream : nullptr;
        }
    }    // namespace

    TEST_CASE( "XTK pathology sweep: planes at the processor boundary of a box, any number of ranks", "[.][XTK][XTK_PATHOLOGY][XTK_PATHOLOGY_MPI]" )
    {
        Grid_Mesh                tBox( { 4, 2, 2 }, { 2.0, 1.0, 1.0 } );
        mtk::Interpolation_Mesh* tMesh = tBox.mMesh;

        // where the processor boundaries are
        Owned_Region tOwned( tMesh );
        real         tOwnedMin[ 3 ] = { MORIS_REAL_MAX, MORIS_REAL_MAX, MORIS_REAL_MAX };
        real         tOwnedMax[ 3 ] = { -MORIS_REAL_MAX, -MORIS_REAL_MAX, -MORIS_REAL_MAX };
        for ( const std::array< real, 6 >& iBox : tOwned.mBoxes )
        {
            for ( uint iDim = 0; iDim < 3; iDim++ )
            {
                tOwnedMin[ iDim ] = std::min( tOwnedMin[ iDim ], iBox[ iDim ] );
                tOwnedMax[ iDim ] = std::max( tOwnedMax[ iDim ], iBox[ 3 + iDim ] );
            }
        }
        std::cout << "rank " << par_rank() << " of " << par_size() << " owns " << tOwned.mBoxes.size() << " background cells in ["
                  << tOwnedMin[ 0 ] << "," << tOwnedMax[ 0 ] << "] x [" << tOwnedMin[ 1 ] << "," << tOwnedMax[ 1 ] << "] x ["
                  << tOwnedMin[ 2 ] << "," << tOwnedMax[ 2 ] << "]\n"
                  << std::flush;

        uint tNumThrew = 0, tNumBadVolume = 0, tNumReference = 0, tNumCellDefects = 0, tNumDisagreements = 0;
        if ( par_rank() == 0 )
        {
            std::cout << "\n=== planes at the processor boundary, " << par_size() << " rank(s) ===\n"
                      << std::setw( 12 ) << "regime" << std::setw( 28 ) << "case" << std::setw( 7 ) << "cells" << std::setw( 6 ) << "subph"
                      << std::setw( 12 ) << "min cell V" << std::setw( 6 ) << "V<=0" << std::setw( 7 ) << "ifacet" << std::setw( 12 ) << "ifacet A"
                      << std::setw( 12 ) << "total V-2" << std::setw( 12 ) << "neg V err" << std::setw( 7 ) << "ids" << std::setw( 7 ) << "shared"
                      << std::setw( 7 ) << "coord" << std::setw( 7 ) << "owner" << std::setw( 7 ) << "label" << std::setw( 7 ) << "dupl" << "\n";
        }

        // MORIS_XTK_PATHOLOGY_MPI_CASE=<n> runs only the n-th case (0-based, in every regime), so a crash can be located
        const char* tOnlyCase = std::getenv( "MORIS_XTK_PATHOLOGY_MPI_CASE" );

        for ( Regime iRegime : { Regime::STORED, Regime::ANALYTIC, Regime::MULTILINEAR } )
        {
            Vector< Parallel_Case > tCases = parallel_cases();
            for ( uint iCaseIndex = 0; iCaseIndex < tCases.size(); iCaseIndex++ )
            {
                const Parallel_Case& iCase = tCases( iCaseIndex );
                if ( tOnlyCase != nullptr and std::atoi( tOnlyCase ) != (int)iCaseIndex )
                {
                    continue;
                }

                // two analytic/multilinear geometries can hit MB-040 (heap corruption) on develop; STORED only
                if ( iCase.mPlanes.size() > 1 and iRegime != Regime::STORED )
                {
                    continue;
                }

                Parallel_Result tResult = decompose_parallel_case( tMesh, iCase.mPlanes, iRegime );

                Vector< real > tReference;
                real           tReferenceError = 0.0;
                for ( uint iGeom = 0; iGeom < iCase.mPlanes.size(); iGeom++ )
                {
                    tReference.push_back( box_reference_negative_volume( iCase.mPlanes( iGeom ) ) );
                    tReferenceError = std::max( tReferenceError, std::abs( tResult.mNegativeVolume( iGeom ) - tReference( iGeom ) ) );
                }

                if ( par_rank() != 0 )
                {
                    continue;
                }

                tNumThrew += tResult.mThrew;
                if ( not tResult.mThrew )
                {
                    tNumBadVolume += std::abs( tResult.mTotalVolume - gBoxVolume ) > gVolumeTolerance;
                    tNumReference += tReferenceError > gVolumeTolerance;
                    tNumCellDefects += tResult.mNumNonPositive + tResult.mNumNoPhase + tResult.mNumMixedSign > 0;
                    tNumDisagreements += tResult.mNumUnassignedIds + tResult.mNumCoordMismatch + tResult.mNumOwnerMismatch + tResult.mNumLabelMismatch > 0;
                }

                std::cout << std::setw( 12 ) << regime_name( iRegime ) << std::setw( 28 ) << iCase.mName;
                if ( tResult.mThrew )
                {
                    std::cout << "  THREW: " << tResult.mMessage.substr( 0, 120 ) << "\n";
                }
                else
                {
                    std::cout << std::setprecision( 3 ) << std::setw( 7 ) << tResult.mNumCells << std::setw( 6 ) << tResult.mNumSubphases
                              << std::setw( 12 ) << tResult.mMinCellVolume << std::setw( 6 ) << tResult.mNumNonPositive
                              << std::setw( 7 ) << tResult.mNumInterfaceFacets << std::setw( 12 ) << tResult.mInterfaceArea
                              << std::setw( 12 ) << tResult.mTotalVolume - gBoxVolume << std::setw( 12 ) << tReferenceError
                              << std::setw( 7 ) << tResult.mNumVertexIds << std::setw( 7 ) << tResult.mNumSharedIds
                              << std::setw( 7 ) << tResult.mNumCoordMismatch << std::setw( 7 ) << tResult.mNumOwnerMismatch
                              << std::setw( 7 ) << tResult.mNumLabelMismatch << std::setw( 7 ) << tResult.mNumDuplicateLocations
                              << std::setprecision( 6 ) << "\n";
                }

                std::ofstream* tStream = parallel_csv_stream();
                if ( tStream != nullptr )
                {
                    std::string tName = iCase.mName;
                    std::replace( tName.begin(), tName.end(), ',', ' ' );
                    *tStream << std::setprecision( 17 ) << par_size() << "," << regime_name( iRegime ) << "," << tName << "," << tResult.mThrew << ","
                             << tResult.mNumCells << "," << tResult.mNumSubphases << "," << tResult.mTotalVolume << ","
                             << tResult.mNegativeVolume( 0 ) << "," << tReference( 0 ) << ","
                             << ( tReference.size() > 1 ? tResult.mNegativeVolume( 1 ) : 0.0 ) << "," << ( tReference.size() > 1 ? tReference( 1 ) : 0.0 ) << ","
                             << tResult.mNumNonPositive << "," << tResult.mNumNoPhase << "," << tResult.mNumMixedSign << ","
                             << tResult.mNumInterfaceFacets << "," << tResult.mInterfaceArea << "," << tResult.mMinCellVolume << ","
                             << tResult.mMinSubphaseVolume << "," << tResult.mNumVertexIds << "," << tResult.mNumSharedIds << ","
                             << tResult.mNumUnassignedIds << "," << tResult.mNumCoordMismatch << "," << tResult.mMaxCoordMismatch << ","
                             << tResult.mNumOwnerMismatch << "," << tResult.mNumLabelMismatch << "," << tResult.mNumDuplicateLocations << "\n"
                             << std::flush;
                }
            }
        }

        if ( par_rank() == 0 )
        {
            std::cout << "\n  threw " << tNumThrew << " | total volume != 2: " << tNumBadVolume << " | negative volume != exact: " << tNumReference
                      << " | cases with collapsed/no-phase/mixed-sign cells: " << tNumCellDefects
                      << " | cases where the ranks disagree on a vertex: " << tNumDisagreements << "\n";
            CHECK( tNumThrew == 0 );
            CHECK( tNumBadVolume == 0 );
            CHECK( tNumDisagreements == 0 );
        }
    }

    namespace
    {
        // ----------------------------------------------------------------------------------
        // aura sweep: a background cell held by several ranks (owner + aura copies) must be decomposed, numbered and
        // enriched identically on each of them, and the enriched meshes and output sets must not depend on the rank count

        // entities held by several ranks: column 0 is the global id, columns aFirst.. must agree between the copies;
        // returns the number of ids whose copies disagree and fills one example
        uint
        count_disagreeing_copies( const Vector< Matrix< DDRMat > >& aGathered, uint aFirst, uint aLast, std::string& aExample, uint& aNumShared )
        {
            std::map< moris_id, std::pair< uint, Matrix< DDRMat > > > tFirst;
            std::set< moris_id >                                      tShared, tBad;
            for ( uint iRank = 0; iRank < aGathered.size(); iRank++ )
            {
                for ( uint iRow = 0; iRow < aGathered( iRank ).n_rows(); iRow++ )
                {
                    moris_id tId = (moris_id)aGathered( iRank )( iRow, 0 );
                    std::map< moris_id, std::pair< uint, Matrix< DDRMat > > >::iterator tIt = tFirst.find( tId );
                    if ( tIt == tFirst.end() )
                    {
                        tFirst[ tId ] = { iRank, aGathered( iRank ).get_row( iRow ) };
                        continue;
                    }
                    tShared.insert( tId );
                    for ( uint iCol = aFirst; iCol <= aLast; iCol++ )
                    {
                        if ( aGathered( iRank )( iRow, iCol ) != tIt->second.second( iCol ) and tBad.insert( tId ).second and aExample.empty() )
                        {
                            std::ostringstream tStream;
                            tStream << "id " << tId << " column " << iCol << ": rank " << tIt->second.first << " has "
                                    << tIt->second.second( iCol ) << ", rank " << iRank << " has " << aGathered( iRank )( iRow, iCol );
                            aExample = tStream.str();
                        }
                    }
                }
            }
            aNumShared = tShared.size();
            return tBad.size();
        }

        // ids that two different entities carry (column 0 = id, column 1 = 1 if owned): owned ids must be unique
        uint
        count_duplicate_owned_ids( const Vector< Matrix< DDRMat > >& aGathered, uint& aNumOwned, uint& aNumInvalid )
        {
            std::map< moris_id, uint > tCount;
            aNumOwned   = 0;
            aNumInvalid = 0;
            for ( uint iRank = 0; iRank < aGathered.size(); iRank++ )
            {
                for ( uint iRow = 0; iRow < aGathered( iRank ).n_rows(); iRow++ )
                {
                    moris_id tId = (moris_id)aGathered( iRank )( iRow, 0 );
                    if ( aGathered( iRank )( iRow, 1 ) == 0.0 )
                    {
                        continue;
                    }
                    aNumOwned++;
                    if ( tId < 0 or tId == MORIS_ID_MAX )
                    {
                        aNumInvalid++;
                        continue;
                    }
                    tCount[ tId ]++;
                }
            }
            uint tDuplicates = 0;
            for ( const std::pair< const moris_id, uint >& iId : tCount )
            {
                tDuplicates += iId.second > 1;
            }
            return tDuplicates;
        }

        struct Aura_Result
        {
            bool        mThrew   = false;
            std::string mMessage = "";

            // copies of one entity on several ranks that disagree
            uint        mSharedBgCells = 0, mBgCellMismatch = 0;
            uint        mSharedIpVertices = 0, mIpVertexMismatch = 0;
            uint        mSharedIpCells = 0, mIpCellMismatch = 0;
            std::string mBgCellExample, mIpVertexExample, mIpCellExample;

            // informational: aura copies of an IP cell whose vertex ids differ from the owner's (rim vertices are local
            // copies), and vertices without a T-matrix
            uint        mIpCellIdMismatch = 0;
            std::string mIpCellIdExample;
            uint        mEmptyOwnedVertices = 0, mEmptyAuraVertices = 0;

            // vertices of IP cells used by a cluster on this rank that have no T-matrix: their contribution is dropped
            // silently in a release build (MSI Pdof_Host asserts only in debug)
            uint mClusterVerticesWithoutBasis = 0;

            // owned ids carried by two entities
            uint mOwnedIgCells = 0, mDuplicateIgCellIds = 0, mInvalidIgCellIds = 0;
            uint mOwnedSubphases = 0, mDuplicateSubphaseIds = 0, mInvalidSubphaseIds = 0;

            // rank-count invariants
            uint mNumEnrichedBasis = 0;    // distinct enriched basis ids referenced by owned enriched vertices
            uint mNumIpVertices    = 0;    // owned enriched interpolation vertices
            uint mNumIpCells       = 0;    // owned enriched interpolation cells

            // per output set: owned clusters and their measure
            std::map< std::string, std::pair< uint, real > > mSets;
        };

        Aura_Result
        aura_case( mtk::Interpolation_Mesh* aMesh, const Vector< Plane >& aPlanes, Regime aRegime )
        {
            Aura_Result tResult;
            uint        tRank = par_rank();

            Vector< Vector< real > > tBgRows, tIgCellRows, tSubphaseRows, tVertexRows, tCellRows, tBasisRows;
            std::map< std::string, std::pair< uint, real > > tLocalSets;
            bool tThrew = false;

            // a rank that throws inside a collective call leaves the others waiting; the sweep then hangs (timeout)
            try
            {
                Vector< std::shared_ptr< gen::Geometry > > tGeometries;
                for ( uint iPlane = 0; iPlane < aPlanes.size(); iPlane++ )
                {
                    gPlanes[ iPlane ] = aPlanes( iPlane );
                    tGeometries.push_back( create_geometry( gPlaneFields[ iPlane ], aRegime ) );
                }
                gen::Geometry_Engine_Parameters tParameters;
                tParameters.mGeometries = tGeometries;
                gen::Geometry_Engine tGeometryEngine( aMesh, tParameters );

                Model tModel( 3, aMesh, &tGeometryEngine );
                tModel.mVerbose = false;
                tModel.decompose( { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 } );
                tModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE, 0 );
                tModel.construct_face_oriented_ghost_penalization_cells();

                // 1. every background cell this rank holds: decomposition, ids and phases
                Cut_Integration_Mesh* tCutMesh = tModel.get_cut_integration_mesh();
                for ( uint iBg = 0; iBg < aMesh->get_num_elems(); iBg++ )
                {
                    mtk::Cell& tBgCell = aMesh->get_mtk_cell( iBg );

                    std::vector< std::vector< moris_id > > tCellVertexIds;
                    Vector< moris_id >             tCellIds;
                    Vector< std::array< moris_id, 3 > > tSubphases;
                    for ( moris_index iSP : tCutMesh->get_parent_cell_subphases( iBg ) )
                    {
                        Vector< mtk::Cell* >& tCells = tCutMesh->get_subphase_ig_cells( iSP )->mIgCellGroup;
                        tSubphases.push_back( { tCutMesh->get_subphase_id( iSP ), tCutMesh->get_subphase_bulk_phase( iSP ), (moris_id)tCells.size() } );
                        bool tOwned = tBgCell.get_owner() == (moris_id)tRank;
                        tSubphaseRows.push_back( { (real)tCutMesh->get_subphase_id( iSP ), (real)tOwned } );
                        for ( mtk::Cell* iCell : tCells )
                        {
                            tCellIds.push_back( iCell->get_id() );
                            tIgCellRows.push_back( { (real)iCell->get_id(), (real)tOwned } );
                            std::vector< moris_id > tIds;
                            for ( mtk::Vertex* iVertex : iCell->get_vertex_pointers() )
                            {
                                tIds.push_back( iVertex->get_id() );
                            }
                            std::sort( tIds.begin(), tIds.end() );
                            tCellVertexIds.push_back( tIds );
                        }
                    }
                    std::sort( tCellVertexIds.begin(), tCellVertexIds.end() );
                    std::sort( tCellIds.begin(), tCellIds.end() );
                    std::sort( tSubphases.begin(), tSubphases.end() );

                    Hash tGeometryHash, tIdHash, tSubphaseHash;
                    for ( const std::vector< moris_id >& iIds : tCellVertexIds )
                    {
                        for ( moris_id iId : iIds )
                        {
                            tGeometryHash.add( iId );
                        }
                        tGeometryHash.add( -1 );
                    }
                    for ( moris_id iId : tCellIds )
                    {
                        tIdHash.add( iId );
                    }
                    for ( const std::array< moris_id, 3 >& iSP : tSubphases )
                    {
                        tSubphaseHash.add( iSP[ 0 ] );
                        tSubphaseHash.add( iSP[ 1 ] );
                        tSubphaseHash.add( iSP[ 2 ] );
                    }
                    // id | owner, #subphases, #IG cells, connectivity, IG cell ids, (subphase id, bulk phase, #cells)
                    tBgRows.push_back( { (real)tBgCell.get_id(), (real)tBgCell.get_owner(), (real)tSubphases.size(), (real)tCellIds.size(),
                            tGeometryHash.high(), tGeometryHash.low(), tIdHash.high(), tIdHash.low(), tSubphaseHash.high(), tSubphaseHash.low() } );
                }

                // 2. enriched interpolation mesh: vertices with their enriched basis, cells with their vertices
                Enriched_Interpolation_Mesh& tEnrIpMesh = tModel.get_enriched_interp_mesh( 0 );
                for ( uint iV = 0; iV < tEnrIpMesh.get_num_entities( mtk::EntityRank::NODE, 0 ); iV++ )
                {
                    mtk::Vertex&               tVertex        = tEnrIpMesh.get_mtk_vertex( iV );
                    mtk::Vertex_Interpolation* tInterpolation = tVertex.get_interpolation( 0 );
                    Matrix< IdMat >            tBasisIds      = tInterpolation->get_ids();
                    Matrix< IdMat >            tBasisOwners   = tInterpolation->get_owners();
                    const Matrix< DDRMat >&    tWeights       = *tInterpolation->get_weights();

                    Vector< std::array< std::int64_t, 3 > > tBasis;
                    for ( uint iB = 0; iB < tBasisIds.numel(); iB++ )
                    {
                        tBasis.push_back( { (std::int64_t)tBasisIds( iB ), (std::int64_t)tBasisOwners( iB ), std::llround( tWeights( iB ) * 1.0e12 ) } );
                    }
                    std::sort( tBasis.begin(), tBasis.end() );
                    Hash tBasisHash;
                    for ( const std::array< std::int64_t, 3 >& iB : tBasis )
                    {
                        tBasisHash.add( iB[ 0 ] );
                        tBasisHash.add( iB[ 1 ] );
                        tBasisHash.add( iB[ 2 ] );
                    }
                    Matrix< DDRMat > tCoords = tVertex.get_coords();
                    bool             tOwned  = tVertex.get_owner() == (moris_id)tRank;

                    // aura vertices get their T-matrix only when a cluster needs it (communicate_T_matrices_for_UIPCs,
                    // ghost): an empty basis is expected there and is counted, not compared
                    if ( tBasis.size() == 0 )
                    {
                        ( tOwned ? tResult.mEmptyOwnedVertices : tResult.mEmptyAuraVertices )++;
                        continue;
                    }

                    // id | owner, x, y, z, #basis, (basis id, basis owner, weight)
                    tVertexRows.push_back( { (real)tVertex.get_id(), (real)tVertex.get_owner(), tCoords( 0 ), tCoords( 1 ), tCoords( 2 ),
                            (real)tBasis.size(), tBasisHash.high(), tBasisHash.low(), (real)tOwned } );
                    if ( tOwned )
                    {
                        for ( const std::array< std::int64_t, 3 >& iB : tBasis )
                        {
                            tBasisRows.push_back( { (real)iB[ 0 ] } );
                        }
                    }
                }
                for ( uint iC = 0; iC < tEnrIpMesh.get_num_entities( mtk::EntityRank::ELEMENT, 0 ); iC++ )
                {
                    mtk::Cell& tCell = tEnrIpMesh.get_mtk_cell( iC );

                    // vertex locations (must agree) and vertex ids (aura rim vertices are local copies: reported only)
                    Vector< std::array< std::int64_t, 3 > > tLocations;
                    Vector< moris_id >                      tIds;
                    for ( mtk::Vertex* iVertex : tCell.get_vertex_pointers() )
                    {
                        Matrix< DDRMat > tX = iVertex->get_coords();
                        tLocations.push_back( { std::llround( tX( 0 ) / gLocationTolerance ), std::llround( tX( 1 ) / gLocationTolerance ),
                                std::llround( tX( 2 ) / gLocationTolerance ) } );
                        tIds.push_back( iVertex->get_id() );
                    }
                    std::sort( tLocations.begin(), tLocations.end() );
                    std::sort( tIds.begin(), tIds.end() );
                    Hash tLocationHash, tIdHash;
                    for ( const std::array< std::int64_t, 3 >& iL : tLocations )
                    {
                        tLocationHash.add( iL[ 0 ] );
                        tLocationHash.add( iL[ 1 ] );
                        tLocationHash.add( iL[ 2 ] );
                    }
                    for ( moris_id iId : tIds )
                    {
                        tIdHash.add( iId );
                    }
                    // id | owner, vertex locations | vertex ids, owned
                    tCellRows.push_back( { (real)tCell.get_id(), (real)tCell.get_owner(), tLocationHash.high(), tLocationHash.low(),
                            tIdHash.high(), tIdHash.low(), (real)( tCell.get_owner() == (moris_id)tRank ) } );
                }

                // MORIS_XTK_PATHOLOGY_AURA_DUMP=<prefix>: write this rank's enriched interpolation cells and vertices
                if ( const char* tDump = std::getenv( "MORIS_XTK_PATHOLOGY_AURA_DUMP" ) )
                {
                    std::ofstream tStream( std::string( tDump ) + "_rank" + std::to_string( tRank ) + ".txt" );
                    tStream << std::setprecision( 6 );
                    for ( uint iC = 0; iC < tEnrIpMesh.get_num_entities( mtk::EntityRank::ELEMENT, 0 ); iC++ )
                    {
                        mtk::Cell& tCell = tEnrIpMesh.get_mtk_cell( iC );
                        Point      tC    = centroid( tCell.get_vertex_pointers() );
                        tStream << "cell index " << iC << " id " << tCell.get_id() << " owner " << tCell.get_owner() << " centroid " << tC[ 0 ] << " "
                                << tC[ 1 ] << " " << tC[ 2 ] << " vertices";
                        for ( mtk::Vertex* iVertex : tCell.get_vertex_pointers() )
                        {
                            tStream << " " << iVertex->get_id();
                        }
                        tStream << "\n";
                    }
                    for ( uint iV = 0; iV < tEnrIpMesh.get_num_entities( mtk::EntityRank::NODE, 0 ); iV++ )
                    {
                        mtk::Vertex&     tVertex = tEnrIpMesh.get_mtk_vertex( iV );
                        Matrix< DDRMat > tX      = tVertex.get_coords();
                        Matrix< IdMat >  tIds    = tVertex.get_interpolation( 0 )->get_ids();
                        tStream << "vertex index " << iV << " id " << tVertex.get_id() << " owner " << tVertex.get_owner() << " at " << tX( 0 ) << " "
                                << tX( 1 ) << " " << tX( 2 ) << " basis";
                        for ( uint iB = 0; iB < tIds.numel(); iB++ )
                        {
                            tStream << " " << tIds( iB );
                        }
                        tStream << "\n";
                    }
                }

                // 3. output sets of the enriched integration mesh: every cluster this rank holds, its measure, and its IP vertices
                Enriched_Integration_Mesh& tEnrIgMesh = tModel.get_enriched_integ_mesh( 0 );
                for ( uint iSet = 0; iSet < tEnrIgMesh.get_num_sets(); iSet++ )
                {
                    mtk::Set*   tSet  = tEnrIgMesh.get_set_by_index( iSet );
                    std::string tName = tSet->get_set_name();
                    std::pair< uint, real >& tEntry = tLocalSets[ tName ];
                    uint tNumAll = 0, tNumLeaderOwned = 0, tNumFollowerOwned = 0;
                    for ( const mtk::Cluster* iCluster : tSet->get_clusters_on_set() )
                    {
                        tNumAll++;
                        bool tLeaderOwned = iCluster->get_interpolation_cell( mtk::Leader_Follower::LEADER ).get_owner() == (moris_id)tRank;
                        tNumLeaderOwned += tLeaderOwned;
                        tNumFollowerOwned += tSet->get_set_type() == mtk::SetType::DOUBLE_SIDED_SIDESET
                                         and iCluster->get_interpolation_cell( mtk::Leader_Follower::FOLLOWER ).get_owner() == (moris_id)tRank;
                    }
                    if ( std::getenv( "MORIS_XTK_PATHOLOGY_AURA_DUMP" ) != nullptr and tNumAll > 0 )
                    {
                        std::cout << "rank " << tRank << " set " << tName << ": clusters " << tNumAll << ", leader owned " << tNumLeaderOwned
                                  << ", follower owned " << tNumFollowerOwned << "\n";
                    }
                    // FEM assembles every cluster a rank holds (no ownership filter, cl_FEM_Set.cpp:266-315): count them all
                    for ( const mtk::Cluster* iCluster : tSet->get_clusters_on_set() )
                    {
                        Vector< mtk::Leader_Follower > tSides = { mtk::Leader_Follower::LEADER };
                        if ( tSet->get_set_type() == mtk::SetType::DOUBLE_SIDED_SIDESET )
                        {
                            tSides.push_back( mtk::Leader_Follower::FOLLOWER );
                        }
                        for ( mtk::Leader_Follower iSide : tSides )
                        {
                            for ( mtk::Vertex* iVertex : iCluster->get_interpolation_cell( iSide ).get_vertex_pointers() )
                            {
                                tResult.mClusterVerticesWithoutBasis += iVertex->get_interpolation( 0 )->get_ids().numel() == 0;
                            }
                        }
                        tEntry.first++;
                        tEntry.second += tSet->get_set_type() == mtk::SetType::BULK
                                               ? iCluster->compute_cluster_cell_measure( mtk::Primary_Void::PRIMARY, mtk::Leader_Follower::LEADER )
                                               : iCluster->compute_cluster_cell_side_measure( mtk::Primary_Void::PRIMARY, mtk::Leader_Follower::LEADER );
                    }
                }
            }
            catch ( const std::exception& aException )
            {
                tThrew           = true;
                tResult.mMessage = aException.what();
            }
            tResult.mThrew                       = sum_all( (uint)tThrew ) > 0;
            tResult.mEmptyOwnedVertices          = sum_all( tResult.mEmptyOwnedVertices );
            tResult.mEmptyAuraVertices           = sum_all( tResult.mEmptyAuraVertices );
            tResult.mClusterVerticesWithoutBasis = sum_all( tResult.mClusterVerticesWithoutBasis );

            Vector< Matrix< DDRMat > > tBg       = gather_rows( tBgRows, 10 );
            Vector< Matrix< DDRMat > > tIgCells  = gather_rows( tIgCellRows, 2 );
            Vector< Matrix< DDRMat > > tSps      = gather_rows( tSubphaseRows, 2 );
            Vector< Matrix< DDRMat > > tVertices = gather_rows( tVertexRows, 9 );
            Vector< Matrix< DDRMat > > tCells    = gather_rows( tCellRows, 7 );
            Vector< Matrix< DDRMat > > tBasis    = gather_rows( tBasisRows, 1 );

            // set names may differ between ranks (a rank without a cut has no interface set): gather names + values
            Vector< Vector< real > > tSetRows;
            Vector< std::string >    tSetNames;
            for ( const std::pair< const std::string, std::pair< uint, real > >& iSet : tLocalSets )
            {
                Hash tNameHash;
                for ( char iChar : iSet.first )
                {
                    tNameHash.add( iChar );
                }
                tSetRows.push_back( { tNameHash.high(), tNameHash.low(), (real)iSet.second.first, iSet.second.second } );
                tSetNames.push_back( iSet.first );
            }
            Vector< Matrix< DDRMat > > tSets = gather_rows( tSetRows, 4 );

            if ( tRank == 0 )
            {
                tResult.mBgCellMismatch   = count_disagreeing_copies( tBg, 1, 9, tResult.mBgCellExample, tResult.mSharedBgCells );
                tResult.mIpVertexMismatch = count_disagreeing_copies( tVertices, 1, 7, tResult.mIpVertexExample, tResult.mSharedIpVertices );
                tResult.mIpCellMismatch   = count_disagreeing_copies( tCells, 1, 3, tResult.mIpCellExample, tResult.mSharedIpCells );
                uint tSharedForIds        = 0;
                tResult.mIpCellIdMismatch = count_disagreeing_copies( tCells, 4, 5, tResult.mIpCellIdExample, tSharedForIds );

                tResult.mDuplicateIgCellIds   = count_duplicate_owned_ids( tIgCells, tResult.mOwnedIgCells, tResult.mInvalidIgCellIds );
                tResult.mDuplicateSubphaseIds = count_duplicate_owned_ids( tSps, tResult.mOwnedSubphases, tResult.mInvalidSubphaseIds );

                std::set< moris_id > tBasisIds, tVertexIds, tCellIds;
                for ( const Matrix< DDRMat >& iRows : tBasis )
                {
                    for ( uint iRow = 0; iRow < iRows.n_rows(); iRow++ )
                    {
                        tBasisIds.insert( (moris_id)iRows( iRow, 0 ) );
                    }
                }
                for ( const Matrix< DDRMat >& iRows : tVertices )
                {
                    for ( uint iRow = 0; iRow < iRows.n_rows(); iRow++ )
                    {
                        if ( iRows( iRow, 8 ) != 0.0 )
                        {
                            tVertexIds.insert( (moris_id)iRows( iRow, 0 ) );
                        }
                    }
                }
                for ( const Matrix< DDRMat >& iRows : tCells )
                {
                    for ( uint iRow = 0; iRow < iRows.n_rows(); iRow++ )
                    {
                        if ( iRows( iRow, 6 ) != 0.0 )
                        {
                            tCellIds.insert( (moris_id)iRows( iRow, 0 ) );
                        }
                    }
                }
                tResult.mNumEnrichedBasis = tBasisIds.size();
                tResult.mNumIpVertices    = tVertexIds.size();
                tResult.mNumIpCells       = tCellIds.size();

                // rank 0 knows its own names; other ranks' names arrive as hashes: map them through rank 0's names, and
                // keep unknown names as "set#<hash>"
                std::map< std::pair< real, real >, std::string > tNameOfHash;
                for ( uint iSet = 0; iSet < tSetNames.size(); iSet++ )
                {
                    tNameOfHash[ { tSetRows( iSet )( 0 ), tSetRows( iSet )( 1 ) } ] = tSetNames( iSet );
                }
                for ( const Matrix< DDRMat >& iRows : tSets )
                {
                    for ( uint iRow = 0; iRow < iRows.n_rows(); iRow++ )
                    {
                        std::pair< real, real > tKey = { iRows( iRow, 0 ), iRows( iRow, 1 ) };
                        std::string             tName;
                        if ( tNameOfHash.count( tKey ) )
                        {
                            tName = tNameOfHash[ tKey ];
                        }
                        else
                        {
                            std::ostringstream tStream;
                            tStream << "set#" << std::setprecision( 12 ) << tKey.first << "_" << tKey.second;
                            tName = tStream.str();
                        }
                        tResult.mSets[ tName ].first += (uint)iRows( iRow, 2 );
                        tResult.mSets[ tName ].second += iRows( iRow, 3 );
                    }
                }
            }
            return tResult;
        }

        std::ofstream*
        aura_csv_stream( const char* aVariable, const std::string& aHeader )
        {
            static std::map< std::string, std::ofstream > sStreams;
            const char* tPath = std::getenv( aVariable );
            if ( tPath == nullptr or par_rank() != 0 )
            {
                return nullptr;
            }
            std::ofstream& tStream = sStreams[ aVariable ];
            if ( not tStream.is_open() )
            {
                tStream.open( tPath );
                tStream << aHeader << "\n";
            }
            return &tStream;
        }
    }    // namespace

    TEST_CASE( "XTK pathology sweep: aura copies, enrichment and output sets agree across ranks", "[.][XTK][XTK_PATHOLOGY][XTK_PATHOLOGY_MPI_AURA]" )
    {
        Grid_Mesh                tBox( { 4, 2, 2 }, { 2.0, 1.0, 1.0 } );
        mtk::Interpolation_Mesh* tMesh = tBox.mMesh;

        const char* tOnlyCase = std::getenv( "MORIS_XTK_PATHOLOGY_MPI_CASE" );

        uint tNumThrew = 0, tNumCopyMismatch = 0, tNumDuplicateIds = 0;
        if ( par_rank() == 0 )
        {
            std::cout << "\n=== aura copies, enrichment and output sets, " << par_size() << " rank(s) ===\n";
        }

        for ( Regime iRegime : { Regime::STORED, Regime::ANALYTIC } )
        {
            Vector< Parallel_Case > tCases = parallel_cases();
            for ( uint iCaseIndex = 0; iCaseIndex < tCases.size(); iCaseIndex++ )
            {
                const Parallel_Case& iCase = tCases( iCaseIndex );
                if ( ( tOnlyCase != nullptr and std::atoi( tOnlyCase ) != (int)iCaseIndex ) or ( iCase.mPlanes.size() > 1 and iRegime != Regime::STORED ) )
                {
                    continue;
                }

                Aura_Result tResult = aura_case( tMesh, iCase.mPlanes, iRegime );
                if ( par_rank() != 0 )
                {
                    continue;
                }

                tNumThrew += tResult.mThrew;
                uint tMismatch  = tResult.mBgCellMismatch + tResult.mIpVertexMismatch + tResult.mIpCellMismatch + tResult.mClusterVerticesWithoutBasis
                               + tResult.mEmptyOwnedVertices;
                uint tDuplicate = tResult.mDuplicateIgCellIds + tResult.mInvalidIgCellIds + tResult.mDuplicateSubphaseIds + tResult.mInvalidSubphaseIds;
                tNumCopyMismatch += tMismatch > 0;
                tNumDuplicateIds += tDuplicate > 0;

                std::cout << std::setw( 9 ) << regime_name( iRegime ) << std::setw( 28 ) << iCase.mName;
                if ( tResult.mThrew )
                {
                    std::cout << "  THREW: " << tResult.mMessage.substr( 0, 120 ) << "\n";
                    continue;
                }
                std::cout << "  bg " << tResult.mBgCellMismatch << "/" << tResult.mSharedBgCells
                          << "  ipv " << tResult.mIpVertexMismatch << "/" << tResult.mSharedIpVertices
                          << "  ipc " << tResult.mIpCellMismatch << "/" << tResult.mSharedIpCells
                          << "  dup ig " << tResult.mDuplicateIgCellIds << "+" << tResult.mInvalidIgCellIds << "/" << tResult.mOwnedIgCells
                          << "  dup sp " << tResult.mDuplicateSubphaseIds << "+" << tResult.mInvalidSubphaseIds << "/" << tResult.mOwnedSubphases
                          << "  basis " << tResult.mNumEnrichedBasis << "  ipv " << tResult.mNumIpVertices << "  ipc " << tResult.mNumIpCells
                          << "  cluster vertices w/o T " << tResult.mClusterVerticesWithoutBasis << "  empty T owned/aura " << tResult.mEmptyOwnedVertices
                          << "/" << tResult.mEmptyAuraVertices << "  ipc id diff " << tResult.mIpCellIdMismatch << "\n";
                for ( const std::string& iExample : { tResult.mBgCellExample, tResult.mIpVertexExample, tResult.mIpCellExample } )
                {
                    if ( not iExample.empty() )
                    {
                        std::cout << "      first mismatch: " << iExample << "\n";
                    }
                }

                std::string tName = iCase.mName;
                std::replace( tName.begin(), tName.end(), ',', ' ' );
                if ( std::ofstream* tStream = aura_csv_stream( "MORIS_XTK_PATHOLOGY_AURA_CSV",
                             "num_ranks,regime,case,threw,shared_bg,bg_mismatch,shared_ipv,ipv_mismatch,shared_ipc,ipc_mismatch,"
                             "owned_ig,dup_ig,invalid_ig,owned_sp,dup_sp,invalid_sp,enriched_basis,ip_vertices,ip_cells,cluster_vertices_without_t,"
                             "empty_t_owned,empty_t_aura,ipc_id_mismatch,example" ) )
                {
                    std::string tExample = tResult.mBgCellExample + " " + tResult.mIpVertexExample + " " + tResult.mIpCellExample;
                    std::replace( tExample.begin(), tExample.end(), ',', ' ' );
                    *tStream << par_size() << "," << regime_name( iRegime ) << "," << tName << "," << tResult.mThrew << ","
                             << tResult.mSharedBgCells << "," << tResult.mBgCellMismatch << "," << tResult.mSharedIpVertices << ","
                             << tResult.mIpVertexMismatch << "," << tResult.mSharedIpCells << "," << tResult.mIpCellMismatch << ","
                             << tResult.mOwnedIgCells << "," << tResult.mDuplicateIgCellIds << "," << tResult.mInvalidIgCellIds << ","
                             << tResult.mOwnedSubphases << "," << tResult.mDuplicateSubphaseIds << "," << tResult.mInvalidSubphaseIds << ","
                             << tResult.mNumEnrichedBasis << "," << tResult.mNumIpVertices << "," << tResult.mNumIpCells << ","
                             << tResult.mClusterVerticesWithoutBasis << "," << tResult.mEmptyOwnedVertices << "," << tResult.mEmptyAuraVertices << ","
                             << tResult.mIpCellIdMismatch << "," << tExample << "\n"
                             << std::flush;
                }
                if ( std::ofstream* tStream = aura_csv_stream( "MORIS_XTK_PATHOLOGY_SETS_CSV", "num_ranks,regime,case,set,owned_clusters,measure" ) )
                {
                    for ( const std::pair< const std::string, std::pair< uint, real > >& iSet : tResult.mSets )
                    {
                        *tStream << par_size() << "," << regime_name( iRegime ) << "," << tName << "," << iSet.first << ","
                                 << iSet.second.first << "," << std::setprecision( 17 ) << iSet.second.second << std::setprecision( 6 ) << "\n";
                    }
                    *tStream << std::flush;
                }
            }
        }

        if ( par_rank() == 0 )
        {
            std::cout << "\n  threw " << tNumThrew << " | cases with disagreeing copies: " << tNumCopyMismatch
                      << " | cases with duplicate or invalid owned ids: " << tNumDuplicateIds << "\n";
            CHECK( tNumThrew == 0 );
            CHECK( tNumCopyMismatch == 0 );
            CHECK( tNumDuplicateIds == 0 );
        }
    }

    TEST_CASE( "XTK pathology: intersected cells on different refinement levels on different ranks", "[.][XTK][XTK_PATHOLOGY][XTK_PATHOLOGY_MPI_LEVELS]" )
    {
        // finding 6: check_intersected_background_cell_levels compares each rank's intersected cells with that rank's
        // first one and combines with all_land, so ranks that are each uniform at different levels pass. The box is
        // refined one level for x > 1.5; the processor face is x = 1 on 2 ranks, and each rank's aura reaches 0.5 past it.
        Grid_Mesh                tBox( { 4, 2, 2 }, { 2.0, 1.0, 1.0 }, 1.5 );
        mtk::Interpolation_Mesh* tMesh = tBox.mMesh;

        struct Level_Case
        {
            std::string     mName;
            Vector< Plane > mPlanes;
        };
        const Point          tXAxis = { 1.0, 0.0, 0.0 };
        Vector< Level_Case > tCases = {
            { "x=1.8 only (fine cells)", { plane_through( tXAxis, { 1.8, 0.5, 0.5 }, 0.0 ) } },
            { "x=0.3 only (coarse cells)", { plane_through( tXAxis, { 0.3, 0.5, 0.5 }, 0.0 ) } },
            { "x=0.3 and x=1.8", { plane_through( tXAxis, { 0.3, 0.5, 0.5 }, 0.0 ), plane_through( tXAxis, { 1.8, 0.5, 0.5 }, 0.0 ) } },
            { "slanted across x=1.5", { plane_through( { 1.0, 0.3, 0.2 }, { 1.5, 0.5, 0.5 }, 0.0 ) } }
        };

        uint tNumInconsistent = 0;
        for ( const Level_Case& iCase : tCases )
        {
            // levels of the cells each plane cuts on this rank (corner signs), owned and aura
            std::set< uint > tLevels;
            for ( uint iCell = 0; iCell < tMesh->get_num_elems(); iCell++ )
            {
                mtk::Cell& tCell = tMesh->get_mtk_cell( iCell );
                for ( const Plane& iPlane : iCase.mPlanes )
                {
                    bool tNegative = false, tPositive = false;
                    for ( mtk::Vertex* iVertex : tCell.get_vertex_pointers() )
                    {
                        Matrix< DDRMat > tX   = iVertex->get_coords();
                        real             tPhi = iPlane[ 0 ] * tX( 0 ) + iPlane[ 1 ] * tX( 1 ) + iPlane[ 2 ] * tX( 2 ) + iPlane[ 3 ];
                        tNegative             = tNegative or tPhi < 0.0;
                        tPositive             = tPositive or tPhi > 0.0;
                    }
                    if ( tNegative and tPositive )
                    {
                        tLevels.insert( tCell.get_level() );
                    }
                }
            }

            // decompose: true if the level check passed on all ranks
            bool tDecomposed = false;
            try
            {
                Vector< std::shared_ptr< gen::Geometry > > tGeometries;
                for ( uint iPlane = 0; iPlane < iCase.mPlanes.size(); iPlane++ )
                {
                    gPlanes[ iPlane ] = iCase.mPlanes( iPlane );
                    tGeometries.push_back( create_geometry( gPlaneFields[ iPlane ], Regime::STORED ) );
                }
                gen::Geometry_Engine_Parameters tParameters;
                tParameters.mGeometries = tGeometries;
                gen::Geometry_Engine tGeometryEngine( tMesh, tParameters );
                Model                tModel( 3, tMesh, &tGeometryEngine );
                tModel.mVerbose = false;
                tDecomposed     = tModel.decompose( { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 } );
            }
            catch ( const std::exception& aException )
            {
                std::cout << "rank " << par_rank() << " " << iCase.mName << " THREW: " << std::string( aException.what() ).substr( 0, 160 ) << "\n";
            }

            // the levels every rank sees, and whether the whole problem has mixed levels
            uint tMinLevel = min_all( tLevels.empty() ? (uint)99 : *tLevels.begin() );
            uint tMaxLevel = max_all( tLevels.empty() ? (uint)0 : *tLevels.rbegin() );
            bool tMixed    = tMinLevel < tMaxLevel;

            std::ostringstream tLine;
            tLine << "rank " << par_rank() << " of " << par_size() << "  " << std::setw( 26 ) << iCase.mName << "  cut levels on this rank {";
            for ( uint iLevel : tLevels )
            {
                tLine << " " << iLevel;
            }
            tLine << " }  globally " << ( tMixed ? "MIXED" : "uniform" ) << "  decomposed " << tDecomposed << "\n";
            std::cout << tLine.str() << std::flush;

            // the level check must refuse every mixed-level problem, on any number of ranks
            if ( par_rank() == 0 and tMixed == tDecomposed )
            {
                tNumInconsistent++;
                std::cout << "  ^ " << ( tDecomposed ? "decomposed although the cut cells are on different levels" : "refused although every cut cell is on one level" ) << "\n";
            }
        }
        if ( par_rank() == 0 )
        {
            CHECK( tNumInconsistent == 0 );
        }
    }

}    // namespace moris::xtk
