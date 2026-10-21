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
 * The sweeps are hidden ([.]) because they characterize current behaviour rather than gate it. Run one with
 *   XTK-test.exe "[XTK_PATHOLOGY_HEX_SIGNS]"
 * Set MORIS_XTK_PATHOLOGY_CSV=<file> to write one row per decomposed case.
 *
 */

#include "catch.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR_Parameters.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_GEN_Level_Set_Geometry.hpp"
#include "cl_GEN_User_Defined_Field.hpp"

namespace moris::xtk
{
    namespace
    {
        using Corner_Values = std::array< real, 8 >;
        using Point         = std::array< real, 3 >;

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

        enum class Regime
        {
            STORED,
            ANALYTIC,
            MULTILINEAR
        };

        std::string
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

        // volume of { f < 0 } in a tet with a linear f given by its vertex values
        real
        tet_negative_volume( const std::array< Point, 4 >& aP, const std::array< real, 4 >& aF )
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

            real tVolume = tet_volume( aP[ 0 ], aP[ 1 ], aP[ 2 ], aP[ 3 ] );
            if ( tNeg.size() == 0 )
            {
                return 0.0;
            }
            if ( not tHasPositive )
            {
                return tVolume;
            }
            if ( tNeg.size() == 1 )
            {
                uint tU = tNeg( 0 );
                return tet_volume( aP[ tU ],
                        zero_on_edge( aP[ tU ], aF[ tU ], aP[ tNonNeg( 0 ) ], aF[ tNonNeg( 0 ) ] ),
                        zero_on_edge( aP[ tU ], aF[ tU ], aP[ tNonNeg( 1 ) ], aF[ tNonNeg( 1 ) ] ),
                        zero_on_edge( aP[ tU ], aF[ tU ], aP[ tNonNeg( 2 ) ], aF[ tNonNeg( 2 ) ] ) );
            }
            if ( tNeg.size() == 3 )
            {
                // the single non-negative vertex is positive (otherwise handled above): cut off its corner tet
                uint tW = tNonNeg( 0 );
                return tVolume
                     - tet_volume( aP[ tW ],
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
            return tet_volume( aP[ tA ], tAC, tAD, aP[ tB ] ) + tet_volume( tAC, tAD, aP[ tB ], tBC ) + tet_volume( tAD, aP[ tB ], tBC, tBD );
        }

        // negative volume of the piecewise-linear interpolant on the 24-tet subdivision of [0,1]^3
        real
        reference_negative_volume( const Corner_Values& aCorners )
        {
            Point tCellCentre = { 0.5, 0.5, 0.5 };
            real  tCellValue  = 0.0;
            for ( real iValue : aCorners )
            {
                tCellValue += iValue / 8.0;
            }

            real tVolume = 0.0;
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
                        tVolume += tet_negative_volume(
                                { corner_point( tFirst ), corner_point( tSecond ), tFaceCentre, tCellCentre },
                                { aCorners[ tFirst ], aCorners[ tSecond ], tFaceValue, tCellValue } );
                    }
                }
            }
            return tVolume;
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
            for ( real tDelta : { 1.0e-13, 1.0e-11 } )
            {
                Sweep_Summary tSummary = single_geometry_sweep( "near_zero", { -1.0, -tDelta, tDelta, 1.0 }, Regime::STORED );
                std::ostringstream tTitle;
                tTitle << "single hex, corner values {-1,-d,+d,1}, d = " << tDelta << ", stored";
                print_summary( tTitle.str(), tSummary );
                check_invariants( tSummary, true, false );
            }
        }
    }

    TEST_CASE( "XTK pathology sweep: single hex, two geometries in both orders", "[.][XTK][XTK_PATHOLOGY][XTK_PATHOLOGY_HEX_TWO_GEOMETRIES]" )
    {
        if ( par_size() == 1 )
        {
            for ( Regime iRegime : { Regime::STORED, Regime::ANALYTIC } )
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

}    // namespace moris::xtk
