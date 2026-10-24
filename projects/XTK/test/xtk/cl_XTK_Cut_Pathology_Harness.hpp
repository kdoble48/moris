/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Cut_Pathology_Harness.hpp
 *
 * The single-hex pathology harness: one decomposition with its invariant checks (Case_Result), the sweep
 * summaries and their printing, the corner-value and two-geometry sweeps, and the multi-plane sweeps.
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
#include "cl_XTK_Cut_Pathology_Reference.hpp"

namespace moris::xtk
{
    namespace cut_pathology_test
    {
        using namespace pathology_test;

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

        inline bool
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
        inline Moments
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
        inline real
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

        inline Case_Result
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

        inline std::string
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
        inline std::ofstream*
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

        inline void
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
        inline bool
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

        inline void
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
        inline Vector< Corner_Values >
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
        inline Sweep_Summary
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
        inline Vector< std::pair< std::string, Corner_Values > >
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
        inline Sweep_Summary
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

        inline void
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
        inline Plane
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

        inline Corner_Values
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

        inline Case_Result
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

        inline void
        print_metrics_header( const std::string& aTitle, const std::string& aExtraColumns = "" )
        {
            std::cout << "\n=== " << aTitle << " ===\n"
                      << std::setw( 26 ) << "case" << std::setw( 7 ) << "cells" << std::setw( 6 ) << "subph"
                      << std::setw( 13 ) << "min cell V" << std::setw( 6 ) << "V<=0" << std::setw( 13 ) << "min subph V"
                      << std::setw( 7 ) << "ifacet" << std::setw( 13 ) << "min ifacet A" << std::setw( 13 ) << "total V - 1"
                      << aExtraColumns << "\n";
        }

        inline void
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
    }    // namespace cut_pathology_test
}    // namespace moris::xtk
