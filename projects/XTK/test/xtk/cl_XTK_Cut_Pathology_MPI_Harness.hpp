/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Cut_Pathology_MPI_Harness.hpp
 *
 * The multi-rank pathology harness on the 4 x 2 x 2 box [0,2] x [0,1] x [0,1]: planes on, near and across the
 * processor boundary (owned regions, gathered integration vertices, the case list), and the aura sweep (owner and
 * aura copies of a background cell decomposed, numbered and enriched identically, output sets independent of the
 * rank count).
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
#include "cl_XTK_Cut_Pathology_Harness.hpp"

namespace moris::xtk
{
    namespace cut_pathology_test
    {
        using namespace pathology_test;

        // ----------------------------------------------------------------------------------
        // multi-rank sweep: planes on, near and across the processor boundary of a 4 x 2 x 2 box [0,2] x [0,1] x [0,1]

        const real gBoxCellSize = 0.5;
        const real gBoxVolume   = 2.0;

        // integration vertices closer than this count as one location
        const real gLocationTolerance = 1.0e-12;

        // exact negative volume of a plane in the box: the 24-tet reference is exact for a linear field, per background cell
        inline real
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

        inline Point
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
        inline void
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
        inline Parallel_Result
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
        inline Vector< Parallel_Case >
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

        // interfaces at the outer edge of the aura rather than on the processor face. With the face at x = 1 on 2 ranks,
        // rank 1's aura ends at x = 0.5 and rank 0's at x = 1.5 both for 4 x 2 x 2 cells with linear B-splines (h = 0.5,
        // one aura layer) and for 8 x 4 x 4 cells with quadratic ones (h = 0.25, two layers). A cut cell in the outermost
        // aura layer lacks its outer neighbours on that rank; its outer vertices have no T-matrix there.
        // aEdges: the aura edges to use (e.g. { 0.5 } on a mesh refined beyond x = 1.5)
        inline Vector< Parallel_Case >
        aura_edge_cases( const Vector< real >& aEdges = { 0.5, 1.5 } )
        {
            const Point tXAxis   = { 1.0, 0.0, 0.0 };
            const Point tSlanted = { 1.0, 0.3, 0.2 };

            Vector< Parallel_Case > tCases;
            for ( real tEdge : aEdges )
            {
                const Point tPoint = { tEdge, 0.5, 0.5 };    // a vertex on the aura edge
                std::ostringstream tPrefix;
                tPrefix << "x=" << tEdge;
                tCases.push_back( { tPrefix.str() + " (aura edge)", { plane_through( tXAxis, tPoint, 0.0 ) } } );
                for ( real tOffset : { 1.0e-11, -1.0e-11, 1.0e-6, -1.0e-6, 0.1, -0.1 } )
                {
                    std::ostringstream tName;
                    tName << tPrefix.str() << ( tOffset > 0.0 ? "+" : "" ) << tOffset;
                    tCases.push_back( { tName.str(), { plane_through( tXAxis, tPoint, tOffset ) } } );
                }
                tCases.push_back( { tPrefix.str() + " slanted", { plane_through( tSlanted, tPoint, 0.0 ) } } );
                tCases.push_back( { tPrefix.str() + " diagonal", { plane_through( { 1.0, 1.0, 0.0 }, tPoint, 0.0 ) } } );
            }
            return tCases;
        }

        // CSV output on rank 0, open only if MORIS_XTK_PATHOLOGY_MPI_CSV is set
        inline std::ofstream*
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

        // ----------------------------------------------------------------------------------
        // aura sweep: a background cell held by several ranks (owner + aura copies) must be decomposed, numbered and
        // enriched identically on each of them, and the enriched meshes and output sets must not depend on the rank count

        // entities held by several ranks: column 0 is the global id, columns aFirst.. must agree between the copies;
        // returns the number of ids whose copies disagree and fills one example
        inline uint
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

        // the converse of count_disagreeing_copies for enriched IP vertices: one vertex, identified by its location and
        // its enriched T-matrix (columns 2-4 and 6-7 of the vertex rows), must carry one id on every rank that holds it.
        // Ids are assigned per unzipped copy before merge_duplicate_interpolation_vertices, and each rank keeps its own
        // first local copy (finding 27 in studies/xtk_holistic_review/findings/15_walkthrough_stops_2_3.md §8.7).
        // Returns the number of vertices with more than one id; aWithinRank counts those where one rank holds two of them
        inline uint
        count_vertices_with_several_ids( const Vector< Matrix< DDRMat > >& aGathered, uint& aWithinRank, uint& aNumKeys, std::string& aExample )
        {
            using Key = std::tuple< std::int64_t, std::int64_t, std::int64_t, real, real >;
            std::map< Key, std::map< moris_id, std::set< uint > > > tIdsPerVertex;
            for ( uint iRank = 0; iRank < aGathered.size(); iRank++ )
            {
                for ( uint iRow = 0; iRow < aGathered( iRank ).n_rows(); iRow++ )
                {
                    const Matrix< DDRMat >& tRows = aGathered( iRank );
                    Key                     tKey  = { std::llround( tRows( iRow, 2 ) / gLocationTolerance ), std::llround( tRows( iRow, 3 ) / gLocationTolerance ), std::llround( tRows( iRow, 4 ) / gLocationTolerance ), tRows( iRow, 6 ), tRows( iRow, 7 ) };
                    tIdsPerVertex[ tKey ][ (moris_id)tRows( iRow, 0 ) ].insert( iRank );
                }
            }
            aNumKeys      = tIdsPerVertex.size();
            aWithinRank   = 0;
            uint tSeveral = 0;
            for ( const std::pair< const Key, std::map< moris_id, std::set< uint > > >& iVertex : tIdsPerVertex )
            {
                if ( iVertex.second.size() < 2 )
                {
                    continue;
                }
                tSeveral++;
                std::map< uint, uint > tIdsOnRank;
                for ( const std::pair< const moris_id, std::set< uint > >& iId : iVertex.second )
                {
                    for ( uint iRank : iId.second )
                    {
                        tIdsOnRank[ iRank ]++;
                    }
                }
                bool tWithin = false;
                for ( const std::pair< const uint, uint >& iRank : tIdsOnRank )
                {
                    tWithin = tWithin or iRank.second > 1;
                }
                aWithinRank += tWithin;
                if ( aExample.empty() )
                {
                    std::ostringstream tStream;
                    tStream << "vertex at (" << std::get< 0 >( iVertex.first ) * gLocationTolerance << ", " << std::get< 1 >( iVertex.first ) * gLocationTolerance
                            << ", " << std::get< 2 >( iVertex.first ) * gLocationTolerance << ") has ids";
                    for ( const std::pair< const moris_id, std::set< uint > >& iId : iVertex.second )
                    {
                        tStream << " " << iId.first << " (rank";
                        for ( uint iRank : iId.second )
                        {
                            tStream << " " << iRank;
                        }
                        tStream << ")";
                    }
                    aExample = tStream.str();
                }
            }
            return tSeveral;
        }

        // corners of IP cells held by several ranks (columns: cell id, x, y, z, T-matrix hash high, low): the enriched
        // T-matrix at a given corner of a given cell must be the same on every rank, since its basis ids are the rows and
        // columns the cell assembles into. Returns the number of (cell, corner) pairs with more than one T-matrix
        inline uint
        count_corners_with_several_t_matrices( const Vector< Matrix< DDRMat > >& aGathered, uint& aNumShared, std::string& aExample )
        {
            using Key = std::tuple< moris_id, std::int64_t, std::int64_t, std::int64_t >;
            std::map< Key, std::map< std::pair< real, real >, std::set< uint > > > tTMatrices;
            std::map< Key, std::set< uint > >                                      tRanks;
            for ( uint iRank = 0; iRank < aGathered.size(); iRank++ )
            {
                const Matrix< DDRMat >& tRows = aGathered( iRank );
                for ( uint iRow = 0; iRow < tRows.n_rows(); iRow++ )
                {
                    Key tKey = { (moris_id)tRows( iRow, 0 ), std::llround( tRows( iRow, 1 ) / gLocationTolerance ), std::llround( tRows( iRow, 2 ) / gLocationTolerance ), std::llround( tRows( iRow, 3 ) / gLocationTolerance ) };
                    tTMatrices[ tKey ][ { tRows( iRow, 4 ), tRows( iRow, 5 ) } ].insert( iRank );
                    tRanks[ tKey ].insert( iRank );
                }
            }
            aNumShared    = 0;
            uint tSeveral = 0;
            for ( const std::pair< const Key, std::map< std::pair< real, real >, std::set< uint > > >& iCorner : tTMatrices )
            {
                aNumShared += tRanks[ iCorner.first ].size() > 1;
                if ( iCorner.second.size() < 2 )
                {
                    continue;
                }
                tSeveral++;
                if ( aExample.empty() )
                {
                    std::ostringstream tStream;
                    tStream << "cell " << std::get< 0 >( iCorner.first ) << " corner (" << std::get< 1 >( iCorner.first ) * gLocationTolerance << ", "
                            << std::get< 2 >( iCorner.first ) * gLocationTolerance << ", " << std::get< 3 >( iCorner.first ) * gLocationTolerance << ") has "
                            << iCorner.second.size() << " different T-matrices on ranks";
                    for ( const std::pair< const std::pair< real, real >, std::set< uint > >& iT : iCorner.second )
                    {
                        tStream << " {";
                        for ( uint iRank : iT.second )
                        {
                            tStream << " " << iRank;
                        }
                        tStream << " }";
                    }
                    aExample = tStream.str();
                }
            }
            return tSeveral;
        }

        // ids that two different entities carry (column 0 = id, column 1 = 1 if owned): owned ids must be unique
        inline uint
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

            // one enriched IP vertex (same location and T-matrix) carrying several ids (finding 27); of these, the ones
            // where a single rank holds two of the ids (an unmerged duplicate) rather than two ranks disagreeing
            uint        mIpVerticesWithSeveralIds = 0, mIpVerticesWithSeveralIdsOnOneRank = 0, mDistinctIpVertices = 0;
            std::string mIpVertexIdExample;

            // a corner of an IP cell (keyed by cell id = subphase id and location, so independent of vertex ids) whose
            // enriched T-matrix (basis ids, owners, weights) differs between the ranks that hold the cell: finding 28
            uint        mSharedCorners = 0, mCornerTMatrixMismatch = 0;
            std::string mCornerTMatrixExample;

            // background cells this rank holds (owned + aura); on rank 0 after aura_case: the min and max over ranks,
            // which give the aura depth (e.g. 8 x 4 x 4 on 2 ranks: 64 owned + 16 per aura layer)
            uint mNumBgCellsOnRank = 0, mMinBgCellsOnRank = 0, mMaxBgCellsOnRank = 0;

            // informational: aura copies of an IP cell whose vertex ids differ from the owner's (the aura copy keeps the
            // owner's pre-merge copy ids, which the owner's merge removed: finding 27), and vertices without a T-matrix
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

        inline Aura_Result
        aura_case_fields( mtk::Interpolation_Mesh* aMesh, const Vector< gen::Field_Function >& aFields, Regime aRegime )
        {
            Aura_Result tResult;
            uint        tRank = par_rank();

            Vector< Vector< real > >                         tBgRows, tIgCellRows, tSubphaseRows, tVertexRows, tCellRows, tBasisRows, tCornerRows;
            std::map< moris_index, std::pair< real, real > > tTMatrixHashOfVertex;    // enriched IP vertex index -> T-matrix hash
            std::map< std::string, std::pair< uint, real > > tLocalSets;
            bool                                             tThrew = false;

            // a rank that throws inside a collective call leaves the others waiting; the sweep then hangs (timeout)
            try
            {
                Vector< std::shared_ptr< gen::Geometry > > tGeometries;
                for ( const gen::Field_Function& iField : aFields )
                {
                    tGeometries.push_back( create_geometry( iField, aRegime ) );
                }
                gen::Geometry_Engine_Parameters tParameters;
                tParameters.mGeometries = tGeometries;
                gen::Geometry_Engine tGeometryEngine( aMesh, tParameters );

                Model tModel( 3, aMesh, &tGeometryEngine );
                tModel.mVerbose = false;
                if ( not tModel.decompose( { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 } ) )
                {
                    throw std::runtime_error( "decomposition refused (e.g. intersected cells on different refinement levels)" );
                }
                tResult.mNumBgCellsOnRank = aMesh->get_num_elems();
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
                    tTMatrixHashOfVertex[ tVertex.get_index() ] = { tBasisHash.high(), tBasisHash.low() };

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

                    // vertex locations (must agree) and vertex ids (reported only: aura copies keep pre-merge ids, finding 27)
                    Vector< std::array< std::int64_t, 3 > > tLocations;
                    Vector< moris_id >                      tIds;
                    for ( mtk::Vertex* iVertex : tCell.get_vertex_pointers() )
                    {
                        Matrix< DDRMat > tX = iVertex->get_coords();
                        tLocations.push_back( { std::llround( tX( 0 ) / gLocationTolerance ), std::llround( tX( 1 ) / gLocationTolerance ),
                                std::llround( tX( 2 ) / gLocationTolerance ) } );
                        tIds.push_back( iVertex->get_id() );

                        // the enriched T-matrix at this corner of this IP cell (cell id = subphase id, rank independent)
                        std::map< moris_index, std::pair< real, real > >::const_iterator tHash = tTMatrixHashOfVertex.find( iVertex->get_index() );
                        if ( tHash != tTMatrixHashOfVertex.end() )
                        {
                            tCornerRows.push_back( { (real)tCell.get_id(), tX( 0 ), tX( 1 ), tX( 2 ), tHash->second.first, tHash->second.second } );
                        }
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
            tResult.mMinBgCellsOnRank            = min_all( tResult.mNumBgCellsOnRank );
            tResult.mMaxBgCellsOnRank            = max_all( tResult.mNumBgCellsOnRank );

            Vector< Matrix< DDRMat > > tBg       = gather_rows( tBgRows, 10 );
            Vector< Matrix< DDRMat > > tIgCells  = gather_rows( tIgCellRows, 2 );
            Vector< Matrix< DDRMat > > tSps      = gather_rows( tSubphaseRows, 2 );
            Vector< Matrix< DDRMat > > tVertices = gather_rows( tVertexRows, 9 );
            Vector< Matrix< DDRMat > > tCells    = gather_rows( tCellRows, 7 );
            Vector< Matrix< DDRMat > > tCorners  = gather_rows( tCornerRows, 6 );
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
                tResult.mIpVerticesWithSeveralIds =
                        count_vertices_with_several_ids( tVertices, tResult.mIpVerticesWithSeveralIdsOnOneRank, tResult.mDistinctIpVertices, tResult.mIpVertexIdExample );
                tResult.mIpCellMismatch   = count_disagreeing_copies( tCells, 1, 3, tResult.mIpCellExample, tResult.mSharedIpCells );
                uint tSharedForIds        = 0;
                tResult.mIpCellIdMismatch = count_disagreeing_copies( tCells, 4, 5, tResult.mIpCellIdExample, tSharedForIds );
                tResult.mCornerTMatrixMismatch =
                        count_corners_with_several_t_matrices( tCorners, tResult.mSharedCorners, tResult.mCornerTMatrixExample );

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

        // one geometry per plane (at most 6, the slots of gPlaneFields)
        inline Aura_Result
        aura_case( mtk::Interpolation_Mesh* aMesh, const Vector< Plane >& aPlanes, Regime aRegime )
        {
            Vector< gen::Field_Function > tFields;
            for ( uint iPlane = 0; iPlane < aPlanes.size(); iPlane++ )
            {
                gPlanes[ iPlane ] = aPlanes( iPlane );
                tFields.push_back( gPlaneFields[ iPlane ] );
            }
            return aura_case_fields( aMesh, tFields, aRegime );
        }

        // ----------------------------------------------------------------------------------
        // thin features (finding 28): a union of axis-aligned boxes, positive inside. On the [0,2] x [0,1] x [0,1] box with
        // 8 x 4 x 4 cells (h = 0.25) a box of half-thickness 0.06 centred on a node row holds only that row of nodes, so it
        // is one node layer thick. A "U" of two such slabs joined by a bridge two cells beyond the processor face (x = 1)
        // makes a quadratic basis next to the face see one connected piece on the rank that holds the bridge cell and two
        // pieces on the rank whose aura stops short of it.

        struct Box
        {
            Point mCentre;
            Point mHalfWidth;
        };

        inline Vector< Box > gBoxes;

        inline real
        boxes_field( const Matrix< DDRMat >& aPoint, const Vector< real >& aParameters )
        {
            real tValue = -MORIS_REAL_MAX;
            for ( const Box& iBox : gBoxes )
            {
                real tInside = MORIS_REAL_MAX;
                for ( uint iDim = 0; iDim < 3; iDim++ )
                {
                    tInside = std::min( tInside, iBox.mHalfWidth[ iDim ] - std::abs( aPoint( iDim ) - iBox.mCentre[ iDim ] ) );
                }
                tValue = std::max( tValue, tInside );
            }
            return tValue;
        }

        struct Thin_Case
        {
            std::string   mName;
            Vector< Box > mBoxes;
        };

        inline Vector< Thin_Case >
        thin_cases()
        {
            const real tThin = 0.06;    // half-thickness: one node layer
            const real tFull = 0.6;     // half-width through z (covers [0,1])

            // a slab along x at node row y = aY, from aX0 to aX1
            auto tSlab = [ & ]( real aY, real aX0, real aX1 ) -> Box {
                return { { 0.5 * ( aX0 + aX1 ), aY, 0.5 }, { 0.5 * ( aX1 - aX0 ), tThin, tFull } };
            };
            // a bridge along y at node column x = aX, joining rows y = 0.25 and y = 0.75
            auto tBridge = [ & ]( real aX ) -> Box {
                return { { aX, 0.5, 0.5 }, { tThin, 0.25 + tThin, tFull } };
            };

            Vector< Thin_Case > tCases;
            tCases.push_back( { "thin slab across the face", { tSlab( 0.25, 0.3, 1.7 ) } } );
            tCases.push_back( { "U, bridge at x=1.5", { tSlab( 0.25, 0.4, 1.5 + tThin ), tSlab( 0.75, 0.4, 1.5 + tThin ), tBridge( 1.5 ) } } );
            tCases.push_back( { "U, bridge at x=0.5", { tSlab( 0.25, 0.5 - tThin, 1.6 ), tSlab( 0.75, 0.5 - tThin, 1.6 ), tBridge( 0.5 ) } } );
            return tCases;
        }

        inline std::ofstream*
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
    }    // namespace cut_pathology_test
}    // namespace moris::xtk
