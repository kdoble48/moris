/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_XTK_Cut_Mesh_Global_Ids.cpp
 *
 * The global ids that the XTK integration meshes report through get_glb_entity_id_from_entity_loc_index() (the ids the
 * exodus writer puts into the node and element maps of xtk_temp.exo) must be the parallel-global ids: on every rank the
 * query returns get_id() and inverts get_loc_entity_ind_from_entity_glb_id(), and across ranks one node id is one
 * location, one location is one node id, and one cell id is one cell (same vertices). A sphere cuts a 8 x 4 x 4 box
 * [0,2] x [0,1] x [0,1] across the processor boundaries; the mesh is decomposed, enriched and ghost-stabilized (the
 * ghost cells get their ids from Enriched_Integration_Mesh::allocate_entity_ids(), seeded by get_max_entity_id()).
 * Runs on any number of ranks: XTK-test.exe "[XTK_CUT_MESH_GLOBAL_IDS]".
 *
 */

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <tuple>
#include <utility>

#include "catch.hpp"

#include "cl_Communication_Tools.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Cut_Pathology_Reference.hpp"
#include "cl_XTK_Pathology_Test_Tools.hpp"

namespace moris::xtk
{
    namespace global_id_test
    {
        using namespace pathology_test;
        using namespace cut_pathology_test;

        // sphere across the processor boundaries of the box
        inline real
        sphere( const Matrix< DDRMat >& aPoint, const Vector< real >& aParameters )
        {
            real tDx = aPoint( 0 ) - 0.93;
            real tDy = aPoint( 1 ) - 0.52;
            real tDz = aPoint( 2 ) - 0.47;
            return std::sqrt( tDx * tDx + tDy * tDy + tDz * tDz ) - 0.36;
        }

        // per-rank checks of one mesh, and the rows (query id, ...) to compare across ranks
        struct Id_Rows
        {
            uint                     mNumQueryNotId  = 0;    // query result != get_id()
            uint                     mNumNoRoundTrip = 0;    // id -> index -> id does not come back
            moris_id                 mMaxNodeId      = 0;    // largest queried node id
            moris_id                 mMaxCellId      = 0;    // largest queried cell id
            Vector< Vector< real > > mVertexRows;            // { id, x, y, z }
            Vector< Vector< real > > mCellRows;              // { id, number of vertices, hash of the sorted vertex ids (high, low) }
        };

        // id -> index lookup that reports a miss instead of throwing (a throw on one rank only would hang the collectives)
        inline bool
        round_trips( mtk::Mesh& aMesh, moris_id aId, mtk::EntityRank aRank, moris_index aIndex )
        {
            try
            {
                return aMesh.get_loc_entity_ind_from_entity_glb_id( aId, aRank ) == aIndex;
            } catch ( const std::exception& )
            {
                return false;
            }
        }

        inline Id_Rows
        collect_id_rows( mtk::Mesh& aMesh )
        {
            Id_Rows tRows;

            uint tNumVertices = aMesh.get_num_entities( mtk::EntityRank::NODE );
            for ( uint iVertex = 0; iVertex < tNumVertices; iVertex++ )
            {
                mtk::Vertex const & tVertex = aMesh.get_mtk_vertex( (moris_index)iVertex );
                moris_id            tId     = aMesh.get_glb_entity_id_from_entity_loc_index( (moris_index)iVertex, mtk::EntityRank::NODE );

                tRows.mNumQueryNotId += tId != tVertex.get_id();
                tRows.mNumNoRoundTrip += !round_trips( aMesh, tId, mtk::EntityRank::NODE, (moris_index)iVertex );
                tRows.mMaxNodeId = std::max( tRows.mMaxNodeId, tId );

                Matrix< DDRMat > tCoords = tVertex.get_coords();
                tRows.mVertexRows.push_back( { (real)tId, tCoords( 0 ), tCoords( 1 ), tCoords( 2 ) } );
            }

            uint tNumCells = aMesh.get_num_entities( mtk::EntityRank::ELEMENT );
            for ( uint iCell = 0; iCell < tNumCells; iCell++ )
            {
                mtk::Cell const & tCell = aMesh.get_mtk_cell( (moris_index)iCell );
                moris_id          tId   = aMesh.get_glb_entity_id_from_entity_loc_index( (moris_index)iCell, mtk::EntityRank::ELEMENT );

                tRows.mNumQueryNotId += tId != tCell.get_id();
                tRows.mNumNoRoundTrip += !round_trips( aMesh, tId, mtk::EntityRank::ELEMENT, (moris_index)iCell );
                tRows.mMaxCellId = std::max( tRows.mMaxCellId, tId );

                // the vertices by their queried ids, sorted (aura copies may order them differently)
                Vector< mtk::Vertex* > tVertices = tCell.get_vertex_pointers();
                Vector< moris_id >     tVertexIds;
                for ( mtk::Vertex* iVertex : tVertices )
                {
                    tVertexIds.push_back( aMesh.get_glb_entity_id_from_entity_loc_index( iVertex->get_index(), mtk::EntityRank::NODE ) );
                }
                std::sort( tVertexIds.begin(), tVertexIds.end() );
                connectivity_check::Hash tHash;
                for ( moris_id iVertexId : tVertexIds )
                {
                    tHash.add( iVertexId );
                }
                tRows.mCellRows.push_back( { (real)tId, (real)tVertexIds.size(), tHash.high(), tHash.low() } );
            }
            return tRows;
        }

        // cross-rank agreement of the gathered rows (rank 0 only)
        struct Id_Agreement
        {
            uint mNumVertexRows        = 0;
            uint mNumVertexIds         = 0;
            uint mNumIdsAtTwoLocations = 0;    // one node id, several locations
            uint mNumLocationsWithIds  = 0;    // one location, different node ids on different ranks
            uint mNumLocalDuplicates   = 0;    // locations with several vertices on one rank (excluded above; not an id defect)
            uint mNumCellIds           = 0;
            uint mNumCellIdsTwoCells   = 0;    // one cell id, different vertices
        };

        inline Id_Agreement
        check_agreement( const Id_Rows& aRows )
        {
            Vector< Matrix< DDRMat > > tVertices = gather_rows( aRows.mVertexRows, 4 );
            Vector< Matrix< DDRMat > > tCells    = gather_rows( aRows.mCellRows, 4 );

            Id_Agreement tAgreement;
            if ( par_rank() != 0 )
            {
                return tAgreement;
            }

            // node id -> first location; rounded location -> first id
            std::map< moris_id, std::tuple< real, real, real > > tLocationOfId;
            std::map< std::tuple< long, long, long >, moris_id > tIdAtLocation;
            std::map< moris_id, bool >                           tIdAtTwoLocations;
            uint                                                 tNumLocationsWithIds = 0;
            for ( Matrix< DDRMat > const & iRank : tVertices )
            {
                // this rank's locations (rounded to 1e-9: cut points are not on a lattice, so a tolerance key is enough here);
                // a location with several vertices on one rank (e.g. a decomposition vertex on a hanging face of a refined
                // background mesh) is a property of the decomposition, not of the ids, and is left out of the cross-rank check
                std::map< std::tuple< long, long, long >, moris_id > tRankIdAtLocation;
                std::map< std::tuple< long, long, long >, bool >     tRankDuplicate;
                for ( uint iRow = 0; iRow < iRank.n_rows(); iRow++ )
                {
                    tAgreement.mNumVertexRows++;
                    moris_id                       tId       = (moris_id)iRank( iRow, 0 );
                    std::tuple< real, real, real > tLocation = { iRank( iRow, 1 ), iRank( iRow, 2 ), iRank( iRow, 3 ) };

                    std::map< moris_id, std::tuple< real, real, real > >::iterator tSeen = tLocationOfId.find( tId );
                    if ( tSeen == tLocationOfId.end() )
                    {
                        tLocationOfId[ tId ] = tLocation;
                    }
                    else
                    {
                        real tDistance = std::abs( std::get< 0 >( tSeen->second ) - std::get< 0 >( tLocation ) )
                                       + std::abs( std::get< 1 >( tSeen->second ) - std::get< 1 >( tLocation ) )
                                       + std::abs( std::get< 2 >( tSeen->second ) - std::get< 2 >( tLocation ) );
                        if ( tDistance > 1.0e-12 )
                        {
                            tIdAtTwoLocations[ tId ] = true;
                        }
                    }

                    std::tuple< long, long, long > tKey = { std::lround( std::get< 0 >( tLocation ) * 1.0e9 ),
                        std::lround( std::get< 1 >( tLocation ) * 1.0e9 ),
                        std::lround( std::get< 2 >( tLocation ) * 1.0e9 ) };
                    if ( tRankIdAtLocation.find( tKey ) == tRankIdAtLocation.end() )
                    {
                        tRankIdAtLocation[ tKey ] = tId;
                    }
                    else
                    {
                        tRankDuplicate[ tKey ] = true;
                    }
                }
                tAgreement.mNumLocalDuplicates += tRankDuplicate.size();

                for ( std::pair< const std::tuple< long, long, long >, moris_id > const & iLocation : tRankIdAtLocation )
                {
                    if ( tRankDuplicate.find( iLocation.first ) != tRankDuplicate.end() )
                    {
                        continue;
                    }
                    std::map< std::tuple< long, long, long >, moris_id >::iterator tAt = tIdAtLocation.find( iLocation.first );
                    if ( tAt == tIdAtLocation.end() )
                    {
                        tIdAtLocation[ iLocation.first ] = iLocation.second;
                    }
                    else if ( tAt->second != iLocation.second )
                    {
                        tNumLocationsWithIds++;
                    }
                }
            }
            tAgreement.mNumVertexIds         = tLocationOfId.size();
            tAgreement.mNumIdsAtTwoLocations = tIdAtTwoLocations.size();
            tAgreement.mNumLocationsWithIds  = tNumLocationsWithIds;

            std::map< moris_id, std::tuple< real, real, real > > tSignatureOfCell;
            std::map< moris_id, bool >                           tCellIdTwoCells;
            for ( Matrix< DDRMat > const & iRank : tCells )
            {
                for ( uint iRow = 0; iRow < iRank.n_rows(); iRow++ )
                {
                    moris_id                       tId        = (moris_id)iRank( iRow, 0 );
                    std::tuple< real, real, real > tSignature = { iRank( iRow, 1 ), iRank( iRow, 2 ), iRank( iRow, 3 ) };

                    std::map< moris_id, std::tuple< real, real, real > >::iterator tSeen = tSignatureOfCell.find( tId );
                    if ( tSeen == tSignatureOfCell.end() )
                    {
                        tSignatureOfCell[ tId ] = tSignature;
                    }
                    else if ( tSeen->second != tSignature )
                    {
                        tCellIdTwoCells[ tId ] = true;
                    }
                }
            }
            tAgreement.mNumCellIds         = tSignatureOfCell.size();
            tAgreement.mNumCellIdsTwoCells = tCellIdTwoCells.size();
            return tAgreement;
        }

        inline void
        require_global_ids( mtk::Mesh& aMesh, const std::string& aName, bool aHasMaxEntityId )
        {
            Id_Rows      tRows      = collect_id_rows( aMesh );
            Id_Agreement tAgreement = check_agreement( tRows );

            // get_max_entity_id() seeds new ids (ghost cells): it must bound every id (only the enriched mesh implements it)
            moris_id tMaxQueriedNodeId = max_all( tRows.mMaxNodeId );
            moris_id tMaxQueriedCellId = max_all( tRows.mMaxCellId );
            moris_id tMaxNodeId        = aHasMaxEntityId ? aMesh.get_max_entity_id( mtk::EntityRank::NODE ) : tMaxQueriedNodeId;
            moris_id tMaxCellId        = aHasMaxEntityId ? aMesh.get_max_entity_id( mtk::EntityRank::ELEMENT ) : tMaxQueriedCellId;

            uint tNumQueryNotId  = sum_all( tRows.mNumQueryNotId );
            uint tNumNoRoundTrip = sum_all( tRows.mNumNoRoundTrip );

            if ( par_rank() == 0 )
            {
                std::cout << aName << " on " << par_size() << " rank(s): " << tAgreement.mNumVertexRows << " vertex rows on " << tAgreement.mNumVertexIds
                          << " node ids, " << tAgreement.mNumIdsAtTwoLocations << " ids at two locations, " << tAgreement.mNumLocationsWithIds
                          << " locations with two ids (" << tAgreement.mNumLocalDuplicates << " rank-local duplicate locations); " << tAgreement.mNumCellIds << " cell ids, " << tAgreement.mNumCellIdsTwoCells
                          << " on two cells; query != get_id(): " << tNumQueryNotId << ", no round trip: " << tNumNoRoundTrip
                          << "; max node / cell id " << tMaxNodeId << " / " << tMaxCellId << " vs queried " << tMaxQueriedNodeId << " / "
                          << tMaxQueriedCellId << "\n"
                          << std::flush;

                CHECK( tNumQueryNotId == 0 );
                CHECK( tNumNoRoundTrip == 0 );
                CHECK( tAgreement.mNumIdsAtTwoLocations == 0 );
                CHECK( tAgreement.mNumLocationsWithIds == 0 );
                CHECK( tAgreement.mNumCellIdsTwoCells == 0 );
                CHECK( tMaxNodeId >= tMaxQueriedNodeId );
                CHECK( tMaxCellId >= tMaxQueriedCellId );
            }
        }
        // MB-107: ids of the enriched interpolation mesh after ghost setup; one id is one cell across ranks, and new ids
        // (allocate_entity_ids( ..., aStartFresh = false ), as the ghost setup requests them) start after the largest id
        inline void
        require_interpolation_cell_ids( Enriched_Interpolation_Mesh& aMesh, const std::string& aName )
        {
            Id_Rows      tRows      = collect_id_rows( aMesh );
            Id_Agreement tAgreement = check_agreement( tRows );

            moris_id tMaxCellId        = aMesh.get_max_entity_id( mtk::EntityRank::ELEMENT );
            moris_id tMaxQueriedCellId = max_all( tRows.mMaxCellId );

            // every rank asks for a few ids; none may collide with an existing cell id
            moris_id tFirstNewId    = aMesh.allocate_entity_ids( 3, mtk::EntityRank::ELEMENT, false );
            uint     tNumCollisions = sum_all( (uint)( tFirstNewId <= tMaxCellId ) );

            if ( par_rank() == 0 )
            {
                std::cout << aName << " on " << par_size() << " rank(s): " << tAgreement.mNumCellIds << " cell ids, " << tAgreement.mNumCellIdsTwoCells
                          << " on two cells; max cell id " << tMaxCellId << " vs queried " << tMaxQueriedCellId << "; first new id on rank 0 "
                          << tFirstNewId << ", ranks whose new ids collide: " << tNumCollisions << "\n"
                          << std::flush;

                CHECK( tAgreement.mNumCellIdsTwoCells == 0 );
                CHECK( tMaxCellId >= tMaxQueriedCellId );
                CHECK( tNumCollisions == 0 );
            }
        }

        // decompose, check the cut mesh, enrich, build the ghost cells, check the enriched mesh
        inline void
        check_box( real aRefineAboveX, const std::string& aName )
        {
            Grid_Mesh                tBox( { 8, 4, 4 }, { 2.0, 1.0, 1.0 }, aRefineAboveX );
            mtk::Interpolation_Mesh* tMesh = tBox.mMesh;

            gen::Geometry_Engine_Parameters tParameters;
            tParameters.mGeometries = { create_geometry( &sphere, Regime::STORED ) };
            gen::Geometry_Engine tGeometryEngine( tMesh, tParameters );

            Model tModel( 3, tMesh, &tGeometryEngine );
            tModel.mVerbose = false;
            REQUIRE( tModel.decompose( { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 } ) );

            // the cut mesh as decomposed
            require_global_ids( *tModel.get_cut_integration_mesh(), aName + ": cut integration mesh", false );

            // the enriched integration mesh (what xtk_temp.exo and FEM see) forwards the query to the cut mesh; the ghost
            // cells at refinement transitions take new ids from get_max_entity_id() + 1
            tModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE, 0 );
            tModel.construct_face_oriented_ghost_penalization_cells();
            require_global_ids( tModel.get_enriched_integ_mesh( 0 ), aName + ": enriched integration mesh with ghost", true );

            // the enriched interpolation mesh after the ghost setup (MB-107)
            require_interpolation_cell_ids( tModel.get_enriched_interp_mesh( 0 ), aName + ": enriched interpolation mesh with ghost" );
        }
    }    // namespace global_id_test

    TEST_CASE( "XTK integration mesh global ids are the parallel ids", "[XTK][XTK_CUT_MESH_GLOBAL_IDS]" )
    {
        using namespace global_id_test;

        // uniform background mesh
        check_box( MORIS_REAL_MAX, "uniform" );

        // background cells at x > 1.4 refined once: the sphere (x up to 1.29) cuts only level-0 cells, and the ghost facets
        // between the cut cells at x in [1.25, 1.5] and their refined neighbours need transition cells with new ids
        check_box( 1.4, "refined" );
    }
}    // namespace moris::xtk
