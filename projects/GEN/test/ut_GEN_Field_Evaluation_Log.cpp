/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_GEN_Field_Evaluation_Log.cpp
 *
 * Diagnostic log of analytic level-set evaluations (log_field_evaluations, log_field_evaluations_box): a circle on a
 * 4 x 4 mesh, cut by the regular-subdivision edges of every cell (corner edges and corner-to-centre spokes) with
 * use_analytic_intersection on.
 *
 */

#include "catch.hpp"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "cl_GEN_Geometry_Engine_Test.hpp"
#include "cl_GEN_Design_Factory.hpp"
#include "cl_GEN_Field_Evaluation_Log.hpp"
#include "cl_GEN_Level_Set_Geometry.hpp"
#include "cl_GEN_Derived_Node.hpp"
#include "fn_GEN_create_simple_mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "fn_PRM_GEN_Parameters.hpp"

namespace moris::gen
{
    namespace
    {
        struct Log_Row
        {
            real        mX, mY, mValue;
            std::string mPurpose, mStage, mOperation;
            sint        mNodeIndex;
        };

        std::vector< Log_Row >
        read_log( const std::string& aFileName )
        {
            std::vector< Log_Row > tRows;
            std::ifstream          tFile( aFileName );
            std::string            tLine;
            std::getline( tFile, tLine );    // header
            while ( std::getline( tFile, tLine ) )
            {
                std::vector< std::string > tItems;
                std::stringstream          tStream( tLine );
                std::string                tItem;
                while ( std::getline( tStream, tItem, ',' ) )
                {
                    tItems.push_back( tItem );
                }
                REQUIRE( tItems.size() == 11 );
                tRows.push_back( { std::stod( tItems[ 0 ] ), std::stod( tItems[ 1 ] ), std::stod( tItems[ 3 ] ), tItems[ 4 ], tItems[ 5 ], tItems[ 6 ], std::stoi( tItems[ 7 ] ) } );
            }
            return tRows;
        }

        struct Cut_Edge
        {
            Matrix< DDRMat > mFirst, mSecond;
        };

        /**
         * Background-node regions (one evaluation per node, stage "ut_background_nodes"), then the regular-subdivision
         * edges of every cell queued for intersection. Returns the local coordinates of the intersections and the
         * cut edges (sign change between the end values).
         */
        Vector< real >
        cut_circle( const std::string& aLogFile, const std::string& aBox, std::vector< Cut_Edge >& aCutEdges, uint& aNumNodesInBox )
        {
            mtk::Interpolation_Mesh* tMesh = create_simple_mesh( 4, 4 );

            Submodule_Parameter_Lists tFieldParameterLists( "FIELDS" );
            tFieldParameterLists.add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
            tFieldParameterLists.set( "center_x", 0.1 );
            tFieldParameterLists.set( "center_y", -0.05 );
            tFieldParameterLists.set( "radius", 0.63 );
            tFieldParameterLists.set( "discretization_mesh_index", -2 );
            tFieldParameterLists.set( "use_analytic_intersection", true );
            if ( not aLogFile.empty() )
            {
                tFieldParameterLists.set( "log_field_evaluations", aLogFile );
                tFieldParameterLists.set( "log_field_evaluations_box", aBox );
            }

            Vector< real > tLocalCoordinates;
            {
                Geometry_Engine_Parameters tGeometryEngineParameters;
                ADV_Manager                tADVManager;
                Design_Factory             tDesignFactory( tFieldParameterLists, tADVManager );
                tGeometryEngineParameters.mGeometries = tDesignFactory.get_geometries();
                Geometry_Engine_Test tGeometryEngine( tMesh, tGeometryEngineParameters );
                Node_Manager&        tNodeManager = tGeometryEngine.get_node_manager();

                // box of the background-node count (as the log box; the whole plane if none)
                std::vector< real > tBox = { -1e300, 1e300, -1e300, 1e300 };
                if ( not aBox.empty() )
                {
                    std::stringstream tStream( aBox );
                    std::string       tItem;
                    for ( uint iItem = 0; iItem < 4 and std::getline( tStream, tItem, ',' ); iItem++ )
                    {
                        tBox[ iItem ] = std::stod( tItem );
                    }
                }

                // each background node once
                aNumNodesInBox = 0;
                {
                    Evaluation_Stage tStage( "ut_background_nodes" );
                    for ( uint iNode = 0; iNode < tMesh->get_num_nodes(); iNode++ )
                    {
                        Matrix< DDRMat > tX = tMesh->get_node_coordinate( iNode );
                        tGeometryEngine.get_geometric_region( 0, iNode, tX );
                        aNumNodesInBox += tX( 0 ) >= tBox[ 0 ] and tX( 0 ) <= tBox[ 1 ] and tX( 1 ) >= tBox[ 2 ] and tX( 1 ) <= tBox[ 3 ];
                    }
                }

                // regular-subdivision edges of every cell: 4 cell edges and 4 spokes to the centre (a derived node)
                Evaluation_Stage tStage( "ut_cut" );
                for ( uint iCell = 0; iCell < tMesh->get_num_elems(); iCell++ )
                {
                    Matrix< IndexMat > tSigned = tMesh->get_nodes_connected_to_element_loc_inds( iCell );
                    Matrix< DDUMat >   tNodes( 4, 1 );
                    Matrix< DDRMat >   tCentre( 1, 2, 0.0 );
                    for ( uint iNode = 0; iNode < 4; iNode++ )
                    {
                        tNodes( iNode ) = tSigned( iNode );
                        tCentre         = tCentre + 0.25 * tMesh->get_node_coordinate( tNodes( iNode ) );
                    }
                    Vector< Matrix< DDRMat > > tXi( 4 );
                    for ( uint iNode = 0; iNode < 4; iNode++ )
                    {
                        Matrix< DDRMat > tX = tMesh->get_node_coordinate( tNodes( iNode ) );
                        tXi( iNode )        = { { tX( 0 ) > tCentre( 0 ) ? 1.0 : -1.0, tX( 1 ) > tCentre( 1 ) ? 1.0 : -1.0 } };
                    }
                    uint tCentreIndex = tNodeManager.get_total_number_of_nodes();
                    tGeometryEngine.create_new_derived_nodes( { tSigned }, { Matrix< DDRMat >( 1, 2, 0.0 ) }, mtk::Geometry_Type::QUAD, mtk::Interpolation_Order::LINEAR );

                    std::vector< std::pair< uint, uint > > tEdges;
                    for ( uint iFirst = 0; iFirst < 4; iFirst++ )
                    {
                        for ( uint iSecond = iFirst + 1; iSecond < 4; iSecond++ )
                        {
                            if ( ( tXi( iFirst )( 0 ) != tXi( iSecond )( 0 ) ) != ( tXi( iFirst )( 1 ) != tXi( iSecond )( 1 ) ) )
                            {
                                tEdges.emplace_back( iFirst, iSecond );
                            }
                        }
                        tEdges.emplace_back( iFirst, 4 );
                    }
                    for ( const auto& [ iFirst, iSecond ] : tEdges )
                    {
                        uint             tSecondNode = iSecond < 4 ? (uint)tNodes( iSecond ) : tCentreIndex;
                        Matrix< DDRMat > tSecondXi   = iSecond < 4 ? tXi( iSecond ) : Matrix< DDRMat >( 1, 2, 0.0 );
                        bool             tCut        = tGeometryEngine.queue_intersection(
                                tNodes( iFirst ), tSecondNode, tXi( iFirst ), tSecondXi, tNodes, mtk::Geometry_Type::QUAD, mtk::Interpolation_Order::LINEAR );
                        if ( tCut and not tGeometryEngine.queued_intersection_first_parent_on_interface()
                                and not tGeometryEngine.queued_intersection_second_parent_on_interface() )
                        {
                            tLocalCoordinates.push_back( tGeometryEngine.get_queued_intersection_local_coordinate() );
                            aCutEdges.push_back( { tMesh->get_node_coordinate( tNodes( iFirst ) ),
                                    iSecond < 4 ? tMesh->get_node_coordinate( tNodes( iSecond ) ) : tCentre } );
                        }
                    }
                }
            }    // the geometry engine and its geometries (and with them the log) are destroyed here

            delete tMesh;
            return tLocalCoordinates;
        }
    }    // namespace

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "GEN field evaluation log", "[gen], [GEN_EVAL_LOG]" )
    {
        if ( par_size() != 1 )
        {
            return;
        }
        Derived_Node::clear_override_linear_interpolation();

        const std::string tLogFile  = "GEN_eval_log_test.csv";
        const std::string tRankFile = Field_Evaluation_Log::rank_file_name( tLogFile );
        CHECK( tRankFile == "GEN_eval_log_test_rank0.csv" );
        std::remove( tRankFile.c_str() );
        const std::string tBoxLogFile  = "GEN_eval_log_test_box.csv";
        const std::string tBoxRankFile = Field_Evaluation_Log::rank_file_name( tBoxLogFile );
        std::remove( tBoxRankFile.c_str() );
        std::remove( ( tBoxRankFile + ".summary.csv" ).c_str() );
        std::remove( ( tRankFile + ".summary.csv" ).c_str() );

        // default off: nothing is written, and the cut is the reference
        std::vector< Cut_Edge > tOffEdges;
        uint                    tOffNodes = 0;
        Vector< real >          tOff      = cut_circle( "", "", tOffEdges, tOffNodes );
        REQUIRE( tOff.size() > 0 );
        CHECK( not std::ifstream( tRankFile ).good() );
        CHECK( not std::ifstream( tRankFile + ".summary.csv" ).good() );

        SECTION( "whole domain: every cut edge has root iterates on it, values unchanged" )
        {
            std::vector< Cut_Edge > tEdges;
            uint                    tNumNodes = 0;
            Vector< real >          tOn       = cut_circle( tLogFile, "", tEdges, tNumNodes );

            // logging does not change the cut
            REQUIRE( tOn.size() == tOff.size() );
            for ( uint iCut = 0; iCut < tOn.size(); iCut++ )
            {
                CHECK( tOn( iCut ) == tOff( iCut ) );
            }

            std::vector< Log_Row > tRows = read_log( tRankFile );
            REQUIRE( tRows.size() > 0 );

            // background nodes: one row per node in the dedicated pass
            uint tBackground = 0;
            for ( const Log_Row& iRow : tRows )
            {
                tBackground += iRow.mStage == "ut_background_nodes" and iRow.mPurpose == "background_node";
            }
            CHECK( tBackground == tNumNodes );
            CHECK( tNumNodes == 25 );

            // centre nodes are evaluated as derived nodes of kind centre_node
            uint tCentre = 0;
            for ( const Log_Row& iRow : tRows )
            {
                tCentre += iRow.mPurpose == "centre_node";
            }
            CHECK( tCentre > 0 );

            // every cut edge has root iterates on it, from the analytic root
            for ( const Cut_Edge& iEdge : tEdges )
            {
                uint tOnEdge = 0;
                real tDx     = iEdge.mSecond( 0 ) - iEdge.mFirst( 0 );
                real tDy     = iEdge.mSecond( 1 ) - iEdge.mFirst( 1 );
                real tLength = std::sqrt( tDx * tDx + tDy * tDy );
                for ( const Log_Row& iRow : tRows )
                {
                    if ( iRow.mPurpose != "root_iterate" )
                    {
                        continue;
                    }
                    CHECK( iRow.mOperation == "analytic_root" );
                    CHECK( iRow.mNodeIndex == -1 );
                    real tCross = ( iRow.mX - iEdge.mFirst( 0 ) ) * tDy - ( iRow.mY - iEdge.mFirst( 1 ) ) * tDx;
                    real tAlong = ( ( iRow.mX - iEdge.mFirst( 0 ) ) * tDx + ( iRow.mY - iEdge.mFirst( 1 ) ) * tDy ) / ( tLength * tLength );
                    tOnEdge += std::abs( tCross ) < 1e-12 * tLength and tAlong > 0.0 and tAlong < 1.0;
                }
                CHECK( tOnEdge > 0 );
            }

            // the root search converges onto the circle: the last iterates have |phi| below the tolerance
            real tMinAbs = 1.0;
            for ( const Log_Row& iRow : tRows )
            {
                if ( iRow.mPurpose == "root_iterate" )
                {
                    tMinAbs = std::min( tMinAbs, std::abs( iRow.mValue ) );
                }
            }
            CHECK( tMinAbs < 1e-12 );

            // the summary counts the rows written
            std::ifstream tSummary( tRankFile + ".summary.csv" );
            REQUIRE( tSummary.good() );
            std::string tLine;
            std::getline( tSummary, tLine );
            luint tSumAll = 0, tSumBox = 0;
            while ( std::getline( tSummary, tLine ) )
            {
                std::vector< std::string > tItems;
                std::stringstream          tStream( tLine );
                std::string                tItem;
                while ( std::getline( tStream, tItem, ',' ) )
                {
                    tItems.push_back( tItem );
                }
                REQUIRE( tItems.size() == 9 );
                tSumAll += std::stoul( tItems[ 4 ] );
                tSumBox += std::stoul( tItems[ 5 ] );
            }
            CHECK( tSumBox == tRows.size() );
            CHECK( tSumAll == tRows.size() );
        }

        SECTION( "box: only evaluations inside it are written" )
        {
            std::vector< Cut_Edge > tEdges;
            uint                    tNumNodes = 0;
            cut_circle( tBoxLogFile, "-0.6,0.1,-1.0,0.6", tEdges, tNumNodes );

            std::vector< Log_Row > tRows       = read_log( tBoxRankFile );
            uint                   tBackground = 0;
            uint                   tRoots      = 0;
            for ( const Log_Row& iRow : tRows )
            {
                CHECK( iRow.mX >= -0.6 );
                CHECK( iRow.mX <= 0.1 );
                CHECK( iRow.mY >= -1.0 );
                CHECK( iRow.mY <= 0.6 );
                tBackground += iRow.mStage == "ut_background_nodes" and iRow.mPurpose == "background_node";
                tRoots += iRow.mPurpose == "root_iterate";
            }
            CHECK( tBackground == tNumNodes );
            CHECK( tNumNodes == 8 );    // x in { -0.5, 0 }, y in { -1, -0.5, 0, 0.5 }
            CHECK( tRoots > 0 );

            // the summary also counts the evaluations outside the box
            std::ifstream tSummary( tBoxRankFile + ".summary.csv" );
            REQUIRE( tSummary.good() );
            std::string tLine;
            std::getline( tSummary, tLine );
            luint tSumAll = 0, tSumBox = 0;
            while ( std::getline( tSummary, tLine ) )
            {
                std::vector< std::string > tItems;
                std::stringstream          tStream( tLine );
                std::string                tItem;
                while ( std::getline( tStream, tItem, ',' ) )
                {
                    tItems.push_back( tItem );
                }
                tSumAll += std::stoul( tItems[ 4 ] );
                tSumBox += std::stoul( tItems[ 5 ] );
            }
            CHECK( tSumBox == tRows.size() );
            CHECK( tSumAll > tSumBox );
        }

        std::remove( tRankFile.c_str() );
        std::remove( tBoxRankFile.c_str() );
        std::remove( ( tBoxRankFile + ".summary.csv" ).c_str() );
        std::remove( ( tRankFile + ".summary.csv" ).c_str() );
    }
}    // namespace moris::gen
