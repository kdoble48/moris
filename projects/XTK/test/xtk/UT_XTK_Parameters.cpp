/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_XTK_Parameters.cpp
 *
 */

#include "catch.hpp"

#include "cl_Parameter_List.hpp"
#include "fn_PRM_XTK_Parameters.hpp"

#define protected public
#define private public
#include "cl_XTK_Integration_Mesh_Generator.hpp"
#undef protected
#undef private

namespace moris::xtk
{
    // ----------------------------------------------------------------------------------

    TEST_CASE( "XTK Parameter List Defaults", "[XTK],[XTK_Parameters]" )
    {
        // create the default XTK parameter list
        moris::Parameter_List tXTKParameters = moris::prm::create_xtk_parameter_list();

        SECTION( "output_cut_ig_mesh defaults to false" )
        {
            REQUIRE( tXTKParameters.get< bool >( "output_cut_ig_mesh" ) == false );
        }

        SECTION( "keep_all_opt_iters defaults to false" )
        {
            REQUIRE( tXTKParameters.get< bool >( "keep_all_opt_iters" ) == false );
        }

        SECTION( "output_path defaults to ./" )
        {
            REQUIRE( tXTKParameters.get< std::string >( "output_path" ) == "./" );
        }

        SECTION( "parameters can be set" )
        {
            tXTKParameters.set( "output_cut_ig_mesh", true );
            tXTKParameters.set( "keep_all_opt_iters", true );
            tXTKParameters.set( "output_path", std::string( "/tmp/test/" ) );

            REQUIRE( tXTKParameters.get< bool >( "output_cut_ig_mesh" ) == true );
            REQUIRE( tXTKParameters.get< bool >( "keep_all_opt_iters" ) == true );
            REQUIRE( tXTKParameters.get< std::string >( "output_path" ) == "/tmp/test/" );
        }
    }

    // ----------------------------------------------------------------------------------

    TEST_CASE( "Integration_Mesh_Generator Default State", "[XTK],[XTK_Parameters],[XTK_CutIgMesh]" )
    {
        Integration_Mesh_Generator tMeshGen;

        SECTION( "mOutputCutIgMesh defaults to false" )
        {
            REQUIRE( tMeshGen.mOutputCutIgMesh == false );
        }

        SECTION( "set_cut_IG_mesh_output setter works" )
        {
            tMeshGen.set_cut_IG_mesh_output( true );
            REQUIRE( tMeshGen.get_cut_IG_mesh_output() == true );
        }
    }

    // ----------------------------------------------------------------------------------

}    // namespace moris::xtk
