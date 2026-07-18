/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_IOS_Library_IO_Deck_Semantics.cpp
 *
 * Safety-net tests pinning the CURRENT input-deck loading contract of Library_IO
 * before any deck-API redesign work. These tests document behavior as-is; several
 * pinned behaviors are surprising (and candidates for later change), but any change
 * to them must be deliberate and show up as a diff to this file.
 *
 * Pinned contract:
 *  (a) a MISSING <MODULE>ParameterList symbol CLEARS that module's parameter lists
 *  (b) a PRESENT-BUT-EMPTY <MODULE>ParameterList function KEEPS the fn_PRM defaults
 *  (c) load_function(): finalize-gate, throw/no-throw on missing symbols, and the
 *      guard against loading *ParameterList symbols as user functions
 *  (d) finalize(): .so applied first, then XML on top; the library locks afterwards;
 *      XML presence RESETS modules that are absent from the XML file to defaults
 *      (including modules the .so had set or cleared)
 *
 * Uses the fixture deck fixtures/deck_semantics_fixture.cpp, compiled to a shared
 * object by CMake and passed in via the IOS_DECK_FIXTURE_SO compile definition.
 */

#include <catch.hpp>

#include <filesystem>
#include <fstream>
#include <string>

#include "cl_Library_IO_Standard.hpp"
#include "cl_Communication_Tools.hpp"
#include "parameters.hpp"

#ifndef IOS_DECK_FIXTURE_SO
#error "IOS_DECK_FIXTURE_SO must be defined by CMake (path to the fixture deck .so)"
#endif

namespace moris
{
    //------------------------------------------------------------------------------------------------------------------

    // signature of the fixture's user callback
    typedef int ( *Fixture_User_Function )();

    //------------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Library_IO deck semantics: missing vs empty module symbols", "[IOS],[Library_IO],[deck_semantics]" )
    {
        Library_IO_Standard tLibrary;
        tLibrary.load_parameter_list( IOS_DECK_FIXTURE_SO, File_Type::SO_FILE );
        tLibrary.finalize( "" );

        SECTION( "present-and-mutating function: .so parameters are applied on top of defaults" )
        {
            auto tOptParams = tLibrary.get_parameters_for_module( Module_Type::OPT );

            // the ctor seeds exactly one optimization-problem parameter list, no algorithms
            REQUIRE( tOptParams( OPT::OPTIMIZATION_PROBLEMS ).size() == 1 );
            CHECK( tOptParams( OPT::ALGORITHMS ).size() == 0 );

            // the fixture flipped this from its default (false)
            CHECK( tOptParams( OPT::OPTIMIZATION_PROBLEMS )( 0 ).get< bool >( "is_optimization_problem" ) == true );
        }

        SECTION( "present-but-empty function: fn_PRM defaults are preserved" )
        {
            auto tHmrParams = tLibrary.get_parameters_for_module( Module_Type::HMR );

            REQUIRE( tHmrParams( HMR::GENERAL ).size() == 1 );
            CHECK( tHmrParams( HMR::GENERAL )( 0 ).get< uint >( "refinement_buffer" ) == 0 );
        }

        SECTION( "missing symbol: the module's parameter lists are cleared" )
        {
            // NOTE: this is the documented trap — omitting a symbol DISABLES the module,
            // while an empty function keeps it enabled with defaults (section above)
            CHECK( tLibrary.get_parameters_for_module( Module_Type::GEN ).size() == 0 );
            CHECK( tLibrary.get_parameters_for_module( Module_Type::FEM ).size() == 0 );
        }
    }

    //------------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Library_IO deck semantics: load_function resolution", "[IOS],[Library_IO],[deck_semantics]" )
    {
        Library_IO_Standard tLibrary;
        tLibrary.load_parameter_list( IOS_DECK_FIXTURE_SO, File_Type::SO_FILE );

        SECTION( "loading a function before finalize() throws" )
        {
            REQUIRE_THROWS( tLibrary.load_function< Fixture_User_Function >( "Deck_Fixture_User_Function" ) );
        }

        SECTION( "after finalize(): existing symbols resolve, missing symbols throw or return nullptr" )
        {
            tLibrary.finalize( "" );

            // an existing user callback resolves and is callable
            Fixture_User_Function tFunc = tLibrary.load_function< Fixture_User_Function >( "Deck_Fixture_User_Function" );
            REQUIRE( tFunc != nullptr );
            CHECK( tFunc() == 42 );

            // a missing symbol throws by default ...
            REQUIRE_THROWS( tLibrary.load_function< Fixture_User_Function >( "This_Function_Does_Not_Exist" ) );

            // ... and returns nullptr when aThrowError = false
            CHECK( tLibrary.load_function< Fixture_User_Function >( "This_Function_Does_Not_Exist", false ) == nullptr );

            // *ParameterList symbols may not be loaded as user functions
            REQUIRE_THROWS( tLibrary.load_function< Parameter_Function >( "OPTParameterList" ) );
        }
    }

    //------------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Library_IO deck semantics: finalize order and locking", "[IOS],[Library_IO],[deck_semantics]" )
    {
        SECTION( "finalize() locks the library against further inputs" )
        {
            Library_IO_Standard tLibrary;
            tLibrary.load_parameter_list( IOS_DECK_FIXTURE_SO, File_Type::SO_FILE );
            tLibrary.finalize( "" );

            REQUIRE_THROWS( tLibrary.load_parameter_list( IOS_DECK_FIXTURE_SO, File_Type::SO_FILE ) );
        }

        SECTION( "finalize( path ) writes a parameter receipt" )
        {
            std::filesystem::path tReceiptPath =
                    std::filesystem::temp_directory_path()
                    / ( "ios_deck_semantics_receipt_" + std::to_string( par_rank() ) + ".xml" );
            std::filesystem::remove( tReceiptPath );

            Library_IO_Standard tLibrary;
            tLibrary.load_parameter_list( IOS_DECK_FIXTURE_SO, File_Type::SO_FILE );
            tLibrary.finalize( tReceiptPath.string() );

            CHECK( std::filesystem::exists( tReceiptPath ) );
            std::filesystem::remove( tReceiptPath );
        }
    }

    //------------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Library_IO deck semantics: XML input semantics on top of a .so", "[IOS],[Library_IO],[deck_semantics]" )
    {
        // per-rank file to avoid concurrent writes under mpirun
        std::filesystem::path tXmlPath =
                std::filesystem::temp_directory_path()
                / ( "ios_deck_semantics_override_" + std::to_string( par_rank() ) + ".xml" );

        // XML overriding a single HMR parameter; no other modules are mentioned
        {
            std::ofstream tXmlFile( tXmlPath );
            tXmlFile << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
                        "<ParameterLists>\n"
                        "    <HMR>\n"
                        "        <GENERAL>\n"
                        "            <General>\n"
                        "                <refinement_buffer>3</refinement_buffer>\n"
                        "            </General>\n"
                        "        </GENERAL>\n"
                        "    </HMR>\n"
                        "</ParameterLists>\n";
        }

        Library_IO_Standard tLibrary;
        tLibrary.load_parameter_list( IOS_DECK_FIXTURE_SO, File_Type::SO_FILE );
        tLibrary.load_parameter_list( tXmlPath.string(), File_Type::XML_FILE );
        tLibrary.finalize( "" );

        SECTION( "XML values are applied after (on top of) the .so pass — but APPENDED behind a fresh default" )
        {
            // NOTE: pinned as-is. The XML pass builds a fresh Module_Parameter_Lists (whose
            // constructor seeds one default list per submodule) and APPENDS the XML-parsed
            // list to it. The XML-set value therefore lands at index 1, behind an untouched
            // default at index 0 — consumers that read index 0 never see XML values.
            auto tHmrParams = tLibrary.get_parameters_for_module( Module_Type::HMR );

            REQUIRE( tHmrParams( HMR::GENERAL ).size() == 2 );
            CHECK( tHmrParams( HMR::GENERAL )( 0 ).get< uint >( "refinement_buffer" ) == 0 );    // ctor-seeded default
            CHECK( tHmrParams( HMR::GENERAL )( 1 ).get< uint >( "refinement_buffer" ) == 3 );    // XML-parsed list
        }

        SECTION( "modules absent from the XML are RESET to defaults, discarding .so values" )
        {
            // NOTE: pinned as-is. The XML pass rebuilds every module that is not named in
            // the XML file from defaults — so the value set by the .so's OPTParameterList
            // is LOST, and the GEN module the .so pass had cleared is re-created with
            // default parameter lists. This makes .so + XML composition module-granular,
            // not key-granular. Deliberately surprising; candidates for redesign.
            auto tOptParams = tLibrary.get_parameters_for_module( Module_Type::OPT );

            REQUIRE( tOptParams( OPT::OPTIMIZATION_PROBLEMS ).size() >= 1 );
            CHECK( tOptParams( OPT::OPTIMIZATION_PROBLEMS )( 0 ).get< bool >( "is_optimization_problem" ) == false );

            // GEN was cleared by the .so pass (missing symbol) but is re-created by the XML pass
            CHECK( tLibrary.get_parameters_for_module( Module_Type::GEN ).size() > 0 );
        }

        std::filesystem::remove( tXmlPath );
    }

    //------------------------------------------------------------------------------------------------------------------

}    // namespace moris
