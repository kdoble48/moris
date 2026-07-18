#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# IOS Tests ---------------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_Logger.cpp
    fn_to_stdio.cpp
    UT_IOS_Parsing_Tool.cpp
    UT_IOS_File_To_Array.cpp
    UT_IOS_Performance_Reporter.cpp
    UT_IOS_Library_IO_Deck_Semantics.cpp
    UT_IOS_Input_Deck_Typed_Layer.cpp )

# Fixture input-deck shared objects for the Library_IO deck-semantics tests.
# Like real input decks they are NOT linked against any moris libraries; their moris
# symbols resolve at dlopen against the (-rdynamic) test executable. The include
# set mirrors the deck-.so convention of dynamic_link_input()/moris_add_example.
set(IOS_FIXTURE_INCLUDES
    ${CMAKE_BINARY_DIR}/generated
    ${MORIS_PACKAGE_DIR}/PRM/src
    ${MORIS_PACKAGE_DIR}/PRM/ENM/src
    ${MORIS_PACKAGE_DIR}/COM/src
    ${MORIS_PACKAGE_DIR}/LINALG/src
    ${MORIS_PACKAGE_DIR}/LINALG/src/${LINALG_IMPLEMENTATION_INCLUDES}
    ${MORIS_BOOST_INCLUDE_DIRS}
    ${MORIS_ARMADILLO_INCLUDE_DIRS})

foreach(FIXTURE IOS-deck-fixture,deck_semantics_fixture
                IOS-deck-v2-fixture,deck_entry_point_fixture
                IOS-deck-mixed-fixture,deck_mixed_fixture)
    string(REPLACE "," ";" FIXTURE ${FIXTURE})
    list(GET FIXTURE 0 FIXTURE_TARGET)
    list(GET FIXTURE 1 FIXTURE_SOURCE)
    add_library(${FIXTURE_TARGET} SHARED fixtures/${FIXTURE_SOURCE}.cpp)
    target_include_directories(${FIXTURE_TARGET} PRIVATE ${IOS_FIXTURE_INCLUDES})
    set_target_properties(${FIXTURE_TARGET} PROPERTIES POSITION_INDEPENDENT_CODE ON)
endforeach()

# List additional includes
include_directories(${MORIS_DIR}/snippets/ios)
include_directories(${CMAKE_BINARY_DIR}/generated/)

# List test dependencies
set(TEST_DEPENDENCIES
    test-libs
    )


# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name for general IOS test
set(IOS_TEST ${IOS}-test)

# Create executable
add_executable(${IOS_TEST} ${TEST_SOURCES})
target_link_libraries(${IOS_TEST} ${TEST_DEPENDENCIES})
set_target_properties(${IOS_TEST} PROPERTIES OUTPUT_NAME ${IOS_TEST}${EXE_EXT})

# Hand the deck-semantics tests the paths to the fixture .so's and make sure they build
target_compile_definitions(${IOS_TEST} PRIVATE
    "IOS_DECK_FIXTURE_SO=\"$<TARGET_FILE:IOS-deck-fixture>\""
    "IOS_DECK_V2_FIXTURE_SO=\"$<TARGET_FILE:IOS-deck-v2-fixture>\""
    "IOS_DECK_MIXED_FIXTURE_SO=\"$<TARGET_FILE:IOS-deck-mixed-fixture>\"")
add_dependencies(${IOS_TEST} IOS-deck-fixture IOS-deck-v2-fixture IOS-deck-mixed-fixture)


# List source files for Query
set(QUERY_SOURCES
    query_main.cpp
    cl_query_test.cpp )

# Set test name for Query
set(IOS_QUERY ${IOS}-query)

# Create executable for Query
add_executable(${IOS_QUERY} ${QUERY_SOURCES})
target_link_libraries(${IOS_QUERY} ${TEST_DEPENDENCIES})
set_target_properties(${IOS_QUERY} PROPERTIES OUTPUT_NAME ${IOS_QUERY}${EXE_EXT})


# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(IOS_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${IOS_TEST_PROCS})
		set(IOS_PARALLEL_TEST ${IOS_TEST}-${PROCS}-procs)
		add_test(NAME ${IOS_PARALLEL_TEST} 
			COMMAND ${MORIS_EXECUTE_COMMAND} -n ${PROCS} ${VALGRIND} ${VALGRIND_OPTIONS} $<TARGET_FILE:${IOS_TEST}>)
	endforeach()
endif()

