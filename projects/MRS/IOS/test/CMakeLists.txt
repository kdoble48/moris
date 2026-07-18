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
    UT_IOS_Library_IO_Deck_Semantics.cpp )

# Fixture input-deck shared object for the Library_IO deck-semantics tests.
# Like a real input deck it is NOT linked against any moris libraries; its moris
# symbols resolve at dlopen against the (-rdynamic) test executable. The include
# set mirrors the deck-.so convention of dynamic_link_input()/moris_add_example.
add_library(IOS-deck-fixture SHARED fixtures/deck_semantics_fixture.cpp)
target_include_directories(IOS-deck-fixture PRIVATE
    ${CMAKE_BINARY_DIR}/generated
    ${MORIS_PACKAGE_DIR}/PRM/src
    ${MORIS_PACKAGE_DIR}/PRM/ENM/src
    ${MORIS_PACKAGE_DIR}/COM/src
    ${MORIS_PACKAGE_DIR}/LINALG/src
    ${MORIS_PACKAGE_DIR}/LINALG/src/${LINALG_IMPLEMENTATION_INCLUDES}
    ${MORIS_BOOST_INCLUDE_DIRS}
    ${MORIS_ARMADILLO_INCLUDE_DIRS})
set_target_properties(IOS-deck-fixture PROPERTIES POSITION_INDEPENDENT_CODE ON)

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

# Hand the deck-semantics tests the path to the fixture .so and make sure it is built
target_compile_definitions(${IOS_TEST} PRIVATE
    "IOS_DECK_FIXTURE_SO=\"$<TARGET_FILE:IOS-deck-fixture>\"")
add_dependencies(${IOS_TEST} IOS-deck-fixture)


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

