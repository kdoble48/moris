#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
# moris_add_example(
#     NAME <ExampleDirId>            # unique; Catch2 tag becomes [EXA_<NAME>]
#     TEST_BASE <name>               # ctest base name == legacy EXAMPLE_FILE (e.g. Laplace)
#     PROCS <n> [<n> ...]            # registers <TEST_BASE>-<n>-procs per entry
#     [NO_PROCS_SUFFIX]              # single test named exactly <TEST_BASE>
#     [SOURCES <files> ...]          # default: example_test_case.cpp
#     [INPUTS <deckbase> ...]        # each: dynamic_link_input(<b> <b> <b>.cpp); default: TEST_BASE
#     [EXTRA_SO_INCLUDES <dirs> ...] # for leaves whose SO_INCLUDES deviate from the standard block
# )
#
# Replaces the legacy per-leaf example boilerplate (see doc/internal/EXA_RUNNER_RFC.md).
# Preserves, byte-for-byte, the dynamic_link_input() targets, their output
# names/locations, and the ctest names/PROCS of the legacy pattern. The leaf's
# test-case TU(s) are contributed to the shared EXA-test runner via the
# EXA_RUNNER_SOURCES global property (consumed by projects/EXA/runner/).

function(moris_add_example)
    set(options NO_PROCS_SUFFIX NO_INPUTS)
    set(oneValueArgs NAME TEST_BASE)
    set(multiValueArgs PROCS SOURCES INPUTS EXTRA_SO_INCLUDES)
    cmake_parse_arguments(EXA "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if(NOT EXA_NAME OR NOT EXA_TEST_BASE)
        message(FATAL_ERROR "moris_add_example: NAME and TEST_BASE are required (in ${CMAKE_CURRENT_SOURCE_DIR})")
    endif()
    if(NOT EXA_SOURCES)
        set(EXA_SOURCES example_test_case.cpp)
    endif()
    if(EXA_NO_INPUTS)
        set(EXA_INPUTS)    # deck-less example (e.g. --meshgen XML input)
    elseif(NOT DEFINED EXA_INPUTS)
        set(EXA_INPUTS ${EXA_TEST_BASE})
    endif()

    # ---- input .so machinery: identical to the legacy leaf boilerplate ------
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

    # The tests' WORKING_DIRECTORY must exist even for NO_INPUTS leaves, where no
    # target output creates it (legacy per-example exes used to). Fresh build
    # trees otherwise leave meshgen tests "Not Run".
    file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

    set(SO_TPLS "trilinos" ${ARMADILLO_EIGEN})

    # get list of INT and MTK subfolders
    get_property(INT_SRC_LIST GLOBAL PROPERTY INT_SRC_LIST)
    get_property(MTK_SRC_LIST GLOBAL PROPERTY MTK_SRC_LIST)

    set(SO_INCLUDES
        ${CMAKE_BINARY_DIR}/generated    # configure-time paths.hpp; legacy leaves had this dir-scoped
        ${MORIS_BOOST_INCLUDE_DIRS}
        ${MORIS_TRILINOS_INCLUDE_DIRS}
        ${MORIS_PACKAGE_DIR}/ALG/src
        ${INT_SRC_LIST}
        ${MTK_SRC_LIST}
        ${MORIS_PACKAGE_DIR}/FEM/MSI/src
        ${MORIS_PACKAGE_DIR}/FEM/VIS/src
        ${MORIS_PACKAGE_DIR}/SOL/DLA/src
        ${MORIS_PACKAGE_DIR}/SOL/TSA/src
        ${MORIS_PACKAGE_DIR}/SOL/NLA/src
        ${MORIS_PACKAGE_DIR}/SOL/SOL_CORE/src
        ${MORIS_PACKAGE_DIR}/LINALG/src
        ${MORIS_PACKAGE_DIR}/LINALG/src/${LINALG_IMPLEMENTATION_INCLUDES}
        ${MORIS_PACKAGE_DIR}/COM/src
        ${MORIS_PACKAGE_DIR}/PRM/ENM/src
        ${MORIS_PACKAGE_DIR}/MTK/src
        ${MORIS_PACKAGE_DIR}/HMR/src
        ${MORIS_PACKAGE_DIR}/XTK/src
        ${MORIS_PACKAGE_DIR}/PRM/src
        ${MORIS_PACKAGE_DIR}/MRS/COR/src
        ${EXA_EXTRA_SO_INCLUDES})

    foreach(TPL ${SO_TPLS})
        string(TOUPPER ${TPL} TPL)
        list(APPEND SO_INCLUDES ${MORIS_${TPL}_INCLUDE_DIRS})
    endforeach()

    foreach(INPUT ${EXA_INPUTS})
        dynamic_link_input(${INPUT} ${INPUT} ${INPUT}.cpp ${SO_INCLUDES})
    endforeach()

    # ---- register this leaf's test-case TU(s) with the shared runner --------
    foreach(SRC ${EXA_SOURCES})
        get_filename_component(SRC_ABS ${SRC} ABSOLUTE)
        set_property(GLOBAL APPEND PROPERTY EXA_RUNNER_SOURCES ${SRC_ABS})
    endforeach()

    # ---- ctest registrations: EXACT legacy names, PROCS, workdir, valgrind --
    if(MORIS_HAVE_PARALLEL_TESTS)
        foreach(PROCS ${EXA_PROCS})
            if(EXA_NO_PROCS_SUFFIX)
                set(_test_name ${EXA_TEST_BASE})
            else()
                set(_test_name ${EXA_TEST_BASE}-${PROCS}-procs)
            endif()
            add_test(NAME ${_test_name}
                COMMAND ${MORIS_EXECUTE_COMMAND} -n ${PROCS} ${VALGRIND} ${VALGRIND_OPTIONS_EXA}
                        $<TARGET_FILE:EXA-test> "[EXA_${EXA_NAME}]"
                WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
        endforeach()
    endif()
endfunction()
