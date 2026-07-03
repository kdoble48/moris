#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
# Consolidated EXA runner (see share/doc/EXA_RUNNER_RFC.md): one executable containing the
# test cases of every migrated example leaf. Leaves register their TUs via the
# EXA_RUNNER_SOURCES global property (moris_add_example.cmake); this directory
# must be added LAST in projects/EXA/CMakeLists.txt.

get_property(EXA_RUNNER_SOURCES GLOBAL PROPERTY EXA_RUNNER_SOURCES)
if(NOT EXA_RUNNER_SOURCES)
    return()   # nothing migrated yet: no runner target, no cost
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Test-case TUs are listed directly in add_executable - never move them into an
# intermediate static library, which would let the linker drop Catch2's
# self-registration objects.
add_executable(EXA-test
    runner_main.cpp
    EXA_Globals.cpp
    ${EXA_RUNNER_SOURCES})

target_include_directories(EXA-test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})   # EXA_Globals.hpp

# Identical to the dependency list every legacy example leaf uses.
target_link_libraries(EXA-test PRIVATE
    ${HMR}-lib
    ${COM}-lib
    ${WRK}-lib
    ${MORIS_PETSC_LIBS}
    ${MORIS_BOOST_LIBS}
    ${MORIS_ACML_LAPACK_MKL_OPENBLAS_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_LDLIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_BASE_LIBS})

set_target_properties(EXA-test PROPERTIES OUTPUT_NAME EXA-test${EXE_EXT})
