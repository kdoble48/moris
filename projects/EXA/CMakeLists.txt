#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
include(${MORIS_CMAKE_DIR}/utilities/moris_add_example.cmake)

add_subdirectory(thermal)
add_subdirectory(structure)
add_subdirectory(fluid)
add_subdirectory(optimization)
add_subdirectory(mesh_generation)

# Consolidated EXA runner (share/doc/EXA_RUNNER_RFC.md) - must stay last: consumes the
# EXA_RUNNER_SOURCES global property that migrated leaves append to.
add_subdirectory(runner)
