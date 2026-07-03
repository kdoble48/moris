#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
# NOTE: this leaf is NOT referenced by the parent CMakeLists (add_subdirectory is
# commented out there); its deck predates the modern parameter-list API. Converted
# to the single-runner pattern so reconnection is a one-line parent change once
# the deck is modernized (share/doc/EXA_RUNNER_RFC.md owner decision 2).
moris_add_example(
    NAME       Homogenization_3D
    TEST_BASE  Homogenization_3D    # legacy add_test name was exactly Homogenization_3D (no -procs suffix)
    PROCS      1
    NO_PROCS_SUFFIX
)
