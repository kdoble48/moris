#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
# NOTE: this leaf is NOT referenced by the parent CMakeLists (no add_subdirectory
# reaches it); its deck predates the modern parameter-list API. Converted to the
# single-runner pattern so reconnection is a one-line parent change once the deck
# is modernized (share/doc/EXA_RUNNER_RFC.md owner decision 2).
moris_add_example(
    NAME       Plane_Strain
    TEST_BASE  Plane_Strain_Problem    # legacy add_test name was exactly Plane_Strain_Problem (no -procs suffix)
    PROCS      1
    NO_PROCS_SUFFIX
    # INPUTS defaults to TEST_BASE: dynamic_link_input(Plane_Strain_Problem ...) as before
)
