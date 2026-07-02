#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Stefans_Problem_cut
    TEST_BASE  Comsol_cut       # legacy ctest name is exactly Comsol_cut (no -procs suffix)
    PROCS      1
    NO_PROCS_SUFFIX
    # INPUTS defaults to TEST_BASE: dynamic_link_input(Comsol_cut Comsol_cut Comsol_cut.cpp ...) as before
)
