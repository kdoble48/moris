#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Stefans_Problem_conform
    TEST_BASE  Comsol_conform       # legacy ctest name is exactly Comsol_conform (no -procs suffix)
    PROCS      1
    NO_PROCS_SUFFIX
    # INPUTS defaults to TEST_BASE: dynamic_link_input(Comsol_conform Comsol_conform Comsol_conform.cpp ...) as before
)
