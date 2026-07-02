#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Fickian_Problem
    TEST_BASE  Fick_Problem       # legacy ctest name is exactly Fick_Problem (no -procs suffix)
    PROCS      1
    NO_PROCS_SUFFIX
    # INPUTS defaults to TEST_BASE: dynamic_link_input(Fick_Problem Fick_Problem Fick_Problem.cpp ...) as before
)
