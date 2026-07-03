#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Two_Bar_Truss
    TEST_BASE  Two_Bar_Truss    # ctest name stays Two_Bar_Truss-1-procs
    PROCS      1
    # INPUTS defaults to TEST_BASE: dynamic_link_input(Two_Bar_Truss Two_Bar_Truss Two_Bar_Truss.cpp ...) as before
)
