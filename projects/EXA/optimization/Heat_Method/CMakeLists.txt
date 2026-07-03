#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Heat_Method
    TEST_BASE  HeatMethod    # ctest name stays HeatMethod-1-procs
    PROCS      1
    # INPUTS defaults to TEST_BASE: dynamic_link_input(HeatMethod HeatMethod HeatMethod.cpp ...) as before
)
