#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Shape_Sensitivity_Sweep
    TEST_BASE  Shape_Sensitivity_Sweep    # ctest name stays Shape_Sensitivity_Sweep-1-procs
    PROCS      1
    # INPUTS defaults to TEST_BASE: dynamic_link_input(Shape_Sensitivity_Sweep ...) as before
)
