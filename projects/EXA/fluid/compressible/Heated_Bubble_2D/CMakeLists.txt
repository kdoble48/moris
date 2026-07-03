#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Heated_Bubble_2D
    TEST_BASE  Heated_Bubble_2D    # ctest name stays Heated_Bubble_2D-1-procs
    PROCS      1
    # INPUTS defaults to TEST_BASE: dynamic_link_input(Heated_Bubble_2D ...) as before
)
