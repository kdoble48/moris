#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Channel_2D_Static
    TEST_BASE  Channel_2D_Static    # ctest names stay Channel_2D_Static-1-procs / Channel_2D_Static-2-procs
    PROCS      1 2
    # INPUTS defaults to TEST_BASE: dynamic_link_input(Channel_2D_Static ...) as before
)
