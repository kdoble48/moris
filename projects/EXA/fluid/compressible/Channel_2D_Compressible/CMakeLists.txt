#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Channel_2D_Compressible
    TEST_BASE  Channel_2D_Compressible    # ctest name stays Channel_2D_Compressible-1-procs
    PROCS      1
    # INPUTS defaults to TEST_BASE: dynamic_link_input(Channel_2D_Compressible ...) as before
)
