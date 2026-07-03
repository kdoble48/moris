#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Couette_Flow_Static
    TEST_BASE  Couette_Flow_Static    # ctest names stay Couette_Flow_Static-1-procs / Couette_Flow_Static-2-procs
    PROCS      1 2
    # INPUTS defaults to TEST_BASE: dynamic_link_input(Couette_Flow_Static ...) as before
)
