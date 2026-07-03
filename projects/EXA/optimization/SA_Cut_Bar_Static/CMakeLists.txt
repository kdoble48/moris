#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       SA_Cut_Bar_Static
    TEST_BASE  SA_Cut_Bar_Static    # ctest name stays SA_Cut_Bar_Static-1-procs
    PROCS      1
    # INPUTS defaults to TEST_BASE: dynamic_link_input(SA_Cut_Bar_Static ...) as before
)
