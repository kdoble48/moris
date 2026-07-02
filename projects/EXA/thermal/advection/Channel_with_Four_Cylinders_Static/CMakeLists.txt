#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Channel_with_Four_Cylinders_Static
    TEST_BASE  Channel_with_Four_Cylinders_Static    # ctest names stay Channel_with_Four_Cylinders_Static-{1,4}-procs
    PROCS      1 4
)
