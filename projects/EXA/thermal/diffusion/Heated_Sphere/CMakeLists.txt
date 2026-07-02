#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Heated_Sphere
    TEST_BASE  Heated_Sphere    # ctest names stay Heated_Sphere-1-procs / Heated_Sphere-4-procs
    PROCS      1 4
)
