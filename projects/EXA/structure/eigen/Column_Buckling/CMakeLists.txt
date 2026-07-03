#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Column_Buckling
    TEST_BASE  Column_Buckling    # ctest names stay Column_Buckling-1-procs / Column_Buckling-4-procs
    PROCS      1 4
    INPUTS     Column_Buckling Column_Buckling_Immersed
)
