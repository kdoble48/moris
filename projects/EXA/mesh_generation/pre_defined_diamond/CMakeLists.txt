#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       pre_defined_diamond
    TEST_BASE  Pre_Defined_Diamond    # ctest name stays Pre_Defined_Diamond-1-procs
    PROCS      1
    NO_INPUTS                         # meshgen example: --meshgen XML input, no deck .so
)
