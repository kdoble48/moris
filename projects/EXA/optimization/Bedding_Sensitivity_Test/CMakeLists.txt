#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Bedding_Sensitivity_Test
    TEST_BASE  Bedding_Sensitivity_Test    # ctest names stay Bedding_Sensitivity_Test-{1,4}-procs
    PROCS      1 4
)
