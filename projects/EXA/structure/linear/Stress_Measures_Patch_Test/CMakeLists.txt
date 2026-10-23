#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# case 0: closed-form checks of the stress measures; cases 1 and 2: one CM named on both sides vs two copies;
# case 3: MAX_NORMAL_STRESS (one defect per case, see the deck header)
moris_add_example(
    NAME       Stress_Measures_Patch_Test
    TEST_BASE  Stress_Measures_Patch_Test
    PROCS      1
    CASES      0 1 2 3
)
