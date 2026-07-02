#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Field_Example
    TEST_BASE  Field_Example    # ctest name stays Field_Example-1-procs
    PROCS      1
    INPUTS     Field_example_write Field_example_read Field_example_compare
)
