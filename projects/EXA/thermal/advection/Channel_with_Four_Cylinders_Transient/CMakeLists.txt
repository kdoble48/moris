#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Channel_with_Four_Cylinders_Transient
    TEST_BASE  Channel_with_Four_Cylinders_Transient    # legacy add_test used the bare name, hardcoded -n 1
    PROCS      1
    NO_PROCS_SUFFIX
)
