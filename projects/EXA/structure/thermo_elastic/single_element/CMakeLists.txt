#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       single_element
    TEST_BASE  single_element
    PROCS      1
    NO_PROCS_SUFFIX    # legacy ctest name is exactly "single_element" (no -1-procs suffix)
)
