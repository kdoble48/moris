#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Grain_Structure
    TEST_BASE  Grain_Structure    # ctest name stays Grain_Structure-1-procs
    PROCS      1
    # INPUTS defaults to TEST_BASE: dynamic_link_input(Grain_Structure Grain_Structure Grain_Structure.cpp ...) as before
)
