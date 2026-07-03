#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Voxel_Grain_Structure
    TEST_BASE  Voxel_Grain_Structure    # ctest name stays Voxel_Grain_Structure-1-procs
    PROCS      1
    # INPUTS defaults to TEST_BASE: dynamic_link_input(Voxel_Grain_Structure Voxel_Grain_Structure Voxel_Grain_Structure.cpp ...) as before
)
