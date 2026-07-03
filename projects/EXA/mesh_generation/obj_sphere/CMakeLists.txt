#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       obj_sphere
    TEST_BASE  OBJ_Sphere      # ctest name stays OBJ_Sphere-1-procs
    PROCS      1
    NO_INPUTS                  # meshgen example: --meshgen XML input, no deck .so
)
