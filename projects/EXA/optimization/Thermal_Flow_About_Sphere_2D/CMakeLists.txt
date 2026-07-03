#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Thermal_Flow_About_Sphere_2D
    TEST_BASE  Thermal_Flow_About_Sphere_2D    # ctest names stay Thermal_Flow_About_Sphere_2D-{1,3}-procs
    PROCS      1 3
)
