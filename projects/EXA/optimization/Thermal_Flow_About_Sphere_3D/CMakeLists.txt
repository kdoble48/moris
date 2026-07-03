#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
# NOTE: this leaf is disconnected (the parent's add_subdirectory is commented out)
# and carries NO example_test_case.cpp in the tree - the legacy CMakeLists already
# referenced a nonexistent TU. Reconnecting requires authoring the test TU first
# (moris_add_example defaults SOURCES to example_test_case.cpp).
moris_add_example(
    NAME       Thermal_Flow_About_Sphere_3D
    TEST_BASE  Thermal_Flow_About_Sphere_3D    # legacy add_test: Thermal_Flow_About_Sphere_3D-{1,4}-procs
    PROCS      1 4
)
