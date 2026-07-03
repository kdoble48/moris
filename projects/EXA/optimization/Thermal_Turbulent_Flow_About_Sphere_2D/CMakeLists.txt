#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Thermal_Turbulent_Flow_About_Sphere_2D
    TEST_BASE  ChannelTurbulent2D    # ctest names stay ChannelTurbulent2D-{1,3}-procs
    PROCS      1 3
    # INPUTS defaults to TEST_BASE: dynamic_link_input(ChannelTurbulent2D ...) as before
)
