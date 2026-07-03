#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Thick_Pressure_Vessel
    TEST_BASE  Thick_Pressure_Vessel    # ctest names stay Thick_Pressure_Vessel-1-procs / -4-procs
    PROCS      1 4
    INPUTS     Pressure_Vessel_3D Pressure_Vessel_3D_Immersed    # deck basenames differ from TEST_BASE
)
