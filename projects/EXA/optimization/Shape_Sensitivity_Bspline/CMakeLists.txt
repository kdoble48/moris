#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Shape_Sensitivity_Bspline
    TEST_BASE  Shape_Sensitivity_Bspline    # ctest name stays Shape_Sensitivity_Bspline-1-procs
    PROCS      1
    INPUTS     Shape_Sensitivity_Bspline_2D Shape_Sensitivity_Bspline_3D
)
