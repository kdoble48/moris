#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Shape_Sensitivity_Circle_Sweep_Transient
    TEST_BASE  Shape_Sensitivity_Circle_Sweep_Transient    # ctest name stays ...-1-procs
    PROCS      1
    INPUTS     Shape_Sensitivity_Circle_Sweep_Transient Shape_Sensitivity_Circle_Sweep_Transient_Vol_2
)
