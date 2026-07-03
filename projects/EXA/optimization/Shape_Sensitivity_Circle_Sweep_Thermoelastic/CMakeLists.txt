#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Shape_Sensitivity_Circle_Sweep_Thermoelastic
    TEST_BASE  Shape_Sensitivity_Circle_Sweep_Thermoelastic    # ctest names stay ...-1-procs / ...-2-procs
    PROCS      1 2
    INPUTS     Shape_Sensitivity_Circle_Sweep_Thermoelastic Shape_Sensitivity_Circle_Sweep_Thermoelastic_Staggered Shape_Sensitivity_Circle_Sweep_Thermoelastic_GCMMA
)
