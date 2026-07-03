#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Axisymmetric
    TEST_BASE  Axisymmetric_Problem    # ctest name stays Axisymmetric_Problem-1-procs
    PROCS      1
    INPUTS     Axisymmetric_Problem Axisymmetric_Problem_ConstantThermalField
)
