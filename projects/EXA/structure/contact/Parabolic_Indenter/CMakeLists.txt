#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Parabolic_Indenter
    TEST_BASE  Parabolic_Indentor    # legacy ctest name typo ("Indentor") preserved - name parity is the hard constraint
    PROCS      1
    INPUTS     Parabolic_Indenter_Linear Parabolic_Indenter_NonLinear
)
