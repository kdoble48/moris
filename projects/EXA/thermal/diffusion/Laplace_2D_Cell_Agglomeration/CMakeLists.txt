#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Laplace_2D_Cell_Agglomeration
    TEST_BASE  Laplace_Cell_Agglomeration    # ctest name: Laplace_Cell_Agglomeration-1-procs
    PROCS      1
)
