#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Laplace_2D
    TEST_BASE  Laplace          # ctest names stay Laplace-1-procs / Laplace-4-procs
    PROCS      1 4
    # INPUTS defaults to TEST_BASE: dynamic_link_input(Laplace Laplace Laplace.cpp ...) as before
)
