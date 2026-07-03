#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Thermoelastic_Varying_Order
    TEST_BASE  Bar_Lin_Disp_Quad_Diff    # ctest name stays Bar_Lin_Disp_Quad_Diff-1-procs
    PROCS      1
    # INPUTS defaults to TEST_BASE: dynamic_link_input(Bar_Lin_Disp_Quad_Diff ...) as before
)
