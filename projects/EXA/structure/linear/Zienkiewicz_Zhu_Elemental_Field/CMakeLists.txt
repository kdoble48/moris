#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Zienkiewicz_Zhu_Elemental_Field
    TEST_BASE  Zienkiewicz_Zhu_Elemental_Field    # ctest name stays Zienkiewicz_Zhu_Elemental_Field-1-procs
    PROCS      1
    # INPUTS defaults to TEST_BASE: dynamic_link_input(Zienkiewicz_Zhu_Elemental_Field ...) as before
)
