#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Beam_Temperature_Field
    TEST_BASE  Beam_Temperature_Field    # ctest name stays Beam_Temperature_Field-1-procs
    PROCS      1
    # INPUTS defaults to TEST_BASE: dynamic_link_input(Beam_Temperature_Field ...) as before
)
