#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Oscillator
    TEST_BASE  Oscillator    # ctest name stays Oscillator-1-procs
    PROCS      1
    # INPUTS defaults to TEST_BASE: dynamic_link_input(Oscillator ...) as before;
    # the deck's paths.hpp include is covered by ${CMAKE_BINARY_DIR}/generated in SO_INCLUDES
)
