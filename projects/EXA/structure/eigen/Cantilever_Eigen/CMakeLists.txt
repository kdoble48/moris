#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Cantilever_Eigen
    TEST_BASE  Cantilever_Eigen    # ctest name stays Cantilever_Eigen-1-procs
    PROCS      1
    # INPUTS defaults to TEST_BASE: dynamic_link_input(Cantilever_Eigen ...) as before
)
