#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
# Level-set B-spline coefficients with and without the Kreissl-Maute filter (adv_filter_radius); one deck, four cases,
# each its own ctest / runner process. Tags: [EXA_LevelSet_ADV_Filter_<0..3>].
#   0, 1: sensitivity sweep of a plate with a hole, without / with the filter (radius 1.5 h)
#   2, 3: 100 MMA iterations of a cantilever, without / with the filter
moris_add_example(
    NAME       LevelSet_ADV_Filter
    TEST_BASE  LevelSet_ADV_Filter
    PROCS      1
    CASES      0 1 2 3
)
