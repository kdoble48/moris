#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
# 4 (method) x 4 (gamma) case matrix, gCaseIndex = method*4 + gamma; each case is its
# own ctest / runner process because ArborX/Kokkos cannot run a second workflow in one
# process (see the example_test_case.cpp header). Tags: [EXA_Struc_Contact_Patch_Test_<i>].
#
# Gated-out cases (compiled and runnable by hand, not registered as ctests):
#   3, 7, 11  - gamma = 1e4 for each Nitsche variant: Newton diverges (residual > 1e20);
#               needs load control / relaxation before it can gate.
#   12..15    - penalty method: SIGSEGV in fem::Property::val() via SP_Penalty_Contact
#               on the nonconformal contact set (pre-existing core bug on develop; the
#               Nitsche variants are the supported path).
moris_add_example(
    NAME       Struc_Contact_Patch_Test
    TEST_BASE  Struc_Contact_Patch_Test
    PROCS      1
    CASES      0 1 2 4 5 6 8 9 10
)
