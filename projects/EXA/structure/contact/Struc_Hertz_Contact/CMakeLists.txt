#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
# 3 (penetration) x 3 (mesh) case matrix; each case is its own ctest / runner process
# because ArborX/Kokkos cannot run a second workflow in one process (see the
# example_test_case.cpp header). Tags: [EXA_Struc_Hertz_Contact_<0..8>].
# All 9 cases pass (2026-07-18; case 8 = deepest penetration x 160^2 mesh, ~78 s).
moris_add_example(
    NAME       Struc_Hertz_Contact
    TEST_BASE  Struc_Hertz_Contact
    PROCS      1
    CASES      0 1 2 3 4 5 6 7 8
)
