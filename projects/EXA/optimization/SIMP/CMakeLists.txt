#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       SIMP
    TEST_BASE  SIMP             # ctest name stays exactly SIMP (legacy add_test had no -procs suffix)
    PROCS      1                # legacy add_test ran mpirun -n 1
    NO_PROCS_SUFFIX
    # INPUTS defaults to TEST_BASE: dynamic_link_input(SIMP SIMP SIMP.cpp ...) as before
    # NOTE: legacy add_test was additionally gated on MORIS_HAVE_GCMMA (default ON);
    #       moris_add_example registers under MORIS_HAVE_PARALLEL_TESTS only.
)
