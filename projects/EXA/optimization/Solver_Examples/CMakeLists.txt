#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       Solver_Examples
    TEST_BASE  Solver_Examples    # ctest name stays exactly Solver_Examples (legacy add_test had no -procs suffix)
    PROCS      1                  # legacy add_test ran mpirun -n 1
    NO_PROCS_SUFFIX
    # INPUTS defaults to TEST_BASE: dynamic_link_input(Solver_Examples ...) as before.
    # One deck, four TEST_CASEs: each runs fn_WRK_Workflow_Main_Interface against the same
    # Solver_Examples.so with different solver-configuration globals (see example_test_case.cpp).
    # NOTE: legacy add_test was additionally gated on MORIS_HAVE_GCMMA (default ON);
    #       moris_add_example registers under MORIS_HAVE_PARALLEL_TESTS only.
)
