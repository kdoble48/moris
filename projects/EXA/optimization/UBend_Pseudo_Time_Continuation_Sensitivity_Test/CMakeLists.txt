#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       UBend_Pseudo_Time_Continuation_Sensitivity_Test
    TEST_BASE  uBend_2D_Optimization    # ctest names stay uBend_2D_Optimization-{1,4}-procs
    PROCS      1 4
    # Legacy leaf added GEN/GEN_CORE/src to SO_INCLUDES (kept for fidelity; the
    # directory no longer exists in the tree - GEN sources moved to GEN/src).
    EXTRA_SO_INCLUDES ${MORIS_PACKAGE_DIR}/GEN/GEN_CORE/src
)
