#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       LevelSet_Boxbeam_Adaptive_Refinement
    TEST_BASE  Levelset_Boxbeam_Adaptive_Refinement    # lowercase 's': ctest names stay Levelset_Boxbeam_Adaptive_Refinement-{1,4}-procs
    PROCS      1 4
    INPUTS     Levelset_Boxbeam_Adaptive_Refinement
               Levelset_Boxbeam_Adaptive_Refinement_Quadratic
    # NOTE: legacy add_test was additionally gated on MORIS_HAVE_GCMMA (default ON);
    #       moris_add_example registers under MORIS_HAVE_PARALLEL_TESTS only.
)
