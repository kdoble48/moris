#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       LevelSet_Boxbeam_Adaptive_Refinement_Restart
    TEST_BASE  Levelset_Boxbeam_Adaptive_Refinement_Restart    # lowercase 's': ctest name stays Levelset_Boxbeam_Adaptive_Refinement_Restart-1-procs
    PROCS      1
    INPUTS     Levelset_Boxbeam_Create_File    # first TEST_CASE writes ADV hdf5 ...
               Levelset_Boxbeam_Restart        # ... second TEST_CASE restarts from it (intra-TU declaration order is load-bearing)
    # NOTE: legacy add_test was additionally gated on MORIS_HAVE_GCMMA (default ON);
    #       moris_add_example registers under MORIS_HAVE_PARALLEL_TESTS only.
)
