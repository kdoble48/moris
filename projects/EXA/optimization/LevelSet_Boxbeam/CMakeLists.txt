#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       LevelSet_Boxbeam
    TEST_BASE  Levelset_Boxbeam    # lowercase 's': ctest names stay Levelset_Boxbeam-1-procs / Levelset_Boxbeam-3-procs
    PROCS      1 3
    # INPUTS defaults to TEST_BASE: dynamic_link_input(Levelset_Boxbeam ...) as before.
    # Levelset_Boxbeam_Lin_Quad_LevelSet.cpp is NOT built (legacy reference was commented out).
    # NOTE: legacy add_test was additionally gated on MORIS_HAVE_GCMMA (default ON);
    #       moris_add_example registers under MORIS_HAVE_PARALLEL_TESTS only.
)
