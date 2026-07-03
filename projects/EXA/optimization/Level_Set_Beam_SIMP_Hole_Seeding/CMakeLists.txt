#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
# Legacy gate: ctest entries are registered only when GCMMA is available (the deck
# runs the GCMMA algorithm); the input .so and the test TU are built either way,
# exactly as the legacy boilerplate did.
set(EXAMPLE_PROCS 1 4)
if(NOT MORIS_HAVE_GCMMA)
    set(EXAMPLE_PROCS)
endif()

moris_add_example(
    NAME       Level_Set_Beam_SIMP_Hole_Seeding
    TEST_BASE  Level_Set_Beam_SIMP_Hole_Seeding    # ctest names stay Level_Set_Beam_SIMP_Hole_Seeding-{1,4}-procs
    PROCS      ${EXAMPLE_PROCS}
)
