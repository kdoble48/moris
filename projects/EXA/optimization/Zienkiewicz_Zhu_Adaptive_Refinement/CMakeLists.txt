#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
# Legacy gate: ctest entries are registered only when GCMMA is available (the deck
# runs the GCMMA algorithm); the input .so and the test TU are built either way,
# exactly as the legacy boilerplate did.
set(EXAMPLE_PROCS 1)
if(NOT MORIS_HAVE_GCMMA)
    set(EXAMPLE_PROCS)
endif()

moris_add_example(
    NAME       Zienkiewicz_Zhu_Adaptive_Refinement
    TEST_BASE  Zienkiewicz_Zhu_Adaptive_Refinement    # ctest name stays Zienkiewicz_Zhu_Adaptive_Refinement-1-procs
    PROCS      ${EXAMPLE_PROCS}
)
