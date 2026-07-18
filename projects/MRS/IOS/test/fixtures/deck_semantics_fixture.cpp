/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * deck_semantics_fixture.cpp
 *
 * Minimal input-deck shared object used by UT_IOS_Library_IO_Deck_Semantics.cpp to pin
 * the Library_IO loading contract. It is compiled as a SHARED library with unresolved
 * moris symbols (like every real input deck) and resolved at dlopen against the test
 * executable's exported symbols.
 *
 * Deliberate structure:
 *   - OPTParameterList: present and MUTATES a default parameter (proves .so applied)
 *   - HMRParameterList: present but EMPTY (pins that defaults are preserved)
 *   - GENParameterList: NOT DEFINED (pins that a missing symbol clears the module)
 *   - Deck_Fixture_User_Function: a user callback resolvable via load_function()
 */

#include "moris_typedefs.hpp"
#include "parameters.hpp"

extern "C" {
namespace moris
{
    //------------------------------------------------------------------------------------------------------------------

    int
    Deck_Fixture_User_Function()
    {
        return 42;
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        // flip a default-seeded parameter so the test can detect that the .so was applied
        aParameterLists.set( "is_optimization_problem", true );
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        // intentionally empty: a present-but-empty function must preserve the
        // fn_PRM-seeded defaults (as opposed to a missing symbol, which clears them)
    }

    //------------------------------------------------------------------------------------------------------------------

    // GENParameterList is intentionally NOT defined: Library_IO clears the GEN module.

    //------------------------------------------------------------------------------------------------------------------

    // Deck-defined version of the well-known "Output_Criterion" callback. Returns the
    // deliberately unusual value FALSE so tests can prove that a deck-defined symbol
    // takes precedence over the built-in implementation (which always returns true).
    // "Func_Const" is intentionally NOT defined so the builtin fallback is exercised.
    bool
    Output_Criterion( void* /* aTimeSolver */ )
    {
        return false;
    }

    //------------------------------------------------------------------------------------------------------------------

}    // namespace moris
}    // extern "C"
