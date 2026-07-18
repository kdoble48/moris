/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * deck_entry_point_fixture.cpp
 *
 * Single-entry-point input deck (MORISInputDeck) used by
 * UT_IOS_Library_IO_Deck_Semantics.cpp to test the Deck API v2 loading path:
 *   - OPT touched and mutated (proves the deck function ran)
 *   - HMR touched but not modified (defaults must be preserved)
 *   - all other modules untouched (must be disabled)
 *   - callbacks registered in-process, including with INTERNAL linkage (no
 *     extern "C", no export) and one shadowing the builtin "Func_Const"
 */

#include "moris_typedefs.hpp"
#include "parameters.hpp"
#include "cl_Input_Deck.hpp"

namespace
{
    // internal linkage on purpose: registered callbacks need no extern "C" export
    int
    Registered_User_Function()
    {
        return 77;
    }

    // registry version of the well-known "Func_Const": doubles the parameter values so
    // tests can distinguish it from the builtin (which copies them unchanged)
    void
    Registered_Func_Const(
            moris::Matrix< moris::DDRMat >&                  aPropMatrix,
            moris::Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            void* /* aFIManager */ )
    {
        aPropMatrix = aParameters( 0 );
        for ( moris::uint iEntry = 0; iEntry < aPropMatrix.numel(); iEntry++ )
        {
            aPropMatrix( iEntry ) *= 2.0;
        }
    }
}    // namespace

MORIS_DECK( aDeck )
{
    // touch OPT and mutate a default parameter
    aDeck.opt().set( "is_optimization_problem", true );

    // touch HMR without modifying anything: defaults must be preserved
    aDeck.hmr();

    // every other module stays untouched and must be disabled

    aDeck.register_function( "Registered_User_Function", &Registered_User_Function );
    aDeck.register_function( "Func_Const", &Registered_Func_Const );
}
