/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * deck_mixed_fixture.cpp
 *
 * Invalid input deck mixing both deck styles: exports MORISInputDeck AND a legacy
 * OPTParameterList symbol. Library_IO must reject it with a hard error.
 */

#include "moris_typedefs.hpp"
#include "parameters.hpp"
#include "cl_Input_Deck.hpp"

MORIS_DECK( aDeck )
{
    aDeck.opt();
}

extern "C" {
namespace moris
{
    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        // intentionally empty
    }
}    // namespace moris
}    // extern "C"
