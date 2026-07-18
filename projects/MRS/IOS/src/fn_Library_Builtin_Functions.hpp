/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_Library_Builtin_Functions.hpp
 *
 */

#pragma once

#include <string>

namespace moris
{
    /**
     * Returns a pointer to a built-in implementation of a well-known input-deck callback,
     * or nullptr if no builtin exists for the given name.
     *
     * Consulted by Library_IO::load_function() only AFTER dlsym on the deck's shared
     * object fails, so a deck-defined symbol always takes precedence over the builtin.
     * This makes the most common copy-paste deck functions (the constant-property
     * callback "Func_Const" and the always-true time-solver output criterion
     * "Output_Criterion") optional: a deck that references them by name no longer
     * needs to define them.
     *
     * @param aName Function name as referenced by a parameter (e.g. "Func_Const")
     * @return Pointer to the builtin implementation, or nullptr
     */
    void* get_builtin_deck_function( const std::string& aName );

}    // namespace moris
