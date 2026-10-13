/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * EXA_Mode.hpp
 *
 * Run mode of an example deck. The EXA examples are regression tests: their
 * optimizations stop after a few iterations so the ctest suite stays fast. The
 * same decks also produce the documentation gallery, which needs the full run.
 *
 *     MORIS_EXA_MODE=showcase   full-length run (documentation gallery)
 *     unset / anything else     regression run (ctest; the default)
 *
 * A deck picks its run length with exa::pick( regression value, showcase value ):
 *
 *     int tMMAMaxIter = exa::pick( 20, 150 );
 *
 * ctest never sets MORIS_EXA_MODE, so the regression tests are unchanged.
 * Header-only: decks are compiled into separate shared objects.
 */

#pragma once

#include <cstdlib>
#include <string>

namespace moris::exa
{
    enum class Mode
    {
        REGRESSION,
        SHOWCASE
    };

    //------------------------------------------------------------------------------

    /** The run mode selected by the MORIS_EXA_MODE environment variable. */
    inline Mode
    mode()
    {
        const char* tValue = std::getenv( "MORIS_EXA_MODE" );
        return ( tValue != nullptr && std::string( tValue ) == "showcase" ) ? Mode::SHOWCASE : Mode::REGRESSION;
    }

    //------------------------------------------------------------------------------

    /** The regression value under ctest, the showcase value for a full-length run. */
    template< typename T >
    T
    pick( T aRegressionValue, T aShowcaseValue )
    {
        return mode() == Mode::SHOWCASE ? aShowcaseValue : aRegressionValue;
    }

}    // namespace moris::exa
