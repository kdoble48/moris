/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Cell.cpp
 *
 */

#include <catch.hpp>
#include <iostream>
#include <type_traits>
#include <utility>

// MORIS project header files.
#include "core.hpp"
#include "cl_Vector.hpp"    // CON/src

using namespace moris;

// ----------------------------------------------------------------------------

TEST_CASE( "moris::Vector" )
{

    SECTION( "Vector operator=" )
    {
        Vector< moris::uint > myval   = { 4 };
        Vector< moris::uint > myCell1 = { 10, 2, 3 };

        // Test with single value
        Vector< moris::uint > val1;
        val1 = myval;

        REQUIRE( val1( 0 ) == 4 );

        // Test with row Cell
        Vector< moris::uint > Cell1;
        Cell1 = myCell1;

        REQUIRE( Cell1( 0 ) == 10 );
        REQUIRE( Cell1( 1 ) == 2 );
        REQUIRE( Cell1( 2 ) == 3 );
    }

    SECTION( "Vector size constructor" )
    {
        Vector< moris::uint > myCell1( 10, 4 );
        Vector< moris::real > myCell2( 10, 3.0 );

        REQUIRE( myCell1( 5 ) == 4 );
        REQUIRE( myCell2( 5 ) == 3.0 );
    }

    SECTION( "Vector.size()",
            "Vector.empty()" )
    {
        Vector< moris::uint > myCell1 = { 10, 2, 3 };

        Vector< moris::uint > myCell2;

        moris::size_t Cell_size = myCell1.size();

        bool isempty = myCell2.empty();

        REQUIRE( Cell_size == 3 );
        REQUIRE( isempty == true );
    }

    SECTION( "Vector.append()" )
    {
        Vector< moris::uint > myCell1 = { 4, 5, 9, 8 };
        Vector< moris::uint > myCell2 = { 10, 2, 3 };

        myCell1.append( myCell2 );

        REQUIRE( myCell1( 0 ) == 4 );
        REQUIRE( myCell1( 1 ) == 5 );
        REQUIRE( myCell1( 2 ) == 9 );
        REQUIRE( myCell1( 3 ) == 8 );
        REQUIRE( myCell1( 4 ) == 10 );
        REQUIRE( myCell1( 5 ) == 2 );
        REQUIRE( myCell1( 6 ) == 3 );
    }

    SECTION( "Vector print_as_row_vector" )
    {
        Vector< moris::uint >                     myCell1 = { 4, 5, 9, 8 };
        Vector< Vector< moris::uint > >           myCell2 = { { 10, 2, 3 }, { 4, 5, 9, 8 } };
        Vector< Vector< Vector< moris::uint > > > myCell3 = { { { 1, 2, 3 }, { 4, 5, 6 } }, { { 7, 8, 9 } } };
        ;

        print_as_row_vector( myCell1 );
        print_as_row_vector( myCell2 );
        print_as_row_vector( myCell3 );
    }

    SECTION( "Vector move semantics" )
    {
        // move constructor transfers contents and buffer, leaves source empty
        Vector< moris::uint > tSource = { 10, 2, 3 };
        moris::uint const *   tBuffer = tSource.memptr();

        Vector< moris::uint > tMoved( std::move( tSource ) );

        REQUIRE( tMoved.size() == 3 );
        REQUIRE( tMoved( 0 ) == 10 );
        REQUIRE( tMoved( 1 ) == 2 );
        REQUIRE( tMoved( 2 ) == 3 );
        REQUIRE( tMoved.memptr() == tBuffer );
        REQUIRE( tSource.empty() );

        // move assignment transfers contents and buffer, leaves source empty
        Vector< moris::uint > tAssigned = { 7 };
        tAssigned                       = std::move( tMoved );

        REQUIRE( tAssigned.size() == 3 );
        REQUIRE( tAssigned( 0 ) == 10 );
        REQUIRE( tAssigned( 2 ) == 3 );
        REQUIRE( tAssigned.memptr() == tBuffer );
        REQUIRE( tMoved.empty() );

        // noexcept, so that std::vector reallocation moves instead of copying
        REQUIRE( std::is_nothrow_move_constructible< Vector< moris::uint > >::value );
        REQUIRE( std::is_nothrow_move_assignable< Vector< moris::uint > >::value );

        // outer reallocation keeps the inner buffers
        Vector< Vector< moris::uint > > tOuter = { { 1, 2, 3 }, { 4, 5 }, { 6 } };

        Vector< moris::uint const * > tInnerBuffers( tOuter.size() );
        for ( moris::uint iI = 0; iI < tOuter.size(); iI++ )
        {
            tInnerBuffers( iI ) = tOuter( iI ).memptr();
        }

        Vector< moris::uint > const * tOuterBuffer = tOuter.memptr();
        tOuter.reserve( 4 * tOuter.capacity() );
        REQUIRE( tOuter.memptr() != tOuterBuffer );    // outer storage was reallocated

        for ( moris::uint iI = 0; iI < tOuter.size(); iI++ )
        {
            REQUIRE( tOuter( iI ).memptr() == tInnerBuffers( iI ) );
        }

        tOuter.shrink_to_fit();
        REQUIRE( tOuter.capacity() == 3 );

        for ( moris::uint iI = 0; iI < tOuter.size(); iI++ )
        {
            REQUIRE( tOuter( iI ).memptr() == tInnerBuffers( iI ) );
        }
        REQUIRE( tOuter( 0 )( 0 ) == 1 );
        REQUIRE( tOuter( 0 )( 2 ) == 3 );
        REQUIRE( tOuter( 1 )( 1 ) == 5 );
        REQUIRE( tOuter( 2 )( 0 ) == 6 );

        // copy is still deep
        Vector< moris::uint > tOriginal = { 1, 2, 3 };
        Vector< moris::uint > tCopy( tOriginal );
        tCopy( 0 ) = 99;

        Vector< moris::uint > tCopyAssigned;
        tCopyAssigned      = tOriginal;
        tCopyAssigned( 1 ) = 98;

        REQUIRE( tCopy.memptr() != tOriginal.memptr() );
        REQUIRE( tOriginal( 0 ) == 1 );
        REQUIRE( tOriginal( 1 ) == 2 );
        REQUIRE( tOriginal( 2 ) == 3 );
    }
}
