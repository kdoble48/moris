/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * op_move.cpp
 *
 */

#include <catch.hpp>

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_move.hpp"
#include "op_equal_equal.hpp" //ALG
#include "fn_all_true.hpp"
#include "fn_isempty.hpp"
#include "fn_print.hpp"
#include <vector>

namespace moris
{
    TEST_CASE("moris::op_move", "[linalgebra],[op_move]" )
    {
        Matrix< DDRMat > Amatrix(50,50,0.0);
        Matrix< DDRMat > Bmatrix(50,50,0.0);
        Matrix< DDRMat > Cmatrix(50,50);

        Amatrix(0,0)=1.0; Amatrix(0,1)=2.0; Amatrix(0,2)=3.0;
        Amatrix(1,0)=4.0; Amatrix(1,1)=5.0; Amatrix(1,2)=6.0;
        Amatrix(2,0)=7.0; Amatrix(2,1)=8.0; Amatrix(2,2)=9.0;

        Bmatrix = Amatrix.copy();

        Cmatrix = move(Bmatrix);

        REQUIRE( all_true(Cmatrix == Amatrix) );

        REQUIRE(isempty(Bmatrix));
    }

    TEST_CASE("moris::Matrix native move constructor", "[linalgebra],[move_semantics]")
    {
        // Create original matrix with known values
        Matrix<DDRMat> original(10, 10, 1.5);
        original(0, 0) = 42.0;
        original(9, 9) = 99.0;

        // capture the heap buffer address (100 elements > any small-matrix
        // preallocation, so the data lives on the heap for both backends)
        const real* tDataPtr = original.data();

        // Use move constructor
        Matrix<DDRMat> moved(std::move(original));

        // a true move steals the buffer; a copy would allocate a new one
        REQUIRE(moved.data() == tDataPtr);

        // Verify moved matrix has correct dimensions and values
        REQUIRE(moved.n_rows() == 10);
        REQUIRE(moved.n_cols() == 10);
        REQUIRE(moved(0, 0) == 42.0);
        REQUIRE(moved(9, 9) == 99.0);
        REQUIRE(moved(5, 5) == 1.5);
    }

    TEST_CASE("moris::Matrix native move assignment", "[linalgebra],[move_semantics]")
    {
        // Create source matrix
        Matrix<DDRMat> source(5, 5, 2.0);
        source(2, 2) = 7.0;

        // capture the heap buffer address (25 elements > small-matrix preallocation)
        const real* tDataPtr = source.data();

        // Move assign to target
        Matrix<DDRMat> target;
        target = std::move(source);

        // a true move steals the buffer; a copy would allocate a new one
        REQUIRE(target.data() == tDataPtr);

        // Verify target has correct dimensions and values
        REQUIRE(target.n_rows() == 5);
        REQUIRE(target.n_cols() == 5);
        REQUIRE(target(2, 2) == 7.0);
        REQUIRE(target(0, 0) == 2.0);
    }

    TEST_CASE("moris::Matrix vector push_back with move semantics", "[linalgebra],[move_semantics]")
    {
        // This test verifies that std::vector can use move semantics
        // when reallocating, avoiding deep copies

        std::vector<Matrix<DDRMat>> vec;
        // Intentionally NOT reserving to force reallocations

        for (int i = 0; i < 100; ++i)
        {
            Matrix<DDRMat> temp(4, 4, static_cast<real>(i));
            vec.push_back(std::move(temp));
        }

        REQUIRE(vec.size() == 100);

        // Verify data integrity after reallocations
        REQUIRE(vec[0](0, 0) == 0.0);
        REQUIRE(vec[50](0, 0) == 50.0);
        REQUIRE(vec[99](0, 0) == 99.0);
    }
}

