/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * EXA_Common.hpp
 *
 * Shared, physics-free boilerplate callbacks for EXA input decks.
 *
 * MORIS resolves the callbacks named in a deck's parameter lists (e.g. a property's
 * "value_function") at runtime via dlsym() using their exact, UNMANGLED names. They
 * are never called from C++ — only referenced by string — so each must have C linkage
 * (for the unmangled symbol) and external, non-inline linkage (so the compiler emits
 * the symbol even though nothing in the translation unit calls it). Each input deck is
 * compiled to its own .so (a single translation unit), so the definitions below yield
 * exactly one symbol per .so when this header is included once.
 *
 * Usage in a deck:
 *   // (optional) opt into the all-ones constraint-type helper, supplying the count:
 *   #define EXA_NUM_CONSTRAINTS 1
 *   // (optional) opt OUT of the default Output_Criterion if the deck needs its own:
 *   #define EXA_CUSTOM_OUTPUT_CRITERION
 *   #include "EXA_Common.hpp"
 * Then delete the deck's local Func_Const / Output_Criterion / get_constraint_types
 * that the header now provides, or the .so link will fail on a duplicate symbol.
 */

#ifndef MORIS_EXA_COMMON_HPP
#define MORIS_EXA_COMMON_HPP

#include "moris_typedefs.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_TSA_Time_Solver.hpp"

#ifdef __cplusplus
extern "C" {
#endif

    //------------------------------------------------------------------------------
    // Property value function: returns the (single) constant parameter unchanged.
    // Byte-for-byte identical across every EXA deck that uses it. Always provided.
    void
    Func_Const(
            moris::Matrix< moris::DDRMat >&                   aPropMatrix,
            moris::Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*          aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    //------------------------------------------------------------------------------
    // Default time-solver output criterion: always output. Opt out with
    // #define EXA_CUSTOM_OUTPUT_CRITERION when a deck needs different behavior.
#ifndef EXA_CUSTOM_OUTPUT_CRITERION
    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }
#endif

    //------------------------------------------------------------------------------
    // All-ones constraint-type vector of length EXA_NUM_CONSTRAINTS. Opt in by
    // #define EXA_NUM_CONSTRAINTS <count> before including; the size must be supplied
    // because it differs per deck. Decks whose constraint types are not all ones (or
    // who prefer their own) simply leave EXA_NUM_CONSTRAINTS undefined and keep theirs.
#ifdef EXA_NUM_CONSTRAINTS
    moris::Matrix< moris::DDSMat >
    get_constraint_types()
    {
        return moris::Matrix< moris::DDSMat >( EXA_NUM_CONSTRAINTS, 1, 1 );
    }
#endif

#ifdef __cplusplus
}
#endif

#endif    // MORIS_EXA_COMMON_HPP
