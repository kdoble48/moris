/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_Library_Builtin_Functions.cpp
 *
 */

#include "fn_Library_Builtin_Functions.hpp"

#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Vector.hpp"

namespace moris
{
    namespace builtin
    {
        //------------------------------------------------------------------------------------------------------------------

        /**
         * Built-in constant-property callback: copies the first parameter group into the
         * property matrix. ABI-compatible with the FEM property function signature
         * void( Matrix< DDRMat >&, Vector< Matrix< DDRMat > >&, fem::Field_Interpolator_Manager* );
         * the field-interpolator-manager argument is opaque here (and unused) because
         * MRS/IOS cannot depend on FEM headers.
         */
        static void
        Func_Const(
                Matrix< DDRMat >&           aPropMatrix,
                Vector< Matrix< DDRMat > >& aParameters,
                void* /* aFIManager */ )
        {
            aPropMatrix = aParameters( 0 );
        }

        //------------------------------------------------------------------------------------------------------------------

        /**
         * Built-in time-solver output criterion: always write output. ABI-compatible
         * with the SOL warehouse Pointer_Function bool( tsa::Time_Solver* ); the
         * time-solver argument is opaque here (and unused).
         */
        static bool
        Output_Criterion( void* /* aTimeSolver */ )
        {
            return true;
        }

        //------------------------------------------------------------------------------------------------------------------

    }    // namespace builtin

    //------------------------------------------------------------------------------------------------------------------

    void*
    get_builtin_deck_function( const std::string& aName )
    {
        if ( aName == "Func_Const" )
        {
            return reinterpret_cast< void* >( &builtin::Func_Const );
        }
        if ( aName == "Output_Criterion" )
        {
            return reinterpret_cast< void* >( &builtin::Output_Criterion );
        }
        return nullptr;
    }

    //------------------------------------------------------------------------------------------------------------------

}    // namespace moris
