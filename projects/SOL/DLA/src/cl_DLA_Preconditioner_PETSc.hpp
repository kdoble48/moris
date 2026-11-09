/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Preconditioner_PETSc.hpp
 *
 */

#pragma once

#include <memory>

#include "cl_MatrixPETSc.hpp"

#include "cl_DLA_Preconditioner.hpp"

namespace moris
{
    namespace sol
    {
        class Dist_Map;
        class Dist_Vector;
        class Dist_Matrix;
    }    // namespace sol

    namespace dla
    {
        class Geometric_Multigrid;
        class Linear_Solver_PETSc;
        class Linear_Problem;
        class Preconditioner_PETSc : public Preconditioner
        {
          private:
            PC mpc;

            // block inverse applied by the "mat" preconditioner; owned here, replaced on every build
            sol::Dist_Matrix* mPreconMat = nullptr;

            dla::Geometric_Multigrid* mGeoMultigrid = nullptr;

          protected:

          public:

            /**
             * Constructor
             *
             * @param aParameterlist Parameter list
             */
            Preconditioner_PETSc(
                    const Parameter_List& aParameterlist );

            //-----------------------------------------------------------------------------------

            /**
             * Destructor: releases the block inverse of the "mat" preconditioner
             */
            ~Preconditioner_PETSc() override;

            //-----------------------------------------------------------------------------------

            // the block inverse is owned by this object
            Preconditioner_PETSc( const Preconditioner_PETSc& )            = delete;
            Preconditioner_PETSc& operator=( const Preconditioner_PETSc& ) = delete;

            //-----------------------------------------------------------------------------------

            void build_ilu_preconditioner( Linear_Problem* aLinearSystem );

            //-----------------------------------------------------------------------------------

            void build_multigrid_preconditioner( Linear_Problem* aLinearSystem );

            //-----------------------------------------------------------------------------------

            void build_schwarz_preconditioner_petsc( Linear_Problem* aLinearSystem, KSP aPetscKSPProblem );

            //-----------------------------------------------------------------------------------

            void build_schwarz_preconditioner( Linear_Problem* aLinearSystem );

            //-----------------------------------------------------------------------------------

            void build_algebaric_multigrid_preconditioner( Linear_Problem* aLinearSystem );

            /**
             * @brief compute the preconditioner based on the type of preconditioner
             *
             * @param aLinearSystem linear problem conraining the matrix
             * @param aPetscKSPProblem KSP problem conext used in petsc
             */

            void build_preconditioner( Linear_Problem* aLinearSystem, KSP aPetscKSPProblem );

            //-----------------------------------------------------------------------------------
        };
    }    // namespace dla
}    // namespace moris
