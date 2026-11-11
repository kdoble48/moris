/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Linear_Solver.cpp
 *
 */
#include "cl_DLA_Linear_Solver.hpp"

#include <algorithm>
#include <utility>

#include "cl_DLA_Solver_Factory.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "cl_DLA_Linear_Problem.hpp"
#include "cl_SOL_Enums.hpp"
#include "cl_SOL_Dist_Vector.hpp"

#include "cl_Communication_Tools.hpp"

#include "HDF5_Tools.hpp"

#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

#include "cl_Stopwatch.hpp"    //CHR/src

using namespace moris;
using namespace dla;

//--------------------------------------------------------------------------------------------------

Linear_Solver::Linear_Solver(
        const Parameter_List& aParameterList )
        : mParameterListLinearSolver( aParameterList )
{
}

//--------------------------------------------------------------------------------------------------
void Linear_Solver::set_linear_algorithm( const std::shared_ptr< Linear_Solver_Algorithm >& aLinSolverAlgorithm )
{
    if ( mCallCounter == 0 )
    {
        // removes all elements from the Cell and destroy them
        mLinearSolverList.clear();

        // Resize the Cell to size = 1
        mLinearSolverList.resize( 1 );

        // Set linear solver on first entry
        mLinearSolverList( 0 ) = aLinSolverAlgorithm;
    }
    else
    {
        // set nonlinear solver on next entry
        mLinearSolverList.push_back( aLinSolverAlgorithm );
    }

    mCallCounter = mCallCounter + 1;
}

//-------------------------------------------------------------------------------------------------------
void Linear_Solver::set_linear_algorithm(
        const moris::uint                          aListEntry,
        std::shared_ptr< Linear_Solver_Algorithm > aLinSolverAlgorithm )
{
    // Check if list is smaller than given entry; never shrink it (setting entry 0 after entry 1 dropped entry 1)
    if ( aListEntry >= mLinearSolverList.size() )
    {
        // Resize to new entry value and set nullptr on new entries
        mLinearSolverList.resize( aListEntry + 1, nullptr );
    }
    // Set linear solver on entry
    mLinearSolverList( aListEntry ) = std::move( aLinSolverAlgorithm );
}

//-------------------------------------------------------------------------------------------------------
void Linear_Solver::solver_linear_system(
        dla::Linear_Problem* aLinearProblem,
        const moris::sint    aIter )
{
    Tracer tTracer( "LinearSolver", LOGGER_NON_SPECIFIC_ENTITY_TYPE, "Solve" );

    moris::sint tErrorStatus = 0;

    moris::sint tMaxNumLinRestarts =
            mParameterListLinearSolver.get< moris::sint >( "DLA_max_lin_solver_restarts" );

    std::string tRHSMatrixType = mParameterListLinearSolver.get< std::string >( "RHS_Matrix_Type" );
    aLinearProblem->set_rhs_matrix_type( tRHSMatrixType );

    // if printing of LHS requested through input file, initialize hdf5 files here
    // and save LHS before and after solve
    if ( !this->get_LHS_output_filename().empty() )
    {
        // get current solution vector
        Matrix< DDRMat > tLHS;
        aLinearProblem->get_free_solver_LHS()->extract_copy( tLHS );

        // generate .hdf5 file
        std::string tHdf5FilePath = this->get_LHS_output_filename() + ".hdf5";
        MORIS_LOG_INFO( "Save LHS to file: %s ", tHdf5FilePath.c_str() );
        hid_t  tFileID = create_hdf5_file( tHdf5FilePath );
        herr_t tStatus = 0;

        // write initial guess to hdf5
        save_matrix_to_hdf5_file( tFileID, "init_guess", tLHS, tStatus );

        // solve system
        tErrorStatus = mLinearSolverList( 0 )->solve_linear_system( aLinearProblem, aIter );

        // get update after solve
        aLinearProblem->get_free_solver_LHS()->extract_copy( tLHS );

        // solve update to hdf5
        save_matrix_to_hdf5_file( tFileID, "update", tLHS, tStatus );

        // close hdf5
        close_hdf5_file( tFileID );
    }
    // otherwise, don't save file and simply solve system
    else
    {
        // solve system
        tErrorStatus = mLinearSolverList( 0 )->solve_linear_system( aLinearProblem, aIter );
    }

    // Restart the linear solver using the current solution as an initial guess if the previous linear solve failed.
    // Restart k uses entry k of the solver list, or its last entry once the list is exhausted (MB-100).
    const moris::sint tNumLinSolvers = static_cast< moris::sint >( mLinearSolverList.size() );

    // on request, rebuild the preconditioner for a restart instead of reusing it (iteration index 1 forces a rebuild)
    const moris::sint tRestartIter =
            mParameterListLinearSolver.get< bool >( "DLA_rebuild_lin_solver_on_fail" ) ? 1 : aIter;

    moris::sint tNumRestarts = 0;

    while ( tErrorStatus != 0 && tNumRestarts < tMaxNumLinRestarts )
    {
        tNumRestarts++;

        const moris::sint tSolverIndex = std::min( tNumRestarts, tNumLinSolvers - 1 );

        // Compute current solution vector norm
        Vector< moris::real > tSolVecNorm = aLinearProblem->get_free_solver_LHS()->vec_norm2();

        MORIS_LOG( " ... Previous linear solve failed. Trying restart %i of %i with linear solver algorithm %i, using current solution with SolVecNorm = %5.15e as an initial guess. ",
                tNumRestarts,
                tMaxNumLinRestarts,
                tSolverIndex,
                tSolVecNorm( 0 ) );

        // Re-solve scaled linear system with current solution as an initial guess
        tErrorStatus = mLinearSolverList( tSolverIndex )->solve_linear_system( aLinearProblem, tRestartIter );
    }

    // a failed solve aborts unless the deck asks to continue with the unconverged solution (MB-101)
    if ( tErrorStatus != 0 )
    {
        MORIS_ERROR( !mParameterListLinearSolver.get< bool >( "DLA_hard_break" ),
                "Linear_Solver::solver_linear_system - the linear solve failed with status %i after %i restart(s). "
                "Raise the iteration limit or loosen the tolerance of the linear solver algorithm, allow restarts "
                "(DLA_max_lin_solver_restarts), or set DLA_hard_break = false to continue with the unconverged solution.",
                tErrorStatus,
                tNumRestarts );

        MORIS_LOG_WARNING( "Linear solver did not converge (status %i after %i restart(s)); continuing with the unconverged solution because DLA_hard_break = false.",
                tErrorStatus,
                tNumRestarts );
    }

    if ( !mParameterListLinearSolver.get< std::string >( "DLA_prec_operator_condition_number_with_moris" ).empty() )
    {
        mLinearSolverList( 0 )->compute_preconditioned_operator_condition_number_with_moris( mParameterListLinearSolver.get< std::string >( "DLA_prec_operator_condition_number_with_moris" ) );
    }

    // write out condition numbers of the matrix with and without preconditioner
    if ( !mParameterListLinearSolver.get< std::string >( "DLA_operator_condition_number_with_moris" ).empty() )
    {
        mLinearSolverList( 0 )->compute_operator_condition_number_with_moris( mParameterListLinearSolver.get< std::string >( "DLA_operator_condition_number_with_moris" ) );
    }
}
