/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Problem_User_Defined.hpp
 *
 */

#pragma once

#include <functional>

#include "cl_OPT_Problem.hpp"
#include "cl_OPT_Criteria_Interface.hpp"
#include "cl_Parameter_List.hpp"
#include "cl_Library_IO.hpp"

namespace moris::opt
{
    // User-defined function types (raw pointers, as resolved from a deck .so)
    typedef Matrix< DDRMat > ( *Objective_Constraint_Function )(
            const Vector< real >&,
            const Vector< real >& );
    typedef Matrix< DDSMat > ( *Constraint_Types_Function )();

    // Stored callback types: std::function so that built-in defaults (lambdas) can be
    // installed for callbacks the deck does not define
    using Objective_Constraint_Functional = std::function< Matrix< DDRMat >( const Vector< real >&, const Vector< real >& ) >;
    using Constraint_Types_Functional     = std::function< Matrix< DDSMat >() >;

    class Problem_User_Defined : public moris::opt::Problem
    {

      public:
        /**
         * Constructor
         *
         * @param aParameterList parameter list for this problem specifying the needed library
         * @param aInterface Interface class written for other module
         * @param aLibrary Optional already-loaded input library. When provided, callbacks
         *        are resolved from it (in-process registry included) instead of
         *        re-opening the .so named by the 'library' parameter.
         */
        Problem_User_Defined(
                Parameter_List&                        aParameterList,
                std::shared_ptr< Criteria_Interface >& aInterface,
                std::shared_ptr< Library_IO >          aLibrary = nullptr );

        /**
         * Alternate constructor where the user-defined functions are provided directly. Used for OPT tests.
         *
         * @param aParameterList Parameter list for the base Problem class
         * @param aInterface Interface class written for other module
         * @param aConstraintTypesFunction Function for getting constraint types
         * @param aObjectiveFunction Objective function
         * @param aConstraintFunction Constraint function
         * @param aObjectiveADVGradientFunction
         * @param aObjectiveCriteriaGradientFunction
         * @param aConstraintADVGradientFunction
         * @param aConstraintCriteriaGradientFunction
         */
        Problem_User_Defined(
                Parameter_List                        aParameterList,
                std::shared_ptr< Criteria_Interface > aInterface,
                Constraint_Types_Function             aConstraintTypesFunction,
                Objective_Constraint_Function         aObjectiveFunction,
                Objective_Constraint_Function         aConstraintFunction,
                Objective_Constraint_Function         aObjectiveADVGradientFunction,
                Objective_Constraint_Function         aObjectiveCriteriaGradientFunction,
                Objective_Constraint_Function         aConstraintADVGradientFunction,
                Objective_Constraint_Function         aConstraintCriteriaGradientFunction );

        /**
         * Gets the constraint types
         *
         * @return vector of integers, 0 = equality constraint, 1 = inequality constraint
         */
        Matrix< DDSMat > get_constraint_types() override;

      protected:
        /**
         * Gets the objective values
         *
         * @return vector of objectives
         */
        Matrix< DDRMat > compute_objectives() override;

        /**
         * Gets the constraint values
         *
         * @return vector of constraints
         */
        Matrix< DDRMat > compute_constraints() override;

        /**
         * Gets the derivative of the objectives with respect to the advs
         *
         * @return matrix d(objective)_i/d(adv)_j
         */
        Matrix< DDRMat > compute_dobjective_dadv() override;

        /**
         * Gets the derivative of the constraints with respect to the advs
         *
         * @return matrix d(constraints)_i/d(adv)_j
         */
        Matrix< DDRMat > compute_dconstraint_dadv() override;

        /**
         * Gets the derivative of the objective with respect to the criteria.
         *
         * @return matrix d(objective)_i/d(criteria)_j
         */
        Matrix< DDRMat > compute_dobjective_dcriteria() override;

        /**
         * Gets the derivative of the constraints with respect to the criteria.
         *
         * @return matrix d(constraint)_i/d(criteria)_j
         */
        Matrix< DDRMat > compute_dconstraint_dcriteria() override;

      private:
        /**
         * Installs built-in defaults for any callback the deck did not define:
         * - get_constraint_types: taken from the OPT 'constraint_types' parameter
         *   (comma-separated 0=equality / 1=inequality per constraint); an error if
         *   neither the function nor the parameter is provided
         * - compute_dobjective_dadv / compute_dconstraint_dadv: zero matrices (the
         *   explicit ADV dependence is assumed zero; criteria sensitivities are
         *   unaffected) — announced with a log message
         * - compute_dobjective_dcriteria / compute_dconstraint_dcriteria: an error
         *   at first use (required for gradient-based optimization)
         *
         * @param aParameterList OPT problem parameter list (source of 'constraint_types')
         */
        void install_default_functions( const Parameter_List& aParameterList );

        Constraint_Types_Functional     get_constraint_types_user_defined;
        Objective_Constraint_Functional compute_objectives_user_defined;
        Objective_Constraint_Functional compute_constraints_user_defined;
        Objective_Constraint_Functional compute_dobjective_dadv_user_defined;
        Objective_Constraint_Functional compute_dobjective_dcriteria_user_defined;
        Objective_Constraint_Functional compute_dconstraint_dadv_user_defined;
        Objective_Constraint_Functional compute_dconstraint_dcriteria_user_defined;
    };
}    // namespace moris::opt
