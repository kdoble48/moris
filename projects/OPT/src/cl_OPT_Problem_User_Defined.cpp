/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Problem_User_Defined.cpp
 *
 */

#include "cl_OPT_Problem_User_Defined.hpp"
#include "cl_Library_Factory.hpp"
#include "cl_OPT_Criteria_Interface.hpp"
#include "fn_Parsing_Tools.hpp"

namespace moris::opt
{

    //--------------------------------------------------------------------------------------------------------------

    Problem_User_Defined::Problem_User_Defined(
            Parameter_List&                        aParameterList,
            std::shared_ptr< Criteria_Interface >& aInterface )
            : Problem( aParameterList, aInterface )
    {
        // Load library
        moris::Library_Factory        tLibraryFactory;
        std::shared_ptr< Library_IO > tLibrary     = tLibraryFactory.create_Library( Library_Type::STANDARD );
        std::string                   tLibraryName = aParameterList.get< std::string >( "library" );
        tLibrary->load_parameter_list( tLibraryName, File_Type::SO_FILE );
        tLibrary->finalize();

        // Set user-defined functions. The objective and constraint evaluations are
        // required; the gradient callbacks and get_constraint_types are optional and
        // receive built-in defaults (see install_default_functions) when absent.
        get_constraint_types_user_defined          = tLibrary->load_function< Constraint_Types_Function >( "get_constraint_types", false );
        compute_objectives_user_defined            = tLibrary->load_function< Objective_Constraint_Function >( "compute_objectives" );
        compute_constraints_user_defined           = tLibrary->load_function< Objective_Constraint_Function >( "compute_constraints" );
        compute_dobjective_dadv_user_defined       = tLibrary->load_function< Objective_Constraint_Function >( "compute_dobjective_dadv", false );
        compute_dobjective_dcriteria_user_defined  = tLibrary->load_function< Objective_Constraint_Function >( "compute_dobjective_dcriteria", false );
        compute_dconstraint_dadv_user_defined      = tLibrary->load_function< Objective_Constraint_Function >( "compute_dconstraint_dadv", false );
        compute_dconstraint_dcriteria_user_defined = tLibrary->load_function< Objective_Constraint_Function >( "compute_dconstraint_dcriteria", false );

        this->install_default_functions( aParameterList );
    }

    //--------------------------------------------------------------------------------------------------------------

    Problem_User_Defined::Problem_User_Defined(
            Parameter_List                        aParameterList,
            std::shared_ptr< Criteria_Interface > aInterface,
            Constraint_Types_Function             aConstraintTypesFunction,
            Objective_Constraint_Function         aObjectiveFunction,
            Objective_Constraint_Function         aConstraintFunction,
            Objective_Constraint_Function         aObjectiveADVGradientFunction,
            Objective_Constraint_Function         aObjectiveCriteriaGradientFunction,
            Objective_Constraint_Function         aConstraintADVGradientFunction,
            Objective_Constraint_Function         aConstraintCriteriaGradientFunction )
            : Problem( aParameterList, aInterface )
            , get_constraint_types_user_defined( aConstraintTypesFunction )
            , compute_objectives_user_defined( aObjectiveFunction )
            , compute_constraints_user_defined( aConstraintFunction )
            , compute_dobjective_dadv_user_defined( aObjectiveADVGradientFunction )
            , compute_dobjective_dcriteria_user_defined( aObjectiveCriteriaGradientFunction )
            , compute_dconstraint_dadv_user_defined( aConstraintADVGradientFunction )
            , compute_dconstraint_dcriteria_user_defined( aConstraintCriteriaGradientFunction )
    {
        this->install_default_functions( aParameterList );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Problem_User_Defined::install_default_functions( const Parameter_List& aParameterList )
    {
        // constraint types: deck-defined function > 'constraint_types' parameter > error
        if ( !get_constraint_types_user_defined )
        {
            std::string tConstraintTypesString;
            if ( aParameterList.exists( "constraint_types" ) )
            {
                tConstraintTypesString = aParameterList.get< std::string >( "constraint_types" );
            }

            MORIS_ERROR( !tConstraintTypesString.empty(),
                    "Problem_User_Defined - the input deck defines neither a get_constraint_types() "
                    "function nor the OPT 'constraint_types' parameter. Provide one of the two, e.g. "
                    "constraint_types = \"1\" for a single inequality constraint." );

            Matrix< DDSMat > tConstraintTypes;
            string_to_matrix( tConstraintTypesString, tConstraintTypes );

            MORIS_LOG( "OPT: no get_constraint_types() defined in the deck; using the 'constraint_types' parameter (%u constraint(s)).",
                    (unsigned int)tConstraintTypes.numel() );

            get_constraint_types_user_defined = [ tConstraintTypes ]() {
                return tConstraintTypes;
            };
        }

        // explicit ADV gradients: default to zeros (the explicit ADV dependence of the
        // objective/constraints is zero in the standard criteria-driven formulation;
        // criteria sensitivities are unaffected)
        if ( !compute_dobjective_dadv_user_defined )
        {
            MORIS_LOG( "OPT: no compute_dobjective_dadv() defined in the deck; defaulting to zeros (no explicit ADV dependence of the objective)." );

            compute_dobjective_dadv_user_defined = []( const Vector< real >& aADVs, const Vector< real >& ) {
                return Matrix< DDRMat >( 1, aADVs.size(), 0.0 );
            };
        }

        if ( !compute_dconstraint_dadv_user_defined )
        {
            MORIS_LOG( "OPT: no compute_dconstraint_dadv() defined in the deck; defaulting to zeros (no explicit ADV dependence of the constraints)." );

            compute_dconstraint_dadv_user_defined = [ this ]( const Vector< real >& aADVs, const Vector< real >& ) {
                return Matrix< DDRMat >( this->get_constraint_types_user_defined().numel(), aADVs.size(), 0.0 );
            };
        }

        // criteria gradients encode the problem-specific weights and cannot be defaulted;
        // error at first use so that they can be omitted when never exercised
        if ( !compute_dobjective_dcriteria_user_defined )
        {
            compute_dobjective_dcriteria_user_defined = []( const Vector< real >&, const Vector< real >& ) -> Matrix< DDRMat > {
                MORIS_ERROR( false,
                        "Problem_User_Defined - the input deck does not define compute_dobjective_dcriteria(), "
                        "which is required for gradient-based optimization." );
                return {};
            };
        }

        if ( !compute_dconstraint_dcriteria_user_defined )
        {
            compute_dconstraint_dcriteria_user_defined = []( const Vector< real >&, const Vector< real >& ) -> Matrix< DDRMat > {
                MORIS_ERROR( false,
                        "Problem_User_Defined - the input deck does not define compute_dconstraint_dcriteria(), "
                        "which is required for gradient-based optimization." );
                return {};
            };
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat >
    Problem_User_Defined::get_constraint_types()
    {
        return this->get_constraint_types_user_defined();
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Problem_User_Defined::compute_objectives()
    {
        return this->compute_objectives_user_defined( mADVs, mCriteria );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Problem_User_Defined::compute_constraints()
    {
        return this->compute_constraints_user_defined( mADVs, mCriteria );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Problem_User_Defined::compute_dobjective_dadv()
    {
        return this->compute_dobjective_dadv_user_defined( mADVs, mCriteria );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Problem_User_Defined::compute_dobjective_dcriteria()
    {
        return this->compute_dobjective_dcriteria_user_defined( mADVs, mCriteria );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Problem_User_Defined::compute_dconstraint_dadv()
    {
        return this->compute_dconstraint_dadv_user_defined( mADVs, mCriteria );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Problem_User_Defined::compute_dconstraint_dcriteria()
    {
        return this->compute_dconstraint_dcriteria_user_defined( mADVs, mCriteria );
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::opt
