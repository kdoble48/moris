/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Input_Deck.hpp
 *
 */

#pragma once

#include <array>

#include "cl_Library_Enums.hpp"
#include "cl_Module_Parameter_Lists.hpp"
#include "cl_Function_Registry.hpp"
#include "cl_Input_Deck_Expressions.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Vector.hpp"

namespace moris
{
    /**
     * The single-entry-point deck API (see doc/internal/DECK_API_RFC.md). Instead of
     * exporting one <MODULE>ParameterList symbol per module, a deck exports exactly one
     * C-linkage symbol:
     *
     *     extern "C" void MORISInputDeck( moris::Input_Deck& aDeck );
     *
     * Inside it, modules are configured through the accessors below (first touch
     * activates a module with its fn_PRM defaults; untouched modules are disabled),
     * and user callbacks are registered directly as function pointers via
     * register_function() — no extern "C", no fixed names, no dlsym string coupling:
     *
     *     aDeck.hmr().set( "number_of_elements_per_dimension", 45u, 15u );
     *     aDeck.register_function( "Func_Neumann", &My_Neumann_Function );
     *
     * The MORIS_DECK convenience macro hides the linkage ceremony:
     *
     *     MORIS_DECK( aDeck ) { ... }
     */
    class Input_Deck
    {
      private:
        Vector< Module_Parameter_Lists >&                        mParameterLists;
        Function_Registry&                                       mRegistry;
        std::array< bool, (size_t)( Module_Type::END_ENUM ) > mTouched{};

        // declared optimization expressions (materialized by finalize_expressions)
        deck::Expression           mObjective;
        Vector< deck::Constraint > mConstraints;
        Vector< std::string >      mCriteriaOrder;

      public:
        /**
         * Constructed by Library_IO over its parameter-list storage and registry;
         * not constructed by user code.
         */
        Input_Deck(
                Vector< Module_Parameter_Lists >& aParameterLists,
                Function_Registry&                aRegistry )
                : mParameterLists( aParameterLists )
                , mRegistry( aRegistry )
        {
        }

        /**
         * Accesses a module's parameter lists. The first touch activates the module
         * (it keeps the fn_PRM-seeded defaults the same way a present-but-empty
         * <MODULE>ParameterList function does in the legacy format); modules never
         * touched are disabled (cleared) after the deck function returns.
         *
         * @param aModule Module to configure
         * @return Module parameter lists, pre-seeded with defaults
         */
        Module_Parameter_Lists&
        module( Module_Type aModule )
        {
            MORIS_ERROR( aModule < Module_Type::END_ENUM,
                    "Input_Deck::module() - invalid module type." );

            mTouched[ (size_t)aModule ] = true;
            return mParameterLists( (uint)aModule );
        }

        // named accessors, mirroring the module structure (and the pymoris vocabulary)
        Module_Parameter_Lists& opt() { return this->module( Module_Type::OPT ); }
        Module_Parameter_Lists& hmr() { return this->module( Module_Type::HMR ); }
        Module_Parameter_Lists& stk() { return this->module( Module_Type::STK ); }
        Module_Parameter_Lists& xtk() { return this->module( Module_Type::XTK ); }
        Module_Parameter_Lists& gen() { return this->module( Module_Type::GEN ); }
        Module_Parameter_Lists& fem() { return this->module( Module_Type::FEM ); }
        Module_Parameter_Lists& sol() { return this->module( Module_Type::SOL ); }
        Module_Parameter_Lists& msi() { return this->module( Module_Type::MSI ); }
        Module_Parameter_Lists& vis() { return this->module( Module_Type::VIS ); }
        Module_Parameter_Lists& mig() { return this->module( Module_Type::MIG ); }
        Module_Parameter_Lists& morisgeneral() { return this->module( Module_Type::MORISGENERAL ); }

        /**
         * Registers a user callback by name; parameters referencing this name
         * (value_function, field_function_name, TSA_Output_Criteria, the OPT
         * callback names, ...) resolve to it with precedence over deck symbols
         * and builtins.
         *
         * @tparam Function_Type Function pointer type (deduced)
         * @param aName Name to register under
         * @param aFunction Function pointer
         */
        template< typename Function_Type >
        void
        register_function( const std::string& aName, Function_Type aFunction )
        {
            mRegistry.register_function( aName, aFunction );
        }

        /**
         * Declares the optimization objective as an expression over named criteria,
         * e.g.  aDeck.objective( deck::criterion( "IQIBulkStrainEnergy" ) / tSE0 );
         *
         * Replaces the compute_objectives / compute_dobjective_dcriteria /
         * compute_dobjective_dadv callbacks; the criteria order (GEN IQI_types) is
         * assigned by the collector in first-appearance order, and the criteria
         * gradient is derived from the expression by reverse-mode differentiation.
         */
        void
        objective( const deck::Expression& aExpression )
        {
            MORIS_ERROR( !mObjective.is_valid(),
                    "Input_Deck::objective() - the objective has already been declared." );
            MORIS_ERROR( aExpression.is_valid(),
                    "Input_Deck::objective() - empty expression." );

            mObjective = aExpression;
            aExpression.collect_criteria( mCriteriaOrder );
        }

        /**
         * Adds an optimization constraint, e.g.
         *   aDeck.constraint( deck::criterion( "IQIBulkVolume" ) / tMaxMass - 1.0 <= 0.0 );
         * Constraint order must list equality constraints before inequality ones
         * (MORIS convention, checked by the OPT problem).
         */
        void
        constraint( const deck::Constraint& aConstraint )
        {
            MORIS_ERROR( aConstraint.mExpression.is_valid(),
                    "Input_Deck::constraint() - empty expression." );

            mConstraints.push_back( aConstraint );
            aConstraint.mExpression.collect_criteria( mCriteriaOrder );
        }

        /**
         * Whether a module was touched by the deck function.
         *
         * @param aModule Module type
         * @return true if the deck accessed this module
         */
        bool
        is_touched( Module_Type aModule ) const
        {
            return mTouched[ (size_t)aModule ];
        }

        /**
         * Marks a module as touched without returning it. Used by the loader to
         * auto-activate OPT (mandatory for workflow selection).
         */
        void
        touch( Module_Type aModule )
        {
            mTouched[ (size_t)aModule ] = true;
        }

        /**
         * Materializes declared objective/constraint expressions. Called by the
         * loader after the deck function returns (before untouched modules are
         * cleared): writes GEN IQI_types in collection order (contract C1), the OPT
         * constraint_types parameter, and registers the four OPT evaluation
         * callbacks (values by tree evaluation, criteria gradients by reverse-mode
         * differentiation) in the function registry. No-op if no objective was
         * declared.
         */
        void
        finalize_expressions()
        {
            if ( !mObjective.is_valid() )
            {
                MORIS_ERROR( mConstraints.empty(),
                        "Input_Deck - constraints were declared but no objective; declare one with aDeck.objective()." );
                return;
            }

            // criteria order -> GEN IQI_types (touches GEN: criteria are GEN-tracked IQIs)
            this->module( Module_Type::GEN )( 0 ).set( "IQI_types", mCriteriaOrder );

            // constraint types -> OPT parameter (consumed by Problem_User_Defined)
            std::string tConstraintTypes;
            for ( uint iConstraint = 0; iConstraint < mConstraints.size(); iConstraint++ )
            {
                tConstraintTypes += ( iConstraint > 0 ? "," : "" ) + std::to_string( mConstraints( iConstraint ).mType );
            }
            this->module( Module_Type::OPT )( 0 ).set( "constraint_types", tConstraintTypes );

            // shared final name -> index map for the evaluation lambdas
            std::shared_ptr< std::map< std::string, uint > > tIndexOf = std::make_shared< std::map< std::string, uint > >();
            for ( uint iCriterion = 0; iCriterion < mCriteriaOrder.size(); iCriterion++ )
            {
                ( *tIndexOf )[ mCriteriaOrder( iCriterion ) ] = iCriterion;
            }

            using Callback = Matrix< DDRMat >( const Vector< real >&, const Vector< real >& );

            deck::Expression           tObjective   = mObjective;
            Vector< deck::Constraint > tConstraints = mConstraints;

            mRegistry.register_functional< Callback >( "compute_objectives",
                    [ tObjective, tIndexOf ]( const Vector< real >&, const Vector< real >& aCriteria ) {
                        return Matrix< DDRMat >{ { tObjective.evaluate( aCriteria, *tIndexOf ) } };
                    } );

            mRegistry.register_functional< Callback >( "compute_dobjective_dcriteria",
                    [ tObjective, tIndexOf ]( const Vector< real >&, const Vector< real >& aCriteria ) {
                        Vector< real > tGradient( aCriteria.size(), 0.0 );
                        tObjective.accumulate_gradient( aCriteria, *tIndexOf, tGradient );
                        Matrix< DDRMat > tResult( 1, aCriteria.size(), 0.0 );
                        for ( uint iCriterion = 0; iCriterion < aCriteria.size(); iCriterion++ )
                        {
                            tResult( 0, iCriterion ) = tGradient( iCriterion );
                        }
                        return tResult;
                    } );

            mRegistry.register_functional< Callback >( "compute_constraints",
                    [ tConstraints, tIndexOf ]( const Vector< real >&, const Vector< real >& aCriteria ) {
                        Matrix< DDRMat > tResult( tConstraints.size(), 1, 0.0 );
                        for ( uint iConstraint = 0; iConstraint < tConstraints.size(); iConstraint++ )
                        {
                            tResult( iConstraint ) = tConstraints( iConstraint ).mExpression.evaluate( aCriteria, *tIndexOf );
                        }
                        return tResult;
                    } );

            mRegistry.register_functional< Callback >( "compute_dconstraint_dcriteria",
                    [ tConstraints, tIndexOf ]( const Vector< real >&, const Vector< real >& aCriteria ) {
                        Matrix< DDRMat > tResult( tConstraints.size(), aCriteria.size(), 0.0 );
                        for ( uint iConstraint = 0; iConstraint < tConstraints.size(); iConstraint++ )
                        {
                            Vector< real > tGradient( aCriteria.size(), 0.0 );
                            tConstraints( iConstraint ).mExpression.accumulate_gradient( aCriteria, *tIndexOf, tGradient );
                            for ( uint iCriterion = 0; iCriterion < aCriteria.size(); iCriterion++ )
                            {
                                tResult( iConstraint, iCriterion ) = tGradient( iCriterion );
                            }
                        }
                        return tResult;
                    } );
        }
    };

    // Deck entry-point function type, resolved by Library_IO under the fixed name
    // "MORISInputDeck". A deck may export EITHER this symbol OR legacy
    // <MODULE>ParameterList symbols — never both.
    typedef void ( *Input_Deck_Function )( Input_Deck& aDeck );

}    // namespace moris

/**
 * Convenience macro hiding the extern "C" ceremony of the deck entry point.
 * Usage:  MORIS_DECK( aDeck ) { aDeck.hmr().set(...); ... }
 */
#define MORIS_DECK( arg ) \
    extern "C" void MORISInputDeck( moris::Input_Deck& arg )
