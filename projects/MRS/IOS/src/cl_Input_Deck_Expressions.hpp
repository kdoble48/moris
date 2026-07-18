/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Input_Deck_Expressions.hpp
 *
 * Criteria expressions for the single-entry-point deck API (Deck API v2, see
 * doc/internal/DECK_API_RFC.md): an objective or constraint is written as a small
 * arithmetic expression over named criteria,
 *
 *     aDeck.objective( criterion( "IQIBulkStrainEnergy" ) / tSE0
 *                    + 0.2 * criterion( "IQIPerimeter" ) / tP0 );
 *     aDeck.constraint( criterion( "IQIBulkVolume" ) / tMaxMass - 1.0 <= 0.0 );
 *
 * replacing the seven hand-written OPT callback functions. The criteria index order
 * (GEN IQI_types, coupling contract C1) is assigned by the expression collector in
 * first-appearance order — it cannot be wired wrong by hand — and the criteria
 * gradients are derived from the expression tree by reverse-mode differentiation.
 */

#pragma once

#include <map>
#include <memory>
#include <string>

#include "assert.hpp"
#include "moris_typedefs.hpp"
#include "cl_Vector.hpp"

namespace moris::deck
{
    //------------------------------------------------------------------------------------------------------------------

    /**
     * A node-based scalar expression over named criteria and constants.
     * Cheap to copy (shared immutable nodes); combine with + - * / and unary minus.
     */
    class Expression
    {
      public:
        enum class Op
        {
            LEAF,     // named criterion
            CONST,    // constant
            ADD,
            SUB,
            MUL,
            DIV,
            NEG
        };

        struct Node
        {
            Op                      mOp;
            real                    mValue = 0.0;    // CONST
            std::string             mName;           // LEAF
            std::shared_ptr< Node > mLeft;
            std::shared_ptr< Node > mRight;
        };

      private:
        std::shared_ptr< Node > mRoot;

      public:
        Expression() = default;

        explicit Expression( std::shared_ptr< Node > aRoot )
                : mRoot( std::move( aRoot ) )
        {
        }

        /* implicit */ Expression( real aValue )
                : mRoot( std::make_shared< Node >( Node{ Op::CONST, aValue, "", nullptr, nullptr } ) )
        {
        }

        bool
        is_valid() const
        {
            return mRoot != nullptr;
        }

        const std::shared_ptr< Node >&
        root() const
        {
            return mRoot;
        }

        /**
         * Collects the names of all criteria in this expression in first-appearance
         * (depth-first, left-to-right) order, appending only names not yet present.
         *
         * @param aNames Criteria name list to append to
         */
        void
        collect_criteria( Vector< std::string >& aNames ) const
        {
            collect( mRoot, aNames );
        }

        /**
         * Evaluates the expression for given criteria values.
         *
         * @param aCriteria Criteria values
         * @param aIndexOf Name -> index map into aCriteria
         * @return Expression value
         */
        real
        evaluate(
                const Vector< real >&                  aCriteria,
                const std::map< std::string, uint >& aIndexOf ) const
        {
            return eval( mRoot, aCriteria, aIndexOf );
        }

        /**
         * Accumulates d(expression)/d(criteria) into aGradient by reverse-mode
         * differentiation of the tree.
         *
         * @param aCriteria Criteria values (needed for product/quotient rules)
         * @param aIndexOf Name -> index map into aCriteria
         * @param aGradient Gradient vector (size = number of criteria), accumulated into
         */
        void
        accumulate_gradient(
                const Vector< real >&                  aCriteria,
                const std::map< std::string, uint >& aIndexOf,
                Vector< real >&                        aGradient ) const
        {
            backprop( mRoot, 1.0, aCriteria, aIndexOf, aGradient );
        }

      private:
        static void
        collect( const std::shared_ptr< Node >& aNode, Vector< std::string >& aNames )
        {
            if ( aNode == nullptr )
            {
                return;
            }
            if ( aNode->mOp == Op::LEAF )
            {
                for ( const std::string& iName : aNames )
                {
                    if ( iName == aNode->mName )
                    {
                        return;
                    }
                }
                aNames.push_back( aNode->mName );
                return;
            }
            collect( aNode->mLeft, aNames );
            collect( aNode->mRight, aNames );
        }

        static real
        eval(
                const std::shared_ptr< Node >&        aNode,
                const Vector< real >&                 aCriteria,
                const std::map< std::string, uint >& aIndexOf )
        {
            MORIS_ERROR( aNode != nullptr, "deck::Expression - evaluating an empty expression." );

            switch ( aNode->mOp )
            {
                case Op::CONST:
                    return aNode->mValue;
                case Op::LEAF:
                {
                    auto tIterator = aIndexOf.find( aNode->mName );
                    MORIS_ERROR( tIterator != aIndexOf.end(),
                            "deck::Expression - criterion '%s' is not tracked.", aNode->mName.c_str() );
                    return aCriteria( tIterator->second );
                }
                case Op::ADD:
                    return eval( aNode->mLeft, aCriteria, aIndexOf ) + eval( aNode->mRight, aCriteria, aIndexOf );
                case Op::SUB:
                    return eval( aNode->mLeft, aCriteria, aIndexOf ) - eval( aNode->mRight, aCriteria, aIndexOf );
                case Op::MUL:
                    return eval( aNode->mLeft, aCriteria, aIndexOf ) * eval( aNode->mRight, aCriteria, aIndexOf );
                case Op::DIV:
                    return eval( aNode->mLeft, aCriteria, aIndexOf ) / eval( aNode->mRight, aCriteria, aIndexOf );
                case Op::NEG:
                    return -eval( aNode->mLeft, aCriteria, aIndexOf );
                default:
                    MORIS_ERROR( false, "deck::Expression - unknown operation." );
                    return 0.0;
            }
        }

        static void
        backprop(
                const std::shared_ptr< Node >&        aNode,
                real                                  aSeed,
                const Vector< real >&                 aCriteria,
                const std::map< std::string, uint >& aIndexOf,
                Vector< real >&                       aGradient )
        {
            MORIS_ERROR( aNode != nullptr, "deck::Expression - differentiating an empty expression." );

            switch ( aNode->mOp )
            {
                case Op::CONST:
                    return;
                case Op::LEAF:
                {
                    auto tIterator = aIndexOf.find( aNode->mName );
                    MORIS_ERROR( tIterator != aIndexOf.end(),
                            "deck::Expression - criterion '%s' is not tracked.", aNode->mName.c_str() );
                    aGradient( tIterator->second ) += aSeed;
                    return;
                }
                case Op::ADD:
                    backprop( aNode->mLeft, aSeed, aCriteria, aIndexOf, aGradient );
                    backprop( aNode->mRight, aSeed, aCriteria, aIndexOf, aGradient );
                    return;
                case Op::SUB:
                    backprop( aNode->mLeft, aSeed, aCriteria, aIndexOf, aGradient );
                    backprop( aNode->mRight, -aSeed, aCriteria, aIndexOf, aGradient );
                    return;
                case Op::MUL:
                {
                    real tLeftValue  = eval( aNode->mLeft, aCriteria, aIndexOf );
                    real tRightValue = eval( aNode->mRight, aCriteria, aIndexOf );
                    backprop( aNode->mLeft, aSeed * tRightValue, aCriteria, aIndexOf, aGradient );
                    backprop( aNode->mRight, aSeed * tLeftValue, aCriteria, aIndexOf, aGradient );
                    return;
                }
                case Op::DIV:
                {
                    real tLeftValue  = eval( aNode->mLeft, aCriteria, aIndexOf );
                    real tRightValue = eval( aNode->mRight, aCriteria, aIndexOf );
                    backprop( aNode->mLeft, aSeed / tRightValue, aCriteria, aIndexOf, aGradient );
                    backprop( aNode->mRight, -aSeed * tLeftValue / ( tRightValue * tRightValue ), aCriteria, aIndexOf, aGradient );
                    return;
                }
                case Op::NEG:
                    backprop( aNode->mLeft, -aSeed, aCriteria, aIndexOf, aGradient );
                    return;
                default:
                    MORIS_ERROR( false, "deck::Expression - unknown operation." );
            }
        }
    };

    //------------------------------------------------------------------------------------------------------------------

    /**
     * Creates a leaf expression referencing a criterion (an IQI tracked by GEN) by name.
     */
    inline Expression
    criterion( const std::string& aIqiName )
    {
        return Expression( std::make_shared< Expression::Node >(
                Expression::Node{ Expression::Op::LEAF, 0.0, aIqiName, nullptr, nullptr } ) );
    }

    //------------------------------------------------------------------------------------------------------------------

    inline Expression
    make_binary( Expression::Op aOp, const Expression& aLeft, const Expression& aRight )
    {
        MORIS_ERROR( aLeft.is_valid() && aRight.is_valid(), "deck::Expression - combining an empty expression." );
        return Expression( std::make_shared< Expression::Node >(
                Expression::Node{ aOp, 0.0, "", aLeft.root(), aRight.root() } ) );
    }

    inline Expression operator+( const Expression& aLeft, const Expression& aRight ) { return make_binary( Expression::Op::ADD, aLeft, aRight ); }
    inline Expression operator-( const Expression& aLeft, const Expression& aRight ) { return make_binary( Expression::Op::SUB, aLeft, aRight ); }
    inline Expression operator*( const Expression& aLeft, const Expression& aRight ) { return make_binary( Expression::Op::MUL, aLeft, aRight ); }
    inline Expression operator/( const Expression& aLeft, const Expression& aRight ) { return make_binary( Expression::Op::DIV, aLeft, aRight ); }

    inline Expression
    operator-( const Expression& aOperand )
    {
        MORIS_ERROR( aOperand.is_valid(), "deck::Expression - negating an empty expression." );
        return Expression( std::make_shared< Expression::Node >(
                Expression::Node{ Expression::Op::NEG, 0.0, "", aOperand.root(), nullptr } ) );
    }

    //------------------------------------------------------------------------------------------------------------------

    /**
     * A constraint: an expression plus its type (0 = equality, 1 = inequality),
     * produced by the comparison sugar below. MORIS convention: inequality
     * constraints are satisfied when the value is <= 0.
     */
    struct Constraint
    {
        Expression mExpression;
        sint       mType;    // 0 = equality, 1 = inequality
    };

    inline Constraint
    operator<=( const Expression& aLeft, real aRight )
    {
        return { aRight == 0.0 ? aLeft : aLeft - Expression( aRight ), 1 };
    }

    inline Constraint
    operator==( const Expression& aLeft, real aRight )
    {
        return { aRight == 0.0 ? aLeft : aLeft - Expression( aRight ), 0 };
    }

    //------------------------------------------------------------------------------------------------------------------

}    // namespace moris::deck
