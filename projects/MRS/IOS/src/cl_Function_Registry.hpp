/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Function_Registry.hpp
 *
 */

#pragma once

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <typeindex>
#include <typeinfo>

#include "assert.hpp"

namespace moris
{
    /**
     * In-process registry of user callbacks for the single-entry-point deck API
     * (see doc/internal/DECK_API_RFC.md). A deck's MORISInputDeck function registers
     * its callbacks here by name instead of exporting them as extern "C" symbols;
     * Library_IO::load_function() consults this registry BEFORE dlsym, so registered
     * callbacks take precedence over both deck symbols and builtins.
     *
     * Entries carry the function-pointer type so that a lookup with the wrong
     * signature fails loudly instead of invoking through a mismatched cast.
     */
    class Function_Registry
    {
      private:
        struct Entry
        {
            void*           mPointer;
            std::type_index mType;
        };

        using Pointer_Map = std::map< std::string, Entry >;

        Pointer_Map mEntries;

      public:
        /**
         * Registers a function pointer under a name. Registering the same name twice
         * is an error (deck bug).
         *
         * @tparam Function_Type Function pointer type (deduced)
         * @param aName Name the callback will be referenced by in parameters
         * @param aFunction Function pointer
         */
        template< typename Function_Type >
        void
        register_function( const std::string& aName, Function_Type aFunction )
        {
            static_assert( std::is_pointer_v< Function_Type >,
                    "Function_Registry::register_function() - a raw function pointer is required." );

            MORIS_ERROR( aFunction != nullptr,
                    "Function_Registry::register_function() - null pointer registered for '%s'.",
                    aName.c_str() );

            bool tInserted = mEntries.emplace( aName, Entry{ reinterpret_cast< void* >( aFunction ), std::type_index( typeid( Function_Type ) ) } ).second;

            MORIS_ERROR( tInserted,
                    "Function_Registry::register_function() - a function named '%s' is already registered.",
                    aName.c_str() );
        }

        /**
         * Looks up a callback by name. Returns nullptr if the name is not registered;
         * errors if the name is registered under a different function type.
         *
         * @tparam Function_Type Function pointer type expected by the caller
         * @param aName Callback name
         * @return Function pointer, or nullptr
         */
        template< typename Function_Type >
        Function_Type
        lookup( const std::string& aName ) const
        {
            Pointer_Map::const_iterator tIterator = mEntries.find( aName );
            if ( tIterator == mEntries.end() )
            {
                return nullptr;
            }

            MORIS_ERROR( tIterator->second.mType == std::type_index( typeid( Function_Type ) ),
                    "Function_Registry::lookup() - the callback '%s' is registered with type '%s' "
                    "but was requested as '%s'. The deck registered a function with the wrong signature.",
                    aName.c_str(),
                    tIterator->second.mType.name(),
                    typeid( Function_Type ).name() );

            return reinterpret_cast< Function_Type >( tIterator->second.mPointer );
        }

        /**
         * Registers a type-erased std::function under a name (capturing lambdas
         * allowed — used by the criteria-expression machinery). Same duplicate rules
         * as register_function().
         *
         * @tparam Signature Function signature, e.g. int( real )
         * @param aName Name the callable will be referenced by
         * @param aFunction Callable
         */
        template< typename Signature >
        void
        register_functional( const std::string& aName, std::function< Signature > aFunction )
        {
            MORIS_ERROR( static_cast< bool >( aFunction ),
                    "Function_Registry::register_functional() - empty function registered for '%s'.",
                    aName.c_str() );

            MORIS_ERROR( mEntries.find( aName ) == mEntries.end(),
                    "Function_Registry::register_functional() - a function pointer named '%s' is already registered.",
                    aName.c_str() );

            bool tInserted = mFunctionals
                                     .emplace( aName,
                                             std::pair< std::shared_ptr< Functional_Holder >, std::type_index >(
                                                     std::make_shared< Functional_Holder_Impl< Signature > >( std::move( aFunction ) ),
                                                     std::type_index( typeid( Signature ) ) ) )
                                     .second;

            MORIS_ERROR( tInserted,
                    "Function_Registry::register_functional() - a functional named '%s' is already registered.",
                    aName.c_str() );
        }

        /**
         * Looks up a type-erased std::function by name. Returns an empty function if
         * the name is not registered; errors on a signature mismatch.
         *
         * @tparam Signature Function signature expected by the caller
         * @param aName Callable name
         * @return The callable, or an empty std::function
         */
        template< typename Signature >
        std::function< Signature >
        lookup_functional( const std::string& aName ) const
        {
            Functional_Map::const_iterator tIterator = mFunctionals.find( aName );
            if ( tIterator == mFunctionals.end() )
            {
                return nullptr;
            }

            MORIS_ERROR( tIterator->second.second == std::type_index( typeid( Signature ) ),
                    "Function_Registry::lookup_functional() - the callable '%s' is registered with signature '%s' "
                    "but was requested as '%s'.",
                    aName.c_str(),
                    tIterator->second.second.name(),
                    typeid( Signature ).name() );

            return static_cast< Functional_Holder_Impl< Signature >* >( tIterator->second.first.get() )->mFunction;
        }

        /**
         * @return Number of registered callbacks
         */
        size_t
        size() const
        {
            return mEntries.size() + mFunctionals.size();
        }

        /**
         * Drops all registered callbacks. MUST be called before dlclosing the deck
         * shared object: registered functionals capture objects allocated by deck
         * code, and destroying them after the .so is unmapped jumps into unmapped
         * deleter code.
         */
        void
        clear()
        {
            mEntries.clear();
            mFunctionals.clear();
        }

      private:
        // Type-erased functional storage. Intentionally NOT std::any: several test
        // translation units '#define private public' before including moris headers,
        // and std::any's access-sensitive internals break under that macro. This
        // plain all-public holder has no such internals.
        struct Functional_Holder
        {
            virtual ~Functional_Holder() = default;
        };

        template< typename Signature >
        struct Functional_Holder_Impl : Functional_Holder
        {
            std::function< Signature > mFunction;

            explicit Functional_Holder_Impl( std::function< Signature > aFunction )
                    : mFunction( std::move( aFunction ) )
            {
            }
        };

        using Functional_Map = std::map< std::string, std::pair< std::shared_ptr< Functional_Holder >, std::type_index > >;

        Functional_Map mFunctionals;
    };

}    // namespace moris
