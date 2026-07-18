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

#include <map>
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

        std::map< std::string, Entry > mEntries;

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
            auto tIterator = mEntries.find( aName );
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
         * @return Number of registered callbacks
         */
        size_t
        size() const
        {
            return mEntries.size();
        }
    };

}    // namespace moris
