/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_Library_Interlink_Checks.cpp
 *
 */

#include "fn_Library_Interlink_Checks.hpp"

#include <algorithm>
#include <set>

#include "cl_Library_Enums.hpp"

namespace moris
{
    namespace
    {
        //------------------------------------------------------------------------------------------------------------------

        // splits "A,B;C,D" on aSeparator, dropping empty tokens
        Vector< std::string >
        split_string( const std::string& aString, char aSeparator )
        {
            Vector< std::string > tTokens;
            size_t                tStart = 0;
            while ( tStart <= aString.size() )
            {
                size_t tEnd = aString.find( aSeparator, tStart );
                if ( tEnd == std::string::npos )
                {
                    tEnd = aString.size();
                }
                if ( tEnd > tStart )
                {
                    tTokens.push_back( aString.substr( tStart, tEnd - tStart ) );
                }
                tStart = tEnd + 1;
            }
            return tTokens;
        }

        //------------------------------------------------------------------------------------------------------------------

        uint
        levenshtein_distance( const std::string& aFirst, const std::string& aSecond )
        {
            Vector< uint > tPrevious( aSecond.size() + 1 );
            Vector< uint > tCurrent( aSecond.size() + 1 );
            for ( uint iCol = 0; iCol <= aSecond.size(); iCol++ )
            {
                tPrevious( iCol ) = iCol;
            }
            for ( uint iRow = 1; iRow <= aFirst.size(); iRow++ )
            {
                tCurrent( 0 ) = iRow;
                for ( uint iCol = 1; iCol <= aSecond.size(); iCol++ )
                {
                    uint tCost      = ( aFirst[ iRow - 1 ] == aSecond[ iCol - 1 ] ) ? 0 : 1;
                    tCurrent( iCol ) = std::min( { tPrevious( iCol ) + 1, tCurrent( iCol - 1 ) + 1, tPrevious( iCol - 1 ) + tCost } );
                }
                tPrevious = tCurrent;
            }
            return tPrevious( aSecond.size() );
        }

        //------------------------------------------------------------------------------------------------------------------

        // "did you mean" suffix for a near-miss name, or "" if nothing is close
        std::string
        suggest( const std::string& aName, const std::set< std::string >& aCandidates )
        {
            for ( const std::string& iCandidate : aCandidates )
            {
                if ( levenshtein_distance( aName, iCandidate ) <= 2 )
                {
                    return " — did you mean '" + iCandidate + "'?";
                }
            }
            return "";
        }

        //------------------------------------------------------------------------------------------------------------------

        // collects the names of a FEM submodule ("property_name" of PROPERTIES, ...)
        std::set< std::string >
        collect_names(
                const Module_Parameter_Lists& aFem,
                uint                          aSubmodule,
                const std::string&            aNameKey )
        {
            std::set< std::string > tNames;
            const auto&             tSubmodule = aFem( aSubmodule );
            for ( uint iList = 0; iList < tSubmodule.size(); iList++ )
            {
                if ( tSubmodule( iList ).exists( aNameKey ) )
                {
                    tNames.insert( tSubmodule( iList ).get< std::string >( aNameKey ) );
                }
            }
            return tNames;
        }

        //------------------------------------------------------------------------------------------------------------------

        // checks a "Name,Role;Name,Role" pair-list parameter against defined names
        void
        check_pair_list_references(
                const Parameter_List&           aList,
                const std::string&              aKey,
                const std::string&              aReferrer,
                const std::string&              aKind,
                const std::set< std::string >&  aDefined,
                Vector< std::string >&          aFindings )
        {
            if ( !aList.exists( aKey ) )
            {
                return;
            }
            for ( const std::string& iPair : split_string( aList.get< std::string >( aKey ), ';' ) )
            {
                Vector< std::string > tParts = split_string( iPair, ',' );
                if ( tParts.size() == 0 )
                {
                    continue;
                }
                const std::string& tName = tParts( 0 );
                if ( aDefined.count( tName ) == 0 )
                {
                    aFindings.push_back( "FEM: " + aReferrer + " references unknown " + aKind + " '" + tName + "'" + suggest( tName, aDefined ) );
                }
            }
        }

        //------------------------------------------------------------------------------------------------------------------

    }    // namespace

    //------------------------------------------------------------------------------------------------------------------

    Vector< std::string >
    collect_deck_interlink_findings( const Vector< Module_Parameter_Lists >& aParameterLists )
    {
        Vector< std::string > tFindings;

        const Module_Parameter_Lists& tFem = aParameterLists( (uint)Module_Type::FEM );
        const Module_Parameter_Lists& tGen = aParameterLists( (uint)Module_Type::GEN );
        const Module_Parameter_Lists& tVis = aParameterLists( (uint)Module_Type::VIS );
        const Module_Parameter_Lists& tXtk = aParameterLists( (uint)Module_Type::XTK );

        // a deck without a FEM block (cleared/disabled) has nothing to check here
        if ( tFem.size() == 0 )
        {
            return tFindings;
        }

        // FEM submodule indices (stable under hack_for_legacy_fem, which only erases index 8)
        constexpr uint tPropertiesIndex = 0;
        constexpr uint tCmIndex         = 1;
        constexpr uint tSpIndex         = 2;
        constexpr uint tIwgIndex        = 3;
        constexpr uint tIqiIndex        = 4;

        std::set< std::string > tProperties = collect_names( tFem, tPropertiesIndex, "property_name" );
        std::set< std::string > tCms        = collect_names( tFem, tCmIndex, "constitutive_name" );
        std::set< std::string > tSps        = collect_names( tFem, tSpIndex, "stabilization_name" );
        std::set< std::string > tIqis       = collect_names( tFem, tIqiIndex, "IQI_name" );

        // ---- C4: FEM name references resolve -------------------------------------
        bool tGhostSetsReferenced = false;

        for ( uint iSubmodule : { tIwgIndex, tIqiIndex } )
        {
            const auto& tSubmodule = tFem( iSubmodule );
            for ( uint iList = 0; iList < tSubmodule.size(); iList++ )
            {
                const Parameter_List& tList = tSubmodule( iList );

                std::string tReferrerName = "?";
                std::string tNameKey      = ( iSubmodule == tIwgIndex ) ? "IWG_name" : "IQI_name";
                if ( tList.exists( tNameKey ) )
                {
                    tReferrerName = tList.get< std::string >( tNameKey );
                }

                check_pair_list_references( tList, "leader_properties", tReferrerName, "property", tProperties, tFindings );
                check_pair_list_references( tList, "follower_properties", tReferrerName, "property", tProperties, tFindings );
                check_pair_list_references( tList, "leader_constitutive_models", tReferrerName, "constitutive model", tCms, tFindings );
                check_pair_list_references( tList, "follower_constitutive_models", tReferrerName, "constitutive model", tCms, tFindings );
                check_pair_list_references( tList, "stabilization_parameters", tReferrerName, "stabilization parameter", tSps, tFindings );

                if ( tList.exists( "mesh_set_names" )
                        && tList.get< std::string >( "mesh_set_names" ).find( "ghost_p" ) != std::string::npos )
                {
                    tGhostSetsReferenced = true;
                }
            }
        }

        // CMs reference properties too
        const auto& tCmSubmodule = tFem( tCmIndex );
        for ( uint iList = 0; iList < tCmSubmodule.size(); iList++ )
        {
            const Parameter_List& tList         = tCmSubmodule( iList );
            std::string           tReferrerName = tList.exists( "constitutive_name" ) ? tList.get< std::string >( "constitutive_name" ) : "?";
            check_pair_list_references( tList, "properties", tReferrerName, "property", tProperties, tFindings );
        }

        // SPs reference properties
        const auto& tSpSubmodule = tFem( tSpIndex );
        for ( uint iList = 0; iList < tSpSubmodule.size(); iList++ )
        {
            const Parameter_List& tList         = tSpSubmodule( iList );
            std::string           tReferrerName = tList.exists( "stabilization_name" ) ? tList.get< std::string >( "stabilization_name" ) : "?";
            check_pair_list_references( tList, "leader_properties", tReferrerName, "property", tProperties, tFindings );
            check_pair_list_references( tList, "follower_properties", tReferrerName, "property", tProperties, tFindings );
        }

        // ---- C1 (producer side): GEN IQI_types name defined FEM IQIs -------------
        if ( tGen.size() > 0 && tGen( 0 ).size() > 0 && tGen( 0 )( 0 ).exists( "IQI_types" ) )
        {
            for ( const std::string& iIqi : tGen( 0 )( 0 ).get< Vector< std::string > >( "IQI_types" ) )
            {
                if ( tIqis.count( iIqi ) == 0 )
                {
                    tFindings.push_back( "GEN: IQI_types references unknown FEM IQI '" + iIqi + "'" + suggest( iIqi, tIqis ) );
                }
            }
        }

        // ---- C8: VIS references and list-length agreement ------------------------
        if ( tVis.size() > 0 && tVis( 0 ).size() > 0 )
        {
            const Parameter_List& tVisList = tVis( 0 )( 0 );

            Vector< std::string > tVisIqis;
            if ( tVisList.exists( "IQI_Names" ) )
            {
                tVisIqis = split_string( tVisList.get< std::string >( "IQI_Names" ), ',' );
                for ( const std::string& iIqi : tVisIqis )
                {
                    if ( tIqis.count( iIqi ) == 0 )
                    {
                        tFindings.push_back( "VIS: IQI_Names references unknown FEM IQI '" + iIqi + "'" + suggest( iIqi, tIqis ) );
                    }
                }
            }

            if ( tVisList.exists( "Field_Names" ) && tVisList.exists( "Field_Type" ) )
            {
                uint tNumFields = split_string( tVisList.get< std::string >( "Field_Names" ), ',' ).size();
                uint tNumTypes  = split_string( tVisList.get< std::string >( "Field_Type" ), ',' ).size();

                if ( tNumFields != tNumTypes || ( tVisIqis.size() > 0 && tNumFields != tVisIqis.size() ) )
                {
                    tFindings.push_back( "VIS: Field_Names (" + std::to_string( tNumFields )
                            + "), Field_Type (" + std::to_string( tNumTypes )
                            + "), and IQI_Names (" + std::to_string( tVisIqis.size() )
                            + ") must have matching lengths" );
                }
            }
        }

        // ---- C5 slice: ghost sets require XTK ghost stabilization ----------------
        if ( tGhostSetsReferenced )
        {
            bool tGhostEnabled = tXtk.size() > 0 && tXtk( 0 ).size() > 0
                              && tXtk( 0 )( 0 ).exists( "ghost_stab" )
                              && tXtk( 0 )( 0 ).get< bool >( "ghost_stab" );
            if ( !tGhostEnabled )
            {
                tFindings.push_back( "XTK: FEM references ghost_p* mesh sets but XTK ghost_stab is off — the ghost sets will not exist" );
            }
        }

        return tFindings;
    }

    //------------------------------------------------------------------------------------------------------------------

}    // namespace moris
