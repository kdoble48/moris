/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Field_Evaluation_Log.cpp
 *
 */

#include "cl_GEN_Field_Evaluation_Log.hpp"

#include <cmath>
#include <iomanip>
#include <limits>

#include "cl_GEN_Derived_Node.hpp"
#include "cl_GEN_Basis_Node.hpp"
#include "cl_GEN_Intersection_Node.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_Logger.hpp"
#include "fn_assert.hpp"

namespace moris::gen
{
    namespace
    {
        // instances per file name in this process: the first truncates the file, later ones append
        std::map< std::string, uint >&
        opened_files()
        {
            static std::map< std::string, uint > tOpened;
            return tOpened;
        }

        Vector< real >
        parse_box( const std::string& aBox )
        {
            Vector< real > tBox;
            std::string    tItem;
            for ( char iChar : aBox + "," )
            {
                if ( iChar == ',' or iChar == ';' )
                {
                    if ( tItem.find_first_not_of( " \t" ) != std::string::npos )
                    {
                        size_t tUsed  = 0;
                        real   tValue = 0.0;
                        try
                        {
                            tValue = std::stod( tItem, &tUsed );
                        } catch ( ... )
                        {
                            tUsed = 0;
                        }
                        MORIS_ERROR( tUsed > 0, "log_field_evaluations_box: cannot read '%s' in '%s'.", tItem.c_str(), aBox.c_str() );
                        tBox.push_back( tValue );
                    }
                    tItem.clear();
                }
                else
                {
                    tItem += iChar;
                }
            }
            MORIS_ERROR( tBox.size() == 0 or tBox.size() == 4 or tBox.size() == 6,
                    "log_field_evaluations_box: '%s' needs 4 (xmin,xmax,ymin,ymax) or 6 (... ,zmin,zmax) numbers.",
                    aBox.c_str() );
            for ( uint iDim = 0; 2 * iDim < tBox.size(); iDim++ )
            {
                MORIS_ERROR( tBox( 2 * iDim ) <= tBox( 2 * iDim + 1 ), "log_field_evaluations_box: min > max in '%s'.", aBox.c_str() );
            }
            return tBox;
        }

        const char*
        tag( const char* aTag )
        {
            return aTag ? aTag : "-";
        }
    }    // namespace

    //--------------------------------------------------------------------------------------------------------------

    std::string
    Field_Evaluation_Log::rank_file_name( const std::string& aFileName )
    {
        std::string tRank  = "_rank" + std::to_string( par_rank() );
        size_t      tSlash = aFileName.find_last_of( '/' );
        size_t      tDot   = aFileName.find_last_of( '.' );
        if ( tDot == std::string::npos or ( tSlash != std::string::npos and tDot < tSlash ) )
        {
            return aFileName + tRank;
        }
        return aFileName.substr( 0, tDot ) + tRank + aFileName.substr( tDot );
    }

    //--------------------------------------------------------------------------------------------------------------

    Field_Evaluation_Log::Field_Evaluation_Log(
            const std::string& aFileName,
            const std::string& aBox,
            uint               aGeometryIndex,
            const std::string& aGeometryName )
            : mFileName( rank_file_name( aFileName ) )
            , mGeometryName( aGeometryName )
            , mGeometryIndex( aGeometryIndex )
            , mBox( parse_box( aBox ) )
    {
        uint& tInstances = opened_files()[ mFileName ];
        mInstance        = tInstances++;
        bool tNewFile    = mInstance == 0 or not std::ifstream( mFileName ).good();
        mFile.open( mFileName, tNewFile ? std::ios::out | std::ios::trunc : std::ios::out | std::ios::app );
        MORIS_ERROR( mFile.good(), "log_field_evaluations: cannot open '%s' for writing.", mFileName.c_str() );
        mFile << std::setprecision( std::numeric_limits< real >::max_digits10 );
        if ( tNewFile )
        {
            mFile << "x,y,z,value,purpose,stage,operation,node_index,pass,geometry,instance\n";
        }
        MORIS_LOG_INFO( "Geometry %u (%s): logging field evaluations to %s%s%s (instance %u)",
                aGeometryIndex,
                aGeometryName.c_str(),
                mFileName.c_str(),
                mBox.size() ? ", box " : "",
                mBox.size() ? aBox.c_str() : "",
                mInstance );
    }

    //--------------------------------------------------------------------------------------------------------------

    Field_Evaluation_Log::~Field_Evaluation_Log()
    {
        try
        {
            this->write_summary();
        } catch ( ... )
        {
            // a destructor must not throw; the rows are already on disk
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    const std::string&
    Field_Evaluation_Log::get_file_name() const
    {
        return mFileName;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Field_Evaluation_Log::write_summary()
    {
        mFile.flush();

        // merge identical tags (equal literals in different translation units can have different addresses)
        std::map< std::tuple< std::string, std::string, std::string, sint >, std::array< luint, 2 > > tMerged;
        for ( const auto& [ iKey, iCount ] : mCounts )
        {
            std::array< luint, 2 >& tCount = tMerged[ { tag( std::get< 0 >( iKey ) ), tag( std::get< 1 >( iKey ) ), tag( std::get< 2 >( iKey ) ), std::get< 3 >( iKey ) } ];
            tCount[ 0 ] += iCount[ 0 ];
            tCount[ 1 ] += iCount[ 1 ];
        }

        std::string   tSummaryName = mFileName + ".summary" + ( mInstance ? "_" + std::to_string( mInstance ) : std::string() ) + ".csv";
        std::ofstream tSummary( tSummaryName, std::ios::out | std::ios::trunc );
        tSummary << "purpose,stage,operation,pass,count_all,count_in_box,geometry,geometry_name,instance\n";
        for ( const auto& [ iKey, iCount ] : tMerged )
        {
            tSummary << std::get< 0 >( iKey ) << ',' << std::get< 1 >( iKey ) << ',' << std::get< 2 >( iKey ) << ','
                     << std::get< 3 >( iKey ) << ',' << iCount[ 0 ] << ',' << iCount[ 1 ] << ',' << mGeometryIndex << ','
                     << mGeometryName << ',' << mInstance << '\n';
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    const char*
    Field_Evaluation_Log::derived_node_kind( const Derived_Node& aDerivedNode )
    {
        if ( dynamic_cast< const Intersection_Node* >( &aDerivedNode ) != nullptr )
        {
            return "intersection_node";
        }

        const Vector< Basis_Node >& tLocators = aDerivedNode.get_locator_nodes();
        uint                        tNumDims  = aDerivedNode.get_global_coordinates().numel();
        uint                        tNum      = tLocators.size();
        bool                        tEqual    = tNum > 1;
        for ( const Basis_Node& iLocator : tLocators )
        {
            tEqual = tEqual and std::abs( iLocator.get_basis() - 1.0 / tNum ) < 1e-12;
        }
        if ( tEqual and tNum == ( 1u << tNumDims ) )
        {
            return "centre_node";
        }
        if ( tEqual and tNumDims == 3 and tNum == 4 )
        {
            return "face_centre_node";
        }
        return "derived_other_node";
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Field_Evaluation_Log::log_background_node( uint aNodeIndex, const Matrix< DDRMat >& aCoordinates, real aValue )
    {
        this->record( "background_node", (sint)aNodeIndex, aCoordinates, aValue );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Field_Evaluation_Log::log_derived_node( const Derived_Node& aDerivedNode, real aValue )
    {
        this->record( derived_node_kind( aDerivedNode ), (sint)aDerivedNode.get_index(), aDerivedNode.get_global_coordinates(), aValue );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Field_Evaluation_Log::log_point( const Matrix< DDRMat >& aCoordinates, real aValue )
    {
        this->record( Field_Evaluation_Context::sPurpose ? Field_Evaluation_Context::sPurpose : "point_untagged", -1, aCoordinates, aValue );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Field_Evaluation_Log::record( const char* aPurpose, sint aNodeIndex, const Matrix< DDRMat >& aCoordinates, real aValue )
    {
        uint tNumCoords = aCoordinates.numel();
        bool tInBox     = true;
        for ( uint iDim = 0; 2 * iDim < mBox.size() and iDim < tNumCoords; iDim++ )
        {
            tInBox = tInBox and aCoordinates( iDim ) >= mBox( 2 * iDim ) and aCoordinates( iDim ) <= mBox( 2 * iDim + 1 );
        }

        std::array< luint, 2 >& tCount = mCounts[ { aPurpose, Field_Evaluation_Context::sStage, Field_Evaluation_Context::sOperation, Field_Evaluation_Context::sPass } ];
        tCount[ 0 ]++;
        if ( not tInBox )
        {
            return;
        }
        tCount[ 1 ]++;

        for ( uint iDim = 0; iDim < 3; iDim++ )
        {
            if ( iDim < tNumCoords )
            {
                mFile << aCoordinates( iDim );
            }
            else
            {
                mFile << 0.0;
            }
            mFile << ',';
        }
        mFile << aValue << ',' << aPurpose << ',' << tag( Field_Evaluation_Context::sStage ) << ','
              << tag( Field_Evaluation_Context::sOperation ) << ',' << aNodeIndex << ',' << Field_Evaluation_Context::sPass << ','
              << mGeometryIndex << ',' << mInstance << '\n';
    }
}    // namespace moris::gen
