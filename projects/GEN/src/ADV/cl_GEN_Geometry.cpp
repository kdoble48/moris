/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Geometry.cpp
 *
 */

#include "cl_GEN_Geometry.hpp"

#include <utility>
#include "cl_GEN_Derived_Node.hpp"

namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------

    Geometry::Geometry( Design_Parameters aParameters )
            : Design( std::move( aParameters ) )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void Geometry::set_node_manager( Node_Manager& aNodeManager )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Geometry::can_change_region_inside_cell() const
    {
        return false;
    }

    //--------------------------------------------------------------------------------------------------------------

    Geometric_Region Geometry::get_geometric_region_at_point( const Matrix< DDRMat >& aCoordinates )
    {
        MORIS_ERROR( false, "Geometry::get_geometric_region_at_point() - not supported by this geometry." );
        return Geometric_Region::INTERFACE;
    }

    //--------------------------------------------------------------------------------------------------------------

}
