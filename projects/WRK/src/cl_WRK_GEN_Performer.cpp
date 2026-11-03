/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_WRK_GEN_Performer.cpp
 *
 */

#include "cl_WRK_GEN_Performer.hpp"

#include <utility>
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_Matrix.hpp"

namespace moris::wrk
{

    // ----------------------------------------------------------------------------

    Gen_Performer::Gen_Performer( std::shared_ptr< moris::gen::Geometry_Engine > aGeometryEngine )
            : mGeometryEngine( std::move( aGeometryEngine ) )
    {
    }

    // ----------------------------------------------------------------------------

    uint
    Gen_Performer::get_num_refinement_fields()
    {
        return mGeometryEngine->get_num_refinement_fields();
    }

    // ----------------------------------------------------------------------------

    gen::Geometric_Region Gen_Performer::get_geometric_region(
            uint                    aFieldIndex,
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        return mGeometryEngine->get_geometric_region( aFieldIndex, aNodeIndex, aCoordinates );
    }

    // ----------------------------------------------------------------------------

    const Vector< uint >&
    Gen_Performer::get_num_refinements( uint aFieldIndex )
    {
        return mGeometryEngine->get_num_refinements( aFieldIndex );
    }

    // ----------------------------------------------------------------------------

    const Vector< uint >&
    Gen_Performer::get_refinement_mesh_indices( uint aFieldIndex )
    {
        return mGeometryEngine->get_refinement_mesh_indices( aFieldIndex );
    }

    // ----------------------------------------------------------------------------

    sint
    Gen_Performer::get_refinement_function_index(
            uint aFieldIndex,
            uint aRefinementIndex )
    {
        return mGeometryEngine->get_refinement_function_index( aFieldIndex, aRefinementIndex );
    }

    // ----------------------------------------------------------------------------

    uint
    Gen_Performer::get_feature_refinement_samples( uint aFieldIndex )
    {
        return mGeometryEngine->get_feature_refinement_samples( aFieldIndex );
    }

    // ----------------------------------------------------------------------------

    bool
    Gen_Performer::can_change_region_inside_cell( uint aFieldIndex )
    {
        return mGeometryEngine->can_change_region_inside_cell( aFieldIndex );
    }

    // ----------------------------------------------------------------------------

    uint
    Gen_Performer::detect_cell_features(
            uint                    aFieldIndex,
            const Matrix< DDRMat >& aCorners,
            uint                    aLagrangeOrder,
            uint                    aSamples )
    {
        return gen::cell_feature_flags( mGeometryEngine->detect_cell_features( aFieldIndex, aCorners, aLagrangeOrder, aSamples ) );
    }

    // ----------------------------------------------------------------------------

    uint
    Gen_Performer::get_num_feature_query_refine_points( uint aFieldIndex )
    {
        return mGeometryEngine->get_feature_query_refine( aFieldIndex ) ? mGeometryEngine->get_num_feature_query_points( aFieldIndex ) : 0;
    }

    // ----------------------------------------------------------------------------

    void
    Gen_Performer::get_feature_query_point(
            uint              aFieldIndex,
            uint              aPointIndex,
            Matrix< DDRMat >& aCoordinates,
            real&             aRadius )
    {
        const gen::Feature_Query_Point& tPoint = mGeometryEngine->get_feature_query_point( aFieldIndex, aPointIndex );

        aCoordinates = tPoint.mCoordinates;
        aRadius      = tPoint.mRadius;
    }

    // ----------------------------------------------------------------------------

    void
    Gen_Performer::get_feature_query_regions(
            uint                    aFieldIndex,
            uint                    aPointIndex,
            const Matrix< DDRMat >& aCorners,
            sint&                   aExpected,
            sint&                   aRepresented )
    {
        aExpected    = mGeometryEngine->get_feature_query_expected_region( aFieldIndex, aPointIndex );
        aRepresented = mGeometryEngine->get_feature_query_represented_region( aFieldIndex, aPointIndex, aCorners );
    }

    // ----------------------------------------------------------------------------
}
