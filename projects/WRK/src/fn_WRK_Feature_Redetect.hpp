/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_WRK_Feature_Redetect.hpp
 *
 * MB-089: re-detection of sub-cell features during an optimization. The GEN sub-cell feature refinement (MB-084) is
 * computed for one design; as the design moves, thin walls can form in cells that were not refined for them, and the
 * XTK geometry resolution check (geometry_resolution_samples) counts the background edges that the interface crosses
 * more often than the cut. This decides, from those counts, whether the workflow redoes the GEN refinement on the
 * current design. A redo cannot refine cells already on the last level (number_of_refinements): features thinner than
 * that level's cells stay hidden, and the post-redo count says so.
 *
 */

#pragma once

#include <algorithm>

#include "moris_typedefs.hpp"
#include "cl_Vector.hpp"

namespace moris::wrk
{
    struct Feature_Redetect_Decision
    {
        bool mRedo        = false;    // redo the GEN refinement on the current design
        uint mGeometry    = 0;        // the geometry that triggered it (largest excess over its threshold)
        uint mHiddenEdges = 0;        // its hidden-edge count on the current cut
        uint mThreshold   = 0;        // the threshold it exceeded: max( deck threshold, count on the first cut after the last refinement )
    };

    /**
     * Whether to redo the GEN refinement on the current design (MB-089).
     *
     * @param aThresholds Per geometry: feature_refinement_redetect_edges (-1 = off)
     * @param aHiddenEdges Per geometry: background edges the interface crosses more often than the cut, from the XTK
     *        geometry resolution check of the current cut (shorter or empty when the check is off: no trigger)
     * @param aBaseline Per geometry: the same count on the first cut after the last (re)refinement, i.e. on the mesh
     *        refined for the design it was cut with (empty before that cut). A redo cannot improve on it, so the
     *        trigger only fires above it; this also keeps nearby designs (finite-difference perturbations, inner
     *        trials) from re-firing on the mesh just refined for them.
     * @param aMeshRefinedForThisDesign The current design is the one the mesh was refined for: a redo would rebuild the
     *        same mesh, never fire
     */
    inline Feature_Redetect_Decision
    decide_feature_redetect(
            const Vector< sint >& aThresholds,
            const Vector< uint >& aHiddenEdges,
            const Vector< uint >& aBaseline,
            bool                  aMeshRefinedForThisDesign )
    {
        Feature_Redetect_Decision tDecision;
        if ( aMeshRefinedForThisDesign )
        {
            return tDecision;
        }

        uint tLargestExcess = 0;
        for ( uint iGeom = 0; iGeom < aThresholds.size() and iGeom < aHiddenEdges.size(); iGeom++ )
        {
            if ( aThresholds( iGeom ) < 0 )
            {
                continue;
            }
            uint tThreshold = uint( aThresholds( iGeom ) );
            if ( iGeom < aBaseline.size() )
            {
                tThreshold = std::max( tThreshold, aBaseline( iGeom ) );
            }
            if ( aHiddenEdges( iGeom ) > tThreshold and aHiddenEdges( iGeom ) - tThreshold > tLargestExcess )
            {
                tLargestExcess          = aHiddenEdges( iGeom ) - tThreshold;
                tDecision.mRedo         = true;
                tDecision.mGeometry     = iGeom;
                tDecision.mHiddenEdges  = aHiddenEdges( iGeom );
                tDecision.mThreshold    = tThreshold;
            }
        }
        return tDecision;
    }
}    // namespace moris::wrk
