/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_WRK_Feature_Redetect.cpp
 *
 * MB-089: the trigger that redoes the GEN sub-cell feature refinement on the current design during an optimization,
 * from the hidden-edge counts of the XTK geometry resolution check (decide_feature_redetect()).
 *
 */

#include "catch.hpp"

#include "fn_WRK_Feature_Redetect.hpp"

namespace moris::wrk
{
    TEST_CASE( "MB-089 feature re-detection trigger", "[WRK],[MB-089]" )
    {
        // geometry 0: the design, threshold 40; geometry 1: fixed pads, re-detection off
        Vector< sint > tThresholds = { 40, -1 };

        SECTION( "off: every threshold -1, or the check off" )
        {
            CHECK( not decide_feature_redetect( { -1, -1 }, { 70, 500 }, {}, false ).mRedo );
            CHECK( not decide_feature_redetect( tThresholds, {}, {}, false ).mRedo );
        }

        SECTION( "fires above the deck threshold, never for a geometry that is off" )
        {
            // the live 3D run (c9_geomfix A_cont): 32 hidden edges at eval 1, 65 at eval 20
            CHECK( not decide_feature_redetect( tThresholds, { 32, 900 }, {}, false ).mRedo );
            CHECK( not decide_feature_redetect( tThresholds, { 40, 900 }, {}, false ).mRedo );

            Feature_Redetect_Decision tDecision = decide_feature_redetect( tThresholds, { 65, 900 }, {}, false );
            CHECK( tDecision.mRedo );
            CHECK( tDecision.mGeometry == 0 );
            CHECK( tDecision.mHiddenEdges == 65 );
            CHECK( tDecision.mThreshold == 40 );
        }

        SECTION( "threshold 0: any hidden edge" )
        {
            CHECK( decide_feature_redetect( { 0 }, { 1 }, {}, false ).mRedo );
            CHECK( not decide_feature_redetect( { 0 }, { 0 }, {}, false ).mRedo );
        }

        SECTION( "never on the mesh refined for this design: a redo would rebuild the same mesh" )
        {
            CHECK( not decide_feature_redetect( tThresholds, { 65, 0 }, {}, true ).mRedo );
        }

        SECTION( "only above the count of the first cut after the last refinement (no re-firing next to a redo)" )
        {
            // after a redo the cut has 50 hidden edges: a perturbed design with 51 does not fire again, 75 does
            Vector< uint > tBaseline = { 50, 0 };
            CHECK( not decide_feature_redetect( tThresholds, { 50, 0 }, tBaseline, false ).mRedo );
            CHECK( not decide_feature_redetect( tThresholds, { 51, 0 }, { 51, 0 }, false ).mRedo );

            Feature_Redetect_Decision tDecision = decide_feature_redetect( tThresholds, { 75, 0 }, tBaseline, false );
            CHECK( tDecision.mRedo );
            CHECK( tDecision.mThreshold == 50 );

            // a baseline below the deck threshold does not lower it
            CHECK( not decide_feature_redetect( tThresholds, { 35, 0 }, { 10, 0 }, false ).mRedo );
        }

        SECTION( "several geometries: the one with the largest excess over its threshold" )
        {
            Feature_Redetect_Decision tDecision = decide_feature_redetect( { 10, 5 }, { 20, 30 }, {}, false );
            CHECK( tDecision.mRedo );
            CHECK( tDecision.mGeometry == 1 );
            CHECK( tDecision.mHiddenEdges == 30 );
            CHECK( tDecision.mThreshold == 5 );
        }
    }
}    // namespace moris::wrk
