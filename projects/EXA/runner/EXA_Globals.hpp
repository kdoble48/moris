/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * EXA_Globals.hpp
 *
 * Deck-visible globals: input .so files resolve these BY NAME against the
 * executable's -rdynamic dynamic symbol table at dlopen time. They must stay
 * at global scope with exactly these names/types, defined once (EXA_Globals.cpp).
 * Do NOT namespace, rename, or make static.
 *
 * Every TEST_CASE must explicitly set each global it or its deck consumes
 * (rule R4, EXA_RUNNER_RFC.md) - never rely on zero/file-scope initialization.
 *
 * Regenerate the union with:
 *   grep -rhP '^\s*extern\s+(?!"C")' projects/EXA --include='*.cpp' \
 *     | grep -v example_test_case | sort -u
 */

#pragma once

#include <string>

#include "moris_typedefs.hpp"    // MRS/COR/src

extern moris::uint gInterpolationOrder;
extern moris::uint gTestCaseIndex;
extern moris::uint gTestIndex;
extern moris::uint gCaseIndex;
extern moris::uint gDim;
extern moris::uint gOrder;
extern moris::uint tGeoModel;
extern moris::uint tDim;
extern moris::uint gLevelSetInterpolationOrder;
extern moris::uint gFEMInterpolationOrder;
extern moris::uint gLagrMeshInterpolationOrder;

extern bool gPrintReferenceValues;
extern bool gUseBspline;
extern bool gUseMixedTimeElements;
extern bool gUseBelosWithILUT;
extern bool gInletVelocityBCFlag;
extern bool gInletPressureBCFlag;
extern bool gHaveStaggeredFA;
extern bool gHaveStaggeredSA;

extern std::string tOrder;
extern std::string tStressType;
extern std::string tOutputFileName;
extern std::string gPrecSolver;
