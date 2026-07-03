/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * EXA_Globals.cpp
 *
 * Single definition site for the deck-visible globals (see EXA_Globals.hpp).
 * Values here are meaningless on purpose: every TEST_CASE sets what it uses
 * (rule R4, share/doc/EXA_RUNNER_RFC.md).
 */

#include "EXA_Globals.hpp"

moris::uint gInterpolationOrder          = 0;
moris::uint gTestCaseIndex               = 0;
moris::uint gTestIndex                   = 0;
moris::uint gCaseIndex                   = 0;
moris::uint gDim                         = 0;
moris::uint gOrder                       = 0;
moris::uint tGeoModel                    = 0;
moris::uint tDim                         = 0;
moris::uint gLevelSetInterpolationOrder  = 0;
moris::uint gFEMInterpolationOrder       = 0;
moris::uint gLagrMeshInterpolationOrder  = 0;

bool gPrintReferenceValues = false;
bool gUseBspline           = false;
bool gUseMixedTimeElements = false;
bool gUseBelosWithILUT     = false;
bool gInletVelocityBCFlag  = false;
bool gInletPressureBCFlag  = false;
bool gHaveStaggeredFA      = false;
bool gHaveStaggeredSA      = false;

std::string tOrder;
std::string tStressType;
std::string tOutputFileName;
std::string gPrecSolver;
