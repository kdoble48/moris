/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_Library_Interlink_Checks.hpp
 *
 */

#pragma once

#include <string>

#include "cl_Module_Parameter_Lists.hpp"
#include "cl_Vector.hpp"

namespace moris
{
    /**
     * Cross-module consistency checks on a finalized deck's parameter lists (Deck API
     * v2, Stage 5 — see doc/internal/DECK_API_RFC.md). Covers the string-coupling
     * contract points that fail silently or cryptically at run time today:
     *
     *  - FEM name references resolve: IWG/IQI leader/follower properties,
     *    constitutive models, and stabilization parameters name a defined
     *    property/CM/SP (contract C4)
     *  - GEN IQI_types name defined FEM IQIs (contract C1, producer side)
     *  - VIS IQI_Names name defined FEM IQIs, and the Field_Names / Field_Type /
     *    IQI_Names lists have matching lengths (contract C8)
     *  - FEM mesh_set_names referencing ghost_p* require XTK ghost_stab (C5 slice)
     *
     * Returns human-readable findings (with did-you-mean suggestions for near-miss
     * names); empty when the deck is consistent. Pure function — the caller decides
     * whether findings are warnings or errors.
     *
     * @param aParameterLists All module parameter lists, indexed by Module_Type
     * @return Findings, empty if consistent
     */
    Vector< std::string >
    collect_deck_interlink_findings( const Vector< Module_Parameter_Lists >& aParameterLists );

}    // namespace moris
