/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Input_Deck_Vocabulary.hpp
 *
 * Typed mesh-set vocabulary for input decks (Deck API v2): Phase / Side handles
 * generate the XTK set-name grammar (HMR_dummy_*, SideSet_*, iside_*, dbl_iside_*,
 * ghost_*) so decks never hand-type — or mistype — set-name strings (coupling
 * contract C5). Mirrors the pymoris fluent vocabulary (MeshSet.bulk_and_cut(i),
 * DOF.DISPLACEMENT_2D, ...) so knowledge transfers between the two layers.
 */

#pragma once

#include <string>

#include "moris_typedefs.hpp"

namespace moris::deck
{
    //------------------------------------------------------------------------------------------------------------------

    /**
     * Domain-boundary sides, numbered per the HMR/Exodus side-set convention
     * (2D quad and 3D hex): 1 = y-min, 2 = x-max, 3 = y-max, 4 = x-min.
     * The 3D z-faces follow the hex convention (5 = z-min, 6 = z-max).
     */
    enum class Side : uint
    {
        Bottom = 1,    // y-min
        Right  = 2,    // x-max
        Top    = 3,    // y-max
        Left   = 4,    // x-min
        Back   = 5,    // z-min
        Front  = 6     // z-max
    };

    //------------------------------------------------------------------------------------------------------------------

    /**
     * A bulk phase of the cut (XTK) mesh, identified by its bulk-phase index.
     * Methods generate the mesh-set name strings FEM/VIS consume.
     */
    struct Phase
    {
        uint mIndex;

        explicit Phase( uint aIndex )
                : mIndex( aIndex )
        {
        }

        /** Bulk sets: non-cut and cut background cells of this phase. */
        std::string
        bulk() const
        {
            std::string tPhase = std::to_string( mIndex );
            return "HMR_dummy_n_p" + tPhase + ",HMR_dummy_c_p" + tPhase;
        }

        /** Domain-boundary side sets of this phase on a given side. */
        std::string
        side( Side aSide ) const
        {
            std::string tSide  = std::to_string( (uint)aSide );
            std::string tPhase = std::to_string( mIndex );
            return "SideSet_" + tSide + "_n_p" + tPhase + ",SideSet_" + tSide + "_c_p" + tPhase;
        }

        /** Ghost facet set of this phase (requires XTK ghost_stab). */
        std::string
        ghost() const
        {
            return "ghost_p" + std::to_string( mIndex );
        }
    };

    //------------------------------------------------------------------------------------------------------------------

    /**
     * Single-sided interface set: facets of aInside's cells facing aOutside
     * (iside_b0_<in>_b1_<out>).
     */
    inline std::string
    interface( const Phase& aInside, const Phase& aOutside )
    {
        return "iside_b0_" + std::to_string( aInside.mIndex ) + "_b1_" + std::to_string( aOutside.mIndex );
    }

    /**
     * Double-sided interface set between a leader and a follower phase
     * (dbl_iside_p0_<leader>_p1_<follower>).
     */
    inline std::string
    between( const Phase& aLeader, const Phase& aFollower )
    {
        return "dbl_iside_p0_" + std::to_string( aLeader.mIndex ) + "_p1_" + std::to_string( aFollower.mIndex );
    }

    /** Comma-joins set-name lists (set union in the mesh-set grammar). */
    inline std::string
    join( const std::string& aFirst, const std::string& aSecond )
    {
        return aFirst + "," + aSecond;
    }

    //------------------------------------------------------------------------------------------------------------------

    /**
     * DOF-list constants, mirroring the pymoris DOF vocabulary.
     */
    namespace Dofs
    {
        inline const std::string Displacement2D = "UX,UY";
        inline const std::string Displacement3D = "UX,UY,UZ";
        inline const std::string Temperature    = "TEMP";
        inline const std::string Pressure       = "P";
    }    // namespace Dofs

    //------------------------------------------------------------------------------------------------------------------

}    // namespace moris::deck
