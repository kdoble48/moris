/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Field_Evaluation_Log.hpp
 *
 */

#pragma once

#include <array>
#include <fstream>
#include <map>
#include <string>
#include <tuple>

#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Vector.hpp"

namespace moris::gen
{
    class Derived_Node;

    /**
     * Diagnostic log of the evaluations of an analytic level-set field (GEN parameter log_field_evaluations).
     *
     * Every value an analytic field computes for its geometry goes through one of Field_Analytic's three entry points
     * (a background node, a derived node, an arbitrary point). With a log attached, each one records the point, the
     * value and three tags:
     * - purpose: what was evaluated. Node evaluations are tagged by the node kind (background_node, centre_node,
     *   face_centre_node, intersection_node, derived_other_node); point evaluations by the Evaluation_Purpose set at
     *   the call site (root_iterate, detector_lattice, query_predictor, ...; point_untagged if none);
     * - stage: why, set by the caller in XTK / WRK (xtk_cell_selection, refine_feature_pass, ...; "-" if none);
     * - operation: the GEN operation that asked (region, analytic_root, linear_crossing, ...; "-" if none);
     * plus the node index (-1 for points), the GEN refinement pass (-1 outside the refinement passes), the geometry
     * index and the log instance (a geometry engine re-created in the same run appends to the same file).
     *
     * Rows are written only inside the optional box; the counts of all evaluations, inside and outside the box, go to
     * <file>.summary.csv when the log is destroyed. One file per processor: "_rank<r>" is inserted before the extension.
     * Off by default; with no log attached the only cost is one pointer test per evaluation.
     */
    class Field_Evaluation_Log
    {
      public:
        /**
         * @param aFileName Output CSV path; "_rank<r>" is inserted before the extension
         * @param aBox Optional "xmin,xmax,ymin,ymax[,zmin,zmax]"; empty = log everywhere
         * @param aGeometryIndex Index of the geometry in the geometry engine
         * @param aGeometryName Name of the geometry (summary only)
         */
        Field_Evaluation_Log(
                const std::string& aFileName,
                const std::string& aBox,
                uint               aGeometryIndex,
                const std::string& aGeometryName );

        ~Field_Evaluation_Log();

        Field_Evaluation_Log( const Field_Evaluation_Log& )            = delete;
        Field_Evaluation_Log& operator=( const Field_Evaluation_Log& ) = delete;

        /**
         * A background (Lagrange) node value.
         */
        void log_background_node( uint aNodeIndex, const Matrix< DDRMat >& aCoordinates, real aValue );

        /**
         * A derived node value (subdivision centre node, intersection node, ...), tagged by the node's kind.
         */
        void log_derived_node( const Derived_Node& aDerivedNode, real aValue );

        /**
         * A value at an arbitrary point, tagged by the current Evaluation_Purpose.
         */
        void log_point( const Matrix< DDRMat >& aCoordinates, real aValue );

        /**
         * Kind of a derived node: intersection_node; centre_node (2^d locators with equal weights: the centre node of a
         * regular subdivision); face_centre_node (3D, 4 locators with equal weights); derived_other_node.
         */
        static const char* derived_node_kind( const Derived_Node& aDerivedNode );

        /**
         * Name of the file this processor writes (aFileName with "_rank<r>" inserted before the extension).
         */
        static std::string rank_file_name( const std::string& aFileName );

        /**
         * Name of the file this processor writes.
         */
        const std::string& get_file_name() const;

        /**
         * Writes the summary now (also done by the destructor) and flushes the rows.
         */
        void write_summary();

      private:
        void record( const char* aPurpose, sint aNodeIndex, const Matrix< DDRMat >& aCoordinates, real aValue );

        std::string    mFileName;
        std::string    mGeometryName;
        uint           mGeometryIndex;
        uint           mInstance = 0;
        Vector< real > mBox;
        std::ofstream  mFile;

        // evaluations per (purpose, stage, operation, pass) tag pointers: { all, inside the box }
        std::map< std::tuple< const char*, const char*, const char*, sint >, std::array< luint, 2 > > mCounts;
    };

    /**
     * The tags of the evaluations currently made (one process, one thread; restored by the guards below).
     */
    struct Field_Evaluation_Context
    {
        static inline const char* sPurpose   = nullptr;
        static inline const char* sStage     = nullptr;
        static inline const char* sOperation = nullptr;
        static inline sint        sPass      = -1;
    };

    /**
     * Sets one tag of Field_Evaluation_Context for its lifetime (string literals only). With aOverride = false the tag
     * is set only if the caller has not set one: a default at a generic GEN site that an XTK / WRK caller can refine.
     */
    class Evaluation_Tag_Guard
    {
      public:
        Evaluation_Tag_Guard( const char*& aSlot, const char* aTag, bool aOverride )
                : mSlot( aSlot )
                , mPrevious( aSlot )
        {
            if ( aOverride or aSlot == nullptr )
            {
                aSlot = aTag;
            }
        }

        ~Evaluation_Tag_Guard()
        {
            mSlot = mPrevious;
        }

        Evaluation_Tag_Guard( const Evaluation_Tag_Guard& )            = delete;
        Evaluation_Tag_Guard& operator=( const Evaluation_Tag_Guard& ) = delete;

      private:
        const char*& mSlot;
        const char*  mPrevious;
    };

    /** Purpose of the point evaluations in scope (see Field_Evaluation_Log). */
    class Evaluation_Purpose : public Evaluation_Tag_Guard
    {
      public:
        explicit Evaluation_Purpose( const char* aTag, bool aOverride = true )
                : Evaluation_Tag_Guard( Field_Evaluation_Context::sPurpose, aTag, aOverride )
        {
        }
    };

    /** Stage (the caller's reason) of the evaluations in scope. */
    class Evaluation_Stage : public Evaluation_Tag_Guard
    {
      public:
        explicit Evaluation_Stage( const char* aTag, bool aOverride = true )
                : Evaluation_Tag_Guard( Field_Evaluation_Context::sStage, aTag, aOverride )
        {
        }
    };

    /** GEN operation of the evaluations in scope. */
    class Evaluation_Operation : public Evaluation_Tag_Guard
    {
      public:
        explicit Evaluation_Operation( const char* aTag, bool aOverride = true )
                : Evaluation_Tag_Guard( Field_Evaluation_Context::sOperation, aTag, aOverride )
        {
        }
    };

    /** GEN refinement pass of the evaluations in scope. */
    class Evaluation_Pass
    {
      public:
        explicit Evaluation_Pass( sint aPass )
                : mPrevious( Field_Evaluation_Context::sPass )
        {
            Field_Evaluation_Context::sPass = aPass;
        }

        ~Evaluation_Pass()
        {
            Field_Evaluation_Context::sPass = mPrevious;
        }

        Evaluation_Pass( const Evaluation_Pass& )            = delete;
        Evaluation_Pass& operator=( const Evaluation_Pass& ) = delete;

      private:
        sint mPrevious;
    };
}    // namespace moris::gen
