/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Connectivity_Checker.hpp
 *
 * An independent check of the three connectivity layers XTK builds on top of the conformal cut, computed from the
 * raw integration-mesh geometry rather than from XTK's facet tables:
 *
 *   1. IG cells -> subphases: union-find over same-phase IG tets of one background cell sharing a triangle;
 *   2. the subphase neighbourhood: two subphases in neighbouring background cells are adjacent iff they have the same
 *      bulk phase and their footprints on the shared face overlap with positive area;
 *   3. enrichment levels: per basis function, the connected components of the subphases in its support.
 *
 * It also checks that every ghost facet joins two subphases that touch. XTK's result must equal the reference as a
 * partition; a level holding two reference components couples disconnected material through one coefficient
 * (cross-talk). A standalone tool: construct a Connectivity_Checker on any decomposed (and enriched, ghost-stabilized)
 * xtk::Model with a 3D HMR background mesh and call check_and_report(). The XTK connectivity tests use it as their
 * reference.
 *
 */

#pragma once

#include <array>
#include <cstdint>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"

namespace moris::xtk
{
    class Model;
    class Cut_Integration_Mesh;
}    // namespace moris::xtk

namespace moris::xtk::connectivity_check
{
    using Point = std::array< real, 3 >;

    struct Hash
    {
        std::uint64_t mValue = 1469598103934665603ULL;

        void
        add( std::int64_t aValue )
        {
            for ( uint iByte = 0; iByte < 8; iByte++ )
            {
                mValue ^= ( static_cast< std::uint64_t >( aValue ) >> ( 8 * iByte ) ) & 0xffULL;
                mValue *= 1099511628211ULL;
            }
        }

        real
        high() const
        {
            return (real)( mValue >> 32 );
        }

        real
        low() const
        {
            return (real)( mValue & 0xffffffffULL );
        }
    };

    struct Union_Find
    {
        Vector< uint > mParent;

        explicit Union_Find( uint aSize )
                : mParent( aSize )
        {
            std::iota( mParent.begin(), mParent.end(), 0 );
        }

        uint
        find( uint aElement )
        {
            while ( mParent( aElement ) != aElement )
            {
                mParent( aElement ) = mParent( mParent( aElement ) );
                aElement            = mParent( aElement );
            }
            return aElement;
        }

        void
        unite( uint aFirst, uint aSecond )
        {
            uint tFirst  = this->find( aFirst );
            uint tSecond = this->find( aSecond );
            if ( tFirst != tSecond )
            {
                mParent( std::max( tFirst, tSecond ) ) = std::min( tFirst, tSecond );
            }
        }

        // component label per element (the smallest element of its component)
        Vector< moris_index >
        labels()
        {
            Vector< moris_index > tLabels( mParent.size() );
            for ( uint iElement = 0; iElement < mParent.size(); iElement++ )
            {
                tLabels( iElement ) = this->find( iElement );
            }
            return tLabels;
        }
    };

    uint
    count_distinct( const Vector< moris_index >& aLabels );

    // two labellings of the same elements: mMerged counts test labels that span several reference labels (cross-talk),
    // mSplit counts reference labels that span several test labels
    struct Partition_Comparison
    {
        uint mMerged = 0;
        uint mSplit  = 0;

        bool
        equal() const
        {
            return mMerged == 0 and mSplit == 0;
        }
    };

    Partition_Comparison
    compare_partitions( const Vector< moris_index >& aReference, const Vector< moris_index >& aTest );

    using Graph = std::vector< std::vector< moris_index > >;

    // neighbour matrix as the enrichment builds it: one row per node, unused columns hold the dummy MORIS_INDEX_MAX
    Matrix< IndexMat >
    graph_matrix( const Graph& aGraph, uint aNumColumns );

    struct Replica_Result
    {
        Vector< moris_index > mLabels;
        uint                  mFrontSize = 0;        // the allocation of flood_fill: 10 N + 1
        uint                  mPeakFront = 0;        // largest position written + 1
        bool                  mOverflow  = false;    // a write at or past mFrontSize
    };

    void
    write_front( std::vector< moris_index >& aFront, uint aPosition, moris_index aValue, Replica_Result& aResult );

    // a statement-by-statement copy of flood_fill (fn_mesh_flood_fill.hpp:38-257, aIncludeAllElements = true, every
    // element included) whose active front grows instead of overflowing, so the peak front can be measured
    Replica_Result
    flood_fill_replica( const Matrix< IndexMat >& aGraph, const Vector< moris_index >& aPhases );

    uint
    max_degree( const Graph& aGraph );

    using Key = std::pair< std::uint64_t, std::uint64_t >;    // ( high, low ) of a hash, as integers

    Key
    key_of( const Hash& aHash );

    struct Box
    {
        Point mMin = { MORIS_REAL_MAX, MORIS_REAL_MAX, MORIS_REAL_MAX };
        Point mMax = { -MORIS_REAL_MAX, -MORIS_REAL_MAX, -MORIS_REAL_MAX };

        real
        size( uint aDim ) const
        {
            return mMax[ aDim ] - mMin[ aDim ];
        }

        real
        volume() const
        {
            return this->size( 0 ) * this->size( 1 ) * this->size( 2 );
        }
    };

    std::array< std::int64_t, 3 >
    rounded_centroid( const Box& aBox );

    Box
    cell_box( const mtk::Cell& aCell );

    // the shared face of two background cells: the axis normal to it, its coordinate, and the overlap rectangle in the
    // two other axes ( ( a + 1 ) % 3, ( a + 2 ) % 3 )
    struct Shared_Face
    {
        uint                  mAxis  = 0;
        real                  mPlane = 0.0;
        std::array< real, 4 > mRectangle = { 0.0, 0.0, 0.0, 0.0 };    // umin, umax, vmin, vmax

        real
        area() const
        {
            return ( mRectangle[ 1 ] - mRectangle[ 0 ] ) * ( mRectangle[ 3 ] - mRectangle[ 2 ] );
        }
    };

    bool
    shared_face( const Box& aFirst, const Box& aSecond, Shared_Face& aFace );

    using Polygon = std::vector< std::array< real, 2 > >;

    real
    polygon_area( const Polygon& aPolygon );

    // counter-clockwise copy
    Polygon
    counter_clockwise( const Polygon& aPolygon );

    // Sutherland-Hodgman: aSubject clipped by the convex counter-clockwise aClip
    Polygon
    clip_polygon( const Polygon& aSubject, const Polygon& aClip );

    Polygon
    rectangle_polygon( const std::array< real, 4 >& aRectangle );

    // the corner vertices of an IG cell (TET4, or the corners of a TET10)
    Vector< mtk::Vertex* >
    corner_vertices( const mtk::Cell& aCell );

    struct Connectivity_Counts
    {
        uint mCases = 0;    // when accumulated over a sweep: decomposed cases

        // layer 1: subphases
        uint mCutCells                  = 0;
        uint mSubphases                 = 0;
        uint mSubphasePartitionMismatch = 0;    // cut cells whose subphases differ from the reference
        uint mCoverage                  = 0;    // IG cells in no subphase, in two, or outside their cell's group
        uint mMixedPhaseSubphase        = 0;    // subphases whose IG cells differ in bulk phase
        uint mVolumeMismatch            = 0;    // background cells whose subphase volumes do not add up to the cell
        real mMaxVolumeError            = 0.0;

        // layer 2: neighbourhood
        uint mLinks                = 0;    // XTK links (directed)
        uint mReferenceLinks       = 0;    // reference links (directed)
        uint mMissingLinks         = 0;    // reference link that XTK lacks
        uint mExtraLinks           = 0;    // XTK link without reference contact (and not marginal)
        uint mExtraNoContact       = 0;    // ... of which the two subphases have no footprint overlap at all
        uint mExtraTransition      = 0;    // ... of which across a refinement transition
        uint mMarginalLinks        = 0;    // XTK links whose contact area is below the tolerance (degenerate facets)
        uint mSelfLinks            = 0;
        uint mDuplicateLinks       = 0;
        uint mAsymmetricLinks      = 0;
        uint mCrossPhaseLinks      = 0;
        uint mSameCellLinks        = 0;
        uint mNonAdjacentCells     = 0;    // XTK link between background cells that share no face
        uint mBadSideOrdinals      = 0;    // side ordinal whose vertices are not on the shared face
        uint mNonconformingFaces   = 0;    // same-level cut faces whose triangles overlap without matching
        uint mDuplicateVertexTriangles = 0;    // identical face triangles (by coordinates) built from different vertices
        uint mTransitionFaces      = 0;    // shared faces between cells of different refinement level

        // layer 3: enrichment
        uint mBases                = 0;
        uint mEnrichedBases        = 0;    // sum over bases of the XTK level count
        uint mReferenceComponents  = 0;    // sum over bases of the reference component count
        uint mCrossTalk            = 0;    // bases with a level spanning several reference components
        uint mSplitComponents      = 0;    // bases with a reference component split over several levels
        uint mLevelCountMismatch   = 0;    // bases whose level count differs from the component count
        uint mMissingBasis         = 0;    // subphase in support without an entry for the basis
        uint mEnrichedTotalMismatch = 0;   // mNumEnrichedBasisFunctions differs from the sum of level counts
        uint mZeroRowBases         = 0;    // enriched bases whose subphases have no positive integration volume (E3)
        uint mTinyRowBases         = 0;    // ... below 1e-12 of a background cell
        uint mMaxSupportSubphases  = 0;
        uint mMaxGraphDegree       = 0;    // widest row of the pruned subphase graph (50 columns allocated in Enrichment::generate_pruned_subphase_graph_in_basis_support)
        uint mFrontOverflow        = 0;    // pruned graphs whose flood-fill front would exceed 10 N + 1
        real mMaxFrontRatio        = 0.0;

        // bodies per phase, on this rank (only meaningful on one rank)
        uint mBodyMismatch = 0;    // phases whose XTK body count differs from the reference

        // ghost (E10)
        uint mGhostPairs         = 0;
        uint mGhostWithoutContact = 0;    // ghost facet between subphases that the reference does not link

        std::string mFirstFailure;

        void
        note( const std::string& aWhat )
        {
            if ( mFirstFailure.empty() )
            {
                mFirstFailure = aWhat;
            }
        }

        uint
        failures() const
        {
            return mSubphasePartitionMismatch + mCoverage + mMixedPhaseSubphase + mVolumeMismatch + mMissingLinks + mExtraLinks + mSelfLinks
                 + mDuplicateLinks + mAsymmetricLinks + mCrossPhaseLinks + mSameCellLinks + mNonAdjacentCells + mBadSideOrdinals + mNonconformingFaces
                 + mCrossTalk + mSplitComponents + mLevelCountMismatch + mMissingBasis + mEnrichedTotalMismatch + mGhostWithoutContact;
        }

        void
        add( const Connectivity_Counts& aOther )
        {
            mCases += 1;
            mCutCells += aOther.mCutCells;
            mSubphases += aOther.mSubphases;
            mSubphasePartitionMismatch += aOther.mSubphasePartitionMismatch;
            mCoverage += aOther.mCoverage;
            mMixedPhaseSubphase += aOther.mMixedPhaseSubphase;
            mVolumeMismatch += aOther.mVolumeMismatch;
            mMaxVolumeError = std::max( mMaxVolumeError, aOther.mMaxVolumeError );
            mLinks += aOther.mLinks;
            mReferenceLinks += aOther.mReferenceLinks;
            mMissingLinks += aOther.mMissingLinks;
            mExtraLinks += aOther.mExtraLinks;
            mExtraNoContact += aOther.mExtraNoContact;
            mExtraTransition += aOther.mExtraTransition;
            mMarginalLinks += aOther.mMarginalLinks;
            mSelfLinks += aOther.mSelfLinks;
            mDuplicateLinks += aOther.mDuplicateLinks;
            mAsymmetricLinks += aOther.mAsymmetricLinks;
            mCrossPhaseLinks += aOther.mCrossPhaseLinks;
            mSameCellLinks += aOther.mSameCellLinks;
            mNonAdjacentCells += aOther.mNonAdjacentCells;
            mBadSideOrdinals += aOther.mBadSideOrdinals;
            mNonconformingFaces += aOther.mNonconformingFaces;
            mDuplicateVertexTriangles += aOther.mDuplicateVertexTriangles;
            mTransitionFaces += aOther.mTransitionFaces;
            mBases += aOther.mBases;
            mEnrichedBases += aOther.mEnrichedBases;
            mReferenceComponents += aOther.mReferenceComponents;
            mCrossTalk += aOther.mCrossTalk;
            mSplitComponents += aOther.mSplitComponents;
            mLevelCountMismatch += aOther.mLevelCountMismatch;
            mMissingBasis += aOther.mMissingBasis;
            mEnrichedTotalMismatch += aOther.mEnrichedTotalMismatch;
            mZeroRowBases += aOther.mZeroRowBases;
            mTinyRowBases += aOther.mTinyRowBases;
            mMaxSupportSubphases = std::max( mMaxSupportSubphases, aOther.mMaxSupportSubphases );
            mMaxGraphDegree      = std::max( mMaxGraphDegree, aOther.mMaxGraphDegree );
            mFrontOverflow += aOther.mFrontOverflow;
            mMaxFrontRatio = std::max( mMaxFrontRatio, aOther.mMaxFrontRatio );
            mBodyMismatch += aOther.mBodyMismatch;
            mGhostPairs += aOther.mGhostPairs;
            mGhostWithoutContact += aOther.mGhostWithoutContact;
        }

        // one line for the log: what was checked and every nonzero failure count
        std::string
        summary() const;

        // the same fields summed over the ranks (maxima reduced as maxima)
        void
        reduce()
        {
            Vector< uint* > tSums = { &mCutCells, &mSubphases, &mSubphasePartitionMismatch, &mCoverage, &mMixedPhaseSubphase, &mVolumeMismatch,
                &mLinks, &mReferenceLinks, &mMissingLinks, &mExtraLinks, &mExtraNoContact, &mExtraTransition, &mMarginalLinks, &mSelfLinks,
                &mDuplicateLinks, &mAsymmetricLinks, &mCrossPhaseLinks, &mSameCellLinks, &mNonAdjacentCells, &mBadSideOrdinals,
                &mNonconformingFaces, &mDuplicateVertexTriangles, &mTransitionFaces, &mBases, &mEnrichedBases, &mReferenceComponents, &mCrossTalk, &mSplitComponents,
                &mLevelCountMismatch, &mMissingBasis, &mEnrichedTotalMismatch, &mZeroRowBases, &mTinyRowBases, &mFrontOverflow, &mBodyMismatch,
                &mGhostPairs, &mGhostWithoutContact };
            for ( uint* iValue : tSums )
            {
                *iValue = sum_all( *iValue );
            }
            mMaxVolumeError      = max_all( mMaxVolumeError );
            mMaxSupportSubphases = max_all( mMaxSupportSubphases );
            mMaxGraphDegree      = max_all( mMaxGraphDegree );
            mMaxFrontRatio       = max_all( mMaxFrontRatio );
        }
    };

    // reference status of a pair of subphases
    enum class Contact
    {
        NONE,        // no overlap of the footprints
        MARGINAL,    // overlap area in ( 0, 1e-12 |face| ], or footprints made of degenerate triangles only
        LINK         // identical triangles on a conforming face, or overlap area above 1e-12 |face|
    };

    // a subphase's piece of a background face: a triangle of one of its IG cells, or the whole face (uncut cell)
    struct Footprint_Piece
    {
        moris_index                  mSubphase = 0;
        Polygon                      mPolygon;
        std::array< moris_index, 3 > mVertices  = { -1, -1, -1 };    // sorted IG vertex indices (triangles only)
        std::array< std::array< std::int64_t, 3 >, 3 > mCoordinates = {};    // sorted rounded vertex coordinates (triangles only)
        bool                         mWholeFace = false;
    };

    // data gathered for the multi-rank comparison
    struct Parallel_Rows
    {
        Vector< Vector< real > > mSubphases;           // subphase id, owned, key high, key low, background cell id, bulk phase
        Vector< Vector< real > > mEdges;               // owned subphase: key (high, low) of it and of the neighbour
        Vector< Vector< real > > mBases;               // owned basis: id, level count, hash of the levels (high, low)
        Vector< Vector< real > > mEnrichedIds;         // subphase id, basis id, enriched basis id, subphase owned
        Vector< Vector< real > > mClusterSubphases;    // subphase id of every cluster's interpolation cell on this rank
        Vector< Vector< real > > mOwnedSubphaseIds;    // owned subphase id
    };

    // the reference check of one decomposed (and optionally enriched, ghost-stabilized) XTK model on a 3D HMR background
    // mesh; the results are in mCounts (local to this rank until mCounts.reduce())
    class Connectivity_Checker
    {
      public:
        mtk::Interpolation_Mesh* mBgMesh    = nullptr;
        Model*                   mModel     = nullptr;
        Cut_Integration_Mesh*    mCutMesh   = nullptr;
        moris_index              mMeshIndex = 0;    // enriched B-spline mesh index currently checked
        uint                     mNumPhases = 0;
        uint                     mMaxDumps  = 0;    // log the vertices of at most this many nonconforming face-triangle pairs

        Connectivity_Counts mCounts;
        Vector< uint >      mReferenceBodies;    // connected components per bulk phase (reference links)
        Vector< uint >      mXtkBodies;          // ... with XTK's links

        Vector< Box >                                                  mBoxes;
        std::map< std::pair< moris_index, moris_index >, Shared_Face > mFaces;      // background cell pair (smaller first)
        std::map< std::pair< moris_index, moris_index >, Contact >     mContact;    // subphase pair (both orders)
        std::set< std::pair< moris_index, moris_index > >              mXtkLinks;
        Vector< Key >                                                  mSubphaseKeys;
        uint                                                           mNumDumped = 0;

        Connectivity_Checker( mtk::Interpolation_Mesh* aBgMesh, Model* aModel, uint aNumPhases, moris_index aMeshIndex );

        // run every check: subphases, neighbourhood, bodies, then (aEnrichment) the enrichment levels of each mesh index
        // and (aGhost) the ghost facets; aRows (may be null) collects canonical rows for a multi-rank comparison
        void
        check_all( const Vector< moris_index >& aMeshIndices, bool aEnrichment, bool aGhost, Parallel_Rows* aRows );

        // check_all, then sum the counts over all ranks, log a one-line summary and each rank's first inconsistency, and
        // warn if there is any (collective); returns the summed counts
        Connectivity_Counts
        check_and_report( const Vector< moris_index >& aMeshIndices, bool aEnrichment, bool aGhost );

        // background cell boxes and the shared faces of face neighbours (from the mesh's neighbour lists)
        void
        compute_boxes_and_faces();

        // layer 1: each cut cell's subphases against union-find over its same-phase IG cells sharing a triangle
        void
        check_subphases();

        // canonical subphase keys: background cell id, bulk phase, rounded IG cell centroids (independent of ranks)
        void
        compute_subphase_keys();

        // layer 2 reference: each subphase's pieces on a shared face, and the contact of every subphase pair
        Vector< Footprint_Piece >
        footprint( moris_index aCell, const Shared_Face& aFace );

        void
        raise_contact( moris_index aFirst, moris_index aSecond, Contact aContact );

        void
        compute_reference_contacts();

        void
        dump_triangles( const std::string& aWhat, moris_index aFirstCell, moris_index aSecondCell, const Footprint_Piece& aA, const Footprint_Piece& aB );

        // one side of the pair is an uncut whole face and the other side's pieces are all of (near) zero area
        bool
        has_only_degenerate_pieces(
                const std::pair< moris_index, moris_index >& aPair,
                const Vector< Footprint_Piece >&             aFirst,
                const Vector< Footprint_Piece >&             aSecond,
                real                                         aFaceArea );

        Contact
        contact( moris_index aFirst, moris_index aSecond ) const;

        bool
        side_on_face( const mtk::Cell& aCell, moris_index aOrdinal, const Shared_Face& aFace );

        // layer 2: XTK's subphase neighbourhood against the reference contacts
        void
        check_neighborhood();

        // bodies: connected components of all subphases, per phase, with the reference and with XTK's links
        void
        check_bodies();

        // layer 3: the enrichment levels of each basis function against the components of its support
        void
        check_enrichment( Parallel_Rows* aRows );

        // every ghost facet joins two subphases that touch
        void
        check_ghost();

        void
        collect_parallel_rows( Parallel_Rows& aRows, bool aEnrichment );
    };
}    // namespace moris::xtk::connectivity_check
