/*
 * Vector-of-vectors growth benchmark for XTK edge connectivity.
 *
 * Replays the container pattern of xtk::Integration_Mesh_Generator::create_edges_from_element_to_node()
 * (cl_XTK_Integration_Mesh_Generator.cpp) on a structured tet mesh (n^3 hexes, 6 Kuhn tets each):
 *   - outer lists reserved to 20% of the maximum edge count, then grown with push_back( value, increment )
 *   - one inner list per edge (vertices, cells, cell edge ordinals), inner reserve of 6
 *   - per-vertex edge lists grown with push_back( edge, 10 ), lookup from the lower-id vertex (edge_exists())
 *   - shrink_to_fit_all() on the outer lists at the end
 * using the real moris::Vector from MRS/CNT/src/cl_Vector.hpp. Which header is picked up decides the variant:
 *   copy          : the header as in MORIS (copy ctor only -> outer reallocations deep-copy every inner list)
 *   move          : patched header with a noexcept move ctor / move assignment (patch_header.py)
 *   move_throwing : patched header compiled with -DBENCH_MOVE_NOT_NOEXCEPT (std::vector still copies)
 * and -DBENCH_CSR builds the same connectivity as flat arrays (fixed width + CSR) for comparison.
 *
 * Usage: bench <n_hex_per_side> <repeats>
 * Output: one CSV line per repeat:
 *   variant,n,tets,edges,build_s,shrink_s,total_s,allocs,peak_live_mb,maxrss_mb,nothrow_move,checksum
 */

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <sys/resource.h>

#include "cl_Vector.hpp"

// ------------------------------------------------------------------------------------------------
// allocation counting (global operator new/delete; malloc underneath, so tcmalloc applies if linked)

static std::size_t gNumAllocs  = 0;
static std::size_t gLiveBytes  = 0;
static std::size_t gPeakBytes  = 0;
static bool        gCountAlloc = false;

static void*
counted_alloc( std::size_t aSize )
{
    // store the size in front of the block so delete can account for it
    void* tBlock = std::malloc( aSize + 16 );
    if ( tBlock == nullptr )
    {
        throw std::bad_alloc();
    }
    *static_cast< std::size_t* >( tBlock ) = aSize;
    if ( gCountAlloc )
    {
        gNumAllocs++;
        gLiveBytes += aSize;
        gPeakBytes = std::max( gPeakBytes, gLiveBytes );
    }
    return static_cast< char* >( tBlock ) + 16;
}

static void
counted_free( void* aPtr )
{
    if ( aPtr == nullptr )
    {
        return;
    }
    void*       tBlock = static_cast< char* >( aPtr ) - 16;
    std::size_t tSize  = *static_cast< std::size_t* >( tBlock );
    if ( gCountAlloc )
    {
        gLiveBytes -= std::min( tSize, gLiveBytes );
    }
    std::free( tBlock );
}

void* operator new( std::size_t aSize ) { return counted_alloc( aSize ); }
void* operator new[]( std::size_t aSize ) { return counted_alloc( aSize ); }
void  operator delete( void* aPtr ) noexcept { counted_free( aPtr ); }
void  operator delete[]( void* aPtr ) noexcept { counted_free( aPtr ); }
void  operator delete( void* aPtr, std::size_t ) noexcept { counted_free( aPtr ); }
void  operator delete[]( void* aPtr, std::size_t ) noexcept { counted_free( aPtr ); }

// ------------------------------------------------------------------------------------------------
// stand-ins for mtk::Vertex / mtk::Cell: the connectivity stores pointers to them, as in XTK

struct Bench_Vertex
{
    int mId;
};

struct Bench_Cell
{
    int           mIndex;
    Bench_Vertex* mVerts[ 4 ];
};

using Clock = std::chrono::steady_clock;

static double
seconds_since( Clock::time_point aStart )
{
    return std::chrono::duration< double >( Clock::now() - aStart ).count();
}

static double
max_rss_mb()
{
    struct rusage tUsage;
    getrusage( RUSAGE_SELF, &tUsage );
    return tUsage.ru_maxrss / 1024.0;
}

// structured n^3 hex grid, each hex split into the 6 Kuhn tets around the diagonal v0-v6 (conforming)
static void
build_tet_mesh(
        int                          aN,
        std::vector< Bench_Vertex >& aVertices,
        std::vector< Bench_Cell >&   aCells )
{
    const int tNp1 = aN + 1;
    aVertices.resize( (std::size_t)tNp1 * tNp1 * tNp1 );
    for ( std::size_t i = 0; i < aVertices.size(); i++ )
    {
        aVertices[ i ].mId = (int)i;
    }

    // hex corner offsets: v0 (0,0,0) v1 (1,0,0) v2 (1,1,0) v3 (0,1,0) v4 (0,0,1) v5 (1,0,1) v6 (1,1,1) v7 (0,1,1)
    const int tCorner[ 8 ][ 3 ] = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 } };
    const int tKuhn[ 6 ][ 4 ]   = { { 0, 1, 2, 6 }, { 0, 2, 3, 6 }, { 0, 3, 7, 6 }, { 0, 7, 4, 6 }, { 0, 4, 5, 6 }, { 0, 5, 1, 6 } };

    aCells.reserve( (std::size_t)6 * aN * aN * aN );
    for ( int k = 0; k < aN; k++ )
    {
        for ( int j = 0; j < aN; j++ )
        {
            for ( int i = 0; i < aN; i++ )
            {
                int tHexVert[ 8 ];
                for ( int c = 0; c < 8; c++ )
                {
                    tHexVert[ c ] = ( i + tCorner[ c ][ 0 ] ) + tNp1 * ( ( j + tCorner[ c ][ 1 ] ) + tNp1 * ( k + tCorner[ c ][ 2 ] ) );
                }
                for ( int t = 0; t < 6; t++ )
                {
                    Bench_Cell tCell;
                    tCell.mIndex = (int)aCells.size();
                    for ( int v = 0; v < 4; v++ )
                    {
                        tCell.mVerts[ v ] = &aVertices[ tHexVert[ tKuhn[ t ][ v ] ] ];
                    }
                    aCells.push_back( tCell );
                }
            }
        }
    }
}

static const int gTetEdgeToNode[ 6 ][ 2 ] = { { 0, 1 }, { 1, 2 }, { 0, 2 }, { 0, 3 }, { 1, 3 }, { 2, 3 } };

struct Result
{
    std::size_t mNumEdges = 0;
    double      mBuild    = 0.0;
    double      mShrink   = 0.0;
    std::size_t mChecksum = 0;
};

#ifndef BENCH_CSR

// ------------------------------------------------------------------------------------------------
// XTK pattern with moris::Vector (copy / move / move_throwing variants)

struct Edge_Connectivity
{
    moris::Vector< moris::Vector< Bench_Vertex* > > mEdgeVertices;
    moris::Vector< moris::Vector< Bench_Cell* > >   mEdgeToCell;
    moris::Vector< moris::Vector< int > >           mEdgeToCellEdgeOrdinal;
    moris::Vector< moris::Vector< int > >           mCellToEdge;
};

// edge_exists(): edge vertices sorted by id, search the edges attached to the lower vertex
static int
edge_exists(
        moris::Vector< Bench_Vertex* >&                        aVerticesOnEdge,
        moris::Vector< moris::Vector< unsigned int > >&        aVertexToEdge,
        moris::Vector< moris::Vector< Bench_Vertex* > > const& aFullEdgeVertices )
{
    if ( aVerticesOnEdge( 1 )->mId < aVerticesOnEdge( 0 )->mId )
    {
        std::swap( aVerticesOnEdge( 0 ), aVerticesOnEdge( 1 ) );
    }
    moris::Vector< unsigned int > const & tEdges = aVertexToEdge( aVerticesOnEdge( 0 )->mId );
    for ( unsigned int iEdge = 0; iEdge < tEdges.size(); iEdge++ )
    {
        if ( aFullEdgeVertices( tEdges( iEdge ) )( 1 ) == aVerticesOnEdge( 1 ) )
        {
            return (int)tEdges( iEdge );
        }
    }
    return -1;
}

static Result
run_variant(
        std::vector< Bench_Vertex >& aVertices,
        std::vector< Bench_Cell >&   aCells )
{
    Result            tResult;
    Edge_Connectivity tConn;

    Clock::time_point tStart = Clock::now();

    const std::size_t tNumEdgePerElem = 6;
    tConn.mCellToEdge.resize( aCells.size() );

    // same capacity policy as XTK: increment = 20% of the maximum number of edges
    const std::size_t tMaxNumEdges = aCells.size() * tNumEdgePerElem;
    const std::size_t tIncNumEdges = tMaxNumEdges * 2 / 10;
    tConn.mEdgeToCell.reserve( tIncNumEdges );
    tConn.mEdgeToCellEdgeOrdinal.reserve( tIncNumEdges );
    tConn.mEdgeVertices.reserve( tIncNumEdges );

    moris::Vector< moris::Vector< unsigned int > > tVertexToEdgeIndex( aVertices.size() );
    moris::Vector< Bench_Vertex* >                 tVerticesOnEdge( 2, nullptr );

    for ( std::size_t i = 0; i < aCells.size(); i++ )
    {
        tConn.mCellToEdge( i ).reserve( tNumEdgePerElem );
        for ( int iEdge = 0; iEdge < 6; iEdge++ )
        {
            tVerticesOnEdge( 0 ) = aCells[ i ].mVerts[ gTetEdgeToNode[ iEdge ][ 0 ] ];
            tVerticesOnEdge( 1 ) = aCells[ i ].mVerts[ gTetEdgeToNode[ iEdge ][ 1 ] ];

            int tEdgeIndex = edge_exists( tVerticesOnEdge, tVertexToEdgeIndex, tConn.mEdgeVertices );

            if ( tEdgeIndex < 0 )
            {
                tEdgeIndex = (int)tConn.mEdgeVertices.size();
                tConn.mEdgeVertices.push_back( tVerticesOnEdge, tIncNumEdges );
                tConn.mEdgeToCell.push_back( moris::Vector< Bench_Cell* >(), tIncNumEdges );
                tConn.mEdgeToCellEdgeOrdinal.push_back( moris::Vector< int >(), tIncNumEdges );
                tConn.mEdgeToCell.back().reserve( tNumEdgePerElem );
                tConn.mEdgeToCellEdgeOrdinal.back().reserve( tNumEdgePerElem );
                tVertexToEdgeIndex( tVerticesOnEdge( 0 )->mId ).push_back( (unsigned int)tEdgeIndex, 10 );
            }

            tConn.mCellToEdge( i ).push_back( tEdgeIndex );
            tConn.mEdgeToCell( tEdgeIndex ).push_back( &aCells[ i ], tNumEdgePerElem );
            tConn.mEdgeToCellEdgeOrdinal( tEdgeIndex ).push_back( iEdge, tNumEdgePerElem );
        }
    }
    tResult.mBuild = seconds_since( tStart );

    tStart = Clock::now();
    shrink_to_fit_all( tConn.mEdgeVertices );
    shrink_to_fit_all( tConn.mEdgeToCell );
    shrink_to_fit_all( tConn.mEdgeToCellEdgeOrdinal );
    tResult.mShrink = seconds_since( tStart );

    // order-independent checksum: every (cell, ordinal) incidence of every edge
    tResult.mNumEdges = tConn.mEdgeVertices.size();
    for ( std::size_t e = 0; e < tConn.mEdgeToCell.size(); e++ )
    {
        for ( std::size_t c = 0; c < tConn.mEdgeToCell( e ).size(); c++ )
        {
            tResult.mChecksum += (std::size_t)tConn.mEdgeToCell( e )( c )->mIndex * 7 + tConn.mEdgeToCellEdgeOrdinal( e )( c );
        }
        tResult.mChecksum += tConn.mEdgeToCell( e ).size() * tConn.mEdgeToCell( e ).size();
    }
    return tResult;
}

#else

// ------------------------------------------------------------------------------------------------
// flat layout: edge vertices fixed width 2, cell-to-edge fixed width 6, edge-to-cell as CSR
// edges are identified by sorting the 6 * n_cells (lower id, upper id) keys once

static Result
run_variant(
        std::vector< Bench_Vertex >& aVertices,
        std::vector< Bench_Cell >&   aCells )
{
    Result            tResult;
    Clock::time_point tStart = Clock::now();

    const std::size_t tNumCellEdges = aCells.size() * 6;

    // one record per (cell, edge ordinal): packed vertex-pair key + cell-edge slot
    std::vector< std::pair< std::uint64_t, std::uint32_t > > tRecords( tNumCellEdges );
    for ( std::size_t i = 0; i < aCells.size(); i++ )
    {
        for ( int iEdge = 0; iEdge < 6; iEdge++ )
        {
            std::uint32_t tA = (std::uint32_t)aCells[ i ].mVerts[ gTetEdgeToNode[ iEdge ][ 0 ] ]->mId;
            std::uint32_t tB = (std::uint32_t)aCells[ i ].mVerts[ gTetEdgeToNode[ iEdge ][ 1 ] ]->mId;
            if ( tB < tA )
            {
                std::swap( tA, tB );
            }
            tRecords[ i * 6 + iEdge ] = { ( (std::uint64_t)tA << 32 ) | tB, (std::uint32_t)( i * 6 + iEdge ) };
        }
    }
    std::sort( tRecords.begin(), tRecords.end() );

    // runs of equal keys are edges; the sorted order is already grouped by edge -> CSR offsets
    std::vector< std::uint32_t > tEdgeToCellOffsets;
    std::vector< Bench_Vertex* > tEdgeVertices;    // 2 per edge
    std::vector< Bench_Cell* >   tEdgeToCell( tNumCellEdges );
    std::vector< int >           tEdgeToCellEdgeOrdinal( tNumCellEdges );
    std::vector< int >           tCellToEdge( tNumCellEdges );    // 6 per cell
    tEdgeToCellOffsets.reserve( tNumCellEdges / 4 + 1 );
    tEdgeVertices.reserve( tNumCellEdges / 2 );

    std::size_t tNumEdges = 0;
    for ( std::size_t r = 0; r < tRecords.size(); r++ )
    {
        if ( r == 0 || tRecords[ r ].first != tRecords[ r - 1 ].first )
        {
            tEdgeToCellOffsets.push_back( (std::uint32_t)r );
            tEdgeVertices.push_back( &aVertices[ tRecords[ r ].first >> 32 ] );
            tEdgeVertices.push_back( &aVertices[ tRecords[ r ].first & 0xffffffffu ] );
            tNumEdges++;
        }
        std::uint32_t tSlot         = tRecords[ r ].second;
        tEdgeToCell[ r ]            = &aCells[ tSlot / 6 ];
        tEdgeToCellEdgeOrdinal[ r ] = (int)( tSlot % 6 );
        tCellToEdge[ tSlot ]        = (int)( tNumEdges - 1 );
    }
    tEdgeToCellOffsets.push_back( (std::uint32_t)tRecords.size() );
    tResult.mBuild = seconds_since( tStart );

    tStart = Clock::now();
    std::vector< std::pair< std::uint64_t, std::uint32_t > >().swap( tRecords );
    tEdgeToCellOffsets.shrink_to_fit();
    tEdgeVertices.shrink_to_fit();
    tResult.mShrink = seconds_since( tStart );

    // same order-independent checksum as the nested variants
    tResult.mNumEdges = tNumEdges;
    for ( std::size_t e = 0; e < tNumEdges; e++ )
    {
        for ( std::uint32_t c = tEdgeToCellOffsets[ e ]; c < tEdgeToCellOffsets[ e + 1 ]; c++ )
        {
            tResult.mChecksum += (std::size_t)tEdgeToCell[ c ]->mIndex * 7 + tEdgeToCellEdgeOrdinal[ c ];
        }
        std::size_t tCount = tEdgeToCellOffsets[ e + 1 ] - tEdgeToCellOffsets[ e ];
        tResult.mChecksum += tCount * tCount;
    }
    (void)tCellToEdge;
    return tResult;
}

#endif

int
main( int argc, char** argv )
{
    if ( argc < 3 )
    {
        std::fprintf( stderr, "usage: %s <n_hex_per_side> <repeats>\n", argv[ 0 ] );
        return 1;
    }
    const int tN       = std::atoi( argv[ 1 ] );
    const int tRepeats = std::atoi( argv[ 2 ] );

#if defined( BENCH_CSR )
    const char* tVariant = "csr";
#elif defined( BENCH_MOVE_NOT_NOEXCEPT )
    const char* tVariant = "move_throwing";
#elif defined( BENCH_PATCHED )
    const char* tVariant = "move";
#else
    const char* tVariant = "copy";
#endif
    const bool tNothrowMove = std::is_nothrow_move_constructible< moris::Vector< int > >::value;

    std::vector< Bench_Vertex > tVertices;
    std::vector< Bench_Cell >   tCells;
    build_tet_mesh( tN, tVertices, tCells );

    for ( int r = 0; r < tRepeats; r++ )
    {
        gNumAllocs  = 0;
        gLiveBytes  = 0;
        gPeakBytes  = 0;
        gCountAlloc = true;
        Result tResult = run_variant( tVertices, tCells );
        gCountAlloc    = false;
        std::printf( "%s,%d,%zu,%zu,%.6f,%.6f,%.6f,%zu,%.1f,%.1f,%d,%zu\n",
                tVariant,
                tN,
                tCells.size(),
                tResult.mNumEdges,
                tResult.mBuild,
                tResult.mShrink,
                tResult.mBuild + tResult.mShrink,
                gNumAllocs,
                gPeakBytes / 1048576.0,
                max_rss_mb(),
                tNothrowMove ? 1 : 0,
                tResult.mChecksum );
        std::fflush( stdout );
    }
    return 0;
}
