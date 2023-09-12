/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Data.cpp
 *
 */

#include "cl_SDF_Data.hpp"
namespace moris
{
    namespace sdf
    {
//-------------------------------------------------------------------------------
        Data::Data( Object & aObject ) :
                  mTriangles( aObject.get_triangles() ),
                  mVertices( aObject.get_vertices() ),
                  mNumberOfTriangles( mTriangles.size() ),
                  mTriangleMinCoordsX(mNumberOfTriangles, 1),
                  mTriangleMinCoordsY(mNumberOfTriangles, 1),
                  mTriangleMinCoordsZ(mNumberOfTriangles, 1),
                  mTriangleMaxCoordsX(mNumberOfTriangles, 1),
                  mTriangleMaxCoordsY(mNumberOfTriangles, 1),
                  mTriangleMaxCoordsZ(mNumberOfTriangles, 1),
#ifdef MORIS_USE_ARMA
                       mCandI(mNumberOfTriangles),
                       mCandJ(mNumberOfTriangles),
                       mCandK(mNumberOfTriangles),
#else
                       mCandJ(mNumberOfTriangles, 1),
#endif
                       mCandidateTriangles(mNumberOfTriangles, 1)

        {
            this->init_triangles();
        }

//-------------------------------------------------------------------------------

        void
        Data::init_triangles()
        {
            // copy triangle bounding box data
            for ( uint k = 0; k < mNumberOfTriangles; ++k)
            {
                // minimum triangle coordinates for lower left point of bounding box
                mTriangleMinCoordsX( k )
                    = mTriangles( k )->get_min_coord( 0 );

                mTriangleMinCoordsY( k )
                    = mTriangles( k )->get_min_coord( 1 );

                mTriangleMinCoordsZ( k )
                    = mTriangles( k )->get_min_coord( 2 );

                // maximum triangle coordinates for upper right point of bounding box
                mTriangleMaxCoordsX( k )
                    = mTriangles( k )->get_max_coord( 0 );

                mTriangleMaxCoordsY( k )
                    = mTriangles( k )->get_max_coord( 1 );

                mTriangleMaxCoordsZ ( k )
                    = mTriangles( k )->get_max_coord( 2 );
            }
        }

//-------------------------------------------------------------------------------
    } /* namespace sdf */
} /* namespace moris */

