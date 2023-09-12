/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Core.hpp
 *
 */

#ifndef PROJECTS_GEN_SDF_SRC_CL_SDF_CORE_HPP_
#define PROJECTS_GEN_SDF_SRC_CL_SDF_CORE_HPP_

#include "typedefs.hpp" // COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_SDF_Object.hpp"
#include "cl_SDF_Mesh.hpp"
#include "cl_SDF_Parameters.hpp"
#include "cl_SDF_Data.hpp"

namespace moris
{
    namespace sdf
    {
//-------------------------------------------------------------------------------

        class Core
        {
                  Mesh          & mMesh;
                  Data          & mData;

                  uint            mCandidateSearchDepth = 1;
                  real            mCandidateSearchDepthEpsilon = 0.01;
                  bool            mVerbose;

//-------------------------------------------------------------------------------
        public :
//-------------------------------------------------------------------------------

            Core(       Mesh & aMesh,
                        Data & aData,
                        bool   aVerbose=false );

//-------------------------------------------------------------------------------

            ~Core(){};

//-------------------------------------------------------------------------------

            void
            set_candidate_search_depth( const uint aCandidateSearchDepth )
            {
                mCandidateSearchDepth = aCandidateSearchDepth;
            }

//-------------------------------------------------------------------------------

            void
            set_candidate_search_epsilon( const real aCandidateSearchEpsilon )
            {
                mCandidateSearchDepthEpsilon = aCandidateSearchEpsilon;
            }

//-------------------------------------------------------------------------------

            void
            calculate_raycast();

//-------------------------------------------------------------------------------

            void
            calculate_raycast(
                    Matrix< IndexMat > & aElementsAtSurface );

//-------------------------------------------------------------------------------
            void
            calculate_raycast(
                    Matrix< IndexMat > & aElementsAtSurface,
                    Matrix< IndexMat > & aElementsInVolume );

//-------------------------------------------------------------------------------

            void
            calculate_raycast_and_sdf( Matrix< DDRMat> & aSDF );

//-------------------------------------------------------------------------------

            void
            calculate_raycast_and_sdf(
                    Matrix< DDRMat>    & aSDF,
                    Matrix< IndexMat > & aElementsAtSurface,
                    Matrix< IndexMat > & aElementsInVolume );

//-------------------------------------------------------------------------------

            void
            save_to_vtk( const std::string & aFilePath );

//-------------------------------------------------------------------------------

            void
            save_unsure_to_vtk( const std::string & aFilePath );

//-------------------------------------------------------------------------------
        private :
//-------------------------------------------------------------------------------

//-------------------------------------------------------------------------------

            void
            voxelize( const uint aAxis );

//-------------------------------------------------------------------------------

            moris::Cell< Vertex * >
            set_candidate_list(  );

//-------------------------------------------------------------------------------

            void
            calculate_udf(moris::Cell< Vertex * > & aCandidateList);

//-------------------------------------------------------------------------------

            /**
             * Kehrwoche :
             * make sure that each vertex is really associated
             * to its closest triangle
             */
            void
            sweep();

//-------------------------------------------------------------------------------

            void
            fill_sdf_with_values( Matrix< DDRMat > & aSDF );

//-------------------------------------------------------------------------------

            void
            preselect_triangles_x( const Matrix< F31RMat >& aPoint );

//-------------------------------------------------------------------------------

            void
            preselect_triangles_y( const Matrix< F31RMat >& aPoint );

//-------------------------------------------------------------------------------

            void
            preselect_triangles_z( const Matrix< F31RMat >& aPoint );

//-------------------------------------------------------------------------------

            void
            intersect_triangles(
                    const uint aAxis,
                    const Matrix< F31RMat >& aPoint );

//-------------------------------------------------------------------------------

            void
            intersect_ray_with_triangles(
                    const uint aAxis,
                    const Matrix< F31RMat >& aPoint,
                    const uint aNodeIndex );

//-------------------------------------------------------------------------------

            void
            check_if_node_is_inside(
                    const uint aAxis,
                    const uint aNodeIndex );
//-------------------------------------------------------------------------------

            void
            calculate_candidate_points_and_buffer_diagonal();

//-------------------------------------------------------------------------------

            void
            get_nodes_withing_bounding_box_of_triangle(
                            Triangle * aTriangle, moris::Cell< Vertex* > & aNodes,
							moris::Cell< Vertex * > & aCandList );

//-------------------------------------------------------------------------------

            /**
             * performs a floodfill in order to fix unsure nodes
             */
            void
            force_unsure_nodes_outside();

//-------------------------------------------------------------------------------

            void
            random_rotation();

//-------------------------------------------------------------------------------

            void
            undo_rotation();

        };

//-------------------------------------------------------------------------------
    } /* namespace sdf */
} /* namespace moris */

#endif /* PROJECTS_GEN_SDF_SRC_CL_SDF_CORE_HPP_ */

