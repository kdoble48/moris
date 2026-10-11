/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integrator.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTEGRATOR_HPP_
#define SRC_MTK_CL_MTK_INTEGRATOR_HPP_

#include "moris_typedefs.hpp"               //MRS/COR/src
#include "cl_Matrix.hpp"                    //LNA/src
#include "cl_MTK_Enums.hpp"                 //MTK/src
#include "cl_MTK_Integration_Rule.hpp"      //MTK/src
#include "cl_MTK_Integration_Coeffs.hpp"    //MTK/src

namespace moris::mtk
{
    //------------------------------------------------------------------------------

    class Integrator
    {
        // pointer to space rule, if specified
        std::unique_ptr< Integration_Coeffs_Base > mSpaceCoeffs;

        // pointer to time rule, if specified
        std::unique_ptr< Integration_Coeffs_Base > mTimeCoeffs;

        // number of points in space
        uint mNumOfSpacePoints;

        // number of points in time
        uint mNumOfTimePoints;

        // matrix with space points
        Matrix< DDRMat > mSpacePoints;

        // matrix with time points
        Matrix< DDRMat > mTimePoints;

        // matrix with space weights
        Matrix< DDRMat > mSpaceWeights;

        // matrix with time weights
        Matrix< DDRMat > mTimeWeights;

        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------

        /**
         * constructs an integrator from an integration rule
         **/
        Integrator( const Integration_Rule &aIntegrationRule );

        //------------------------------------------------------------------------------
        /**
         * get the number of integration points
         **/
        uint get_number_of_points() const;

        //------------------------------------------------------------------------------
        /**
         * @brief get the integration points as a (d x n) matrix, where d is the dimension of the space
         * and n is the number of integration points (each column contains one point)
         * @param aIntegrationPoints
         */
        void get_points( Matrix< DDRMat > &aIntegrationPoints ) const;

        //------------------------------------------------------------------------------
        /**
         * @brief get the integration points as a (d x n) matrix, where d is the dimension of the space
         * and n is the number of integration points (each column contains one point)
         * @return
         */
        Matrix< DDRMat > get_points() const;

        //------------------------------------------------------------------------------
        /**
         * get the integration point weights
         **/
        void get_weights( Matrix< DDRMat > &aIntegrationWeights ) const;

        Matrix< DDRMat > get_weights() const;

        //------------------------------------------------------------------------------
        /**
         * get the space-only integration points as a (d x n) matrix
         **/
        const Matrix< DDRMat >&
        get_space_points() const
        {
            return mSpacePoints;
        }

        //------------------------------------------------------------------------------
        /**
         * get the space-only integration weights as a (1 x n) matrix
         **/
        const Matrix< DDRMat >&
        get_space_weights() const
        {
            return mSpaceWeights;
        }

        //------------------------------------------------------------------------------
        /**
         * get the time-only integration points as a (1 x n) matrix
         **/
        const Matrix< DDRMat >&
        get_time_points() const
        {
            return mTimePoints;
        }

        //------------------------------------------------------------------------------
        /**
         * get the time-only integration weights as a (1 x n) matrix
         **/
        const Matrix< DDRMat >&
        get_time_weights() const
        {
            return mTimeWeights;
        }

        //------------------------------------------------------------------------------

    };    // class Integrator

    //------------------------------------------------------------------------------

}    // namespace moris::mtk

#endif /* SRC_MTK_CL_MTK_INTEGRATOR_HPP_ */
