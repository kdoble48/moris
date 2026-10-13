/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Cylindrical_Stress.hpp
 *
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_CYLINDRICAL_STRESS_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_CYLINDRICAL_STRESS_HPP_

#include <map>

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"         //MRS/CNT/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_IQI.hpp"    //FEM/INT/src

namespace moris::fem
{
    //------------------------------------------------------------------------------

    /**
     * IQI evaluating the in-plane cylindrical (polar) stress components about a
     * user-specified center point ( cx, cy ). The Cartesian stress is obtained from
     * the leader linear-isotropic-elasticity constitutive model and rotated into the
     * polar frame defined by the angle theta = atan2( y - cy, x - cx ).
     *
     * The component returned is selected via mIQITypeIndex:
     *   0 -> sigma_rr      (radial)
     *   1 -> sigma_thetatheta (hoop)
     *   2 -> sigma_rtheta  (shear)
     */
    class IQI_Cylindrical_Stress : public IQI
    {

        //------------------------------------------------------------------------------

      private:
        //------------------------------------------------------------------------------

        // flux type to evaluate
        enum CM_Function_Type mFluxType;

        enum class IQI_Constitutive_Type
        {
            ELAST_LIN_ISO,
            MAX_ENUM
        };

        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------
        /*
         * constructor
         * @param[ in ] aFluxType flux type to evaluate for the CM to evaluate stress in IQI
         */
        explicit IQI_Cylindrical_Stress(
                enum CM_Function_Type aFluxType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IQI_Cylindrical_Stress() override {};

        //------------------------------------------------------------------------------

      private:
        //------------------------------------------------------------------------------
        /**
         * gets the stress vector from the constitutive model and sorts it into a
         * standardized format (independent of 2D-plane strain, 2D-plane stress, and 3D)
         * @param[ out ] aStressVector vector with 6 entries containing the normal and shear stress values
         */
        void get_stress_vector( Matrix< DDRMat >& aStressVector );

        //------------------------------------------------------------------------------
        /**
         * evaluate the requested cylindrical/polar stress component about the center
         * point provided through mParameters, using the Cartesian stress vector
         * provided by the constitutive model and the evaluation-point location.
         * @param[ out ] aStressValue the value of the requested cylindrical stress component
         */
        void eval_cylindrical_stress( Matrix< DDRMat >& aStressValue );

        //------------------------------------------------------------------------------
        /**
         * compute the quantity of interest
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void compute_QI( real aWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * Evaluate the quantity of interest and fill aQI with value
         * @param[ in ] aQI IQI value at evaluation point
         */
        void compute_QI( Matrix< DDRMat >& aQI ) override;

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of the quantity of interest wrt dof types
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void compute_dQIdu( real aWStar ) override
        {
            MORIS_ERROR( false, "IQI_Cylindrical_Stress::compute_dQIdu - not implemented." );
        }

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of the quantity of interest wrt dof types
         * @param[ in ] aDofType group of dof types wrt which derivatives are evaluated
         * @param[ in ] adQIdu   derivative of quantity of interest matrix to fill
         */
        void compute_dQIdu(
                Vector< MSI::Dof_Type >& aDofType,
                Matrix< DDRMat >&        adQIdu ) override
        {
            MORIS_ERROR( false, "IQI_Cylindrical_Stress::compute_dQIdu() - not implemented." );
        }

        //------------------------------------------------------------------------------
        /**
         * compute matrix dimension of the IQI
         * @param[ out ] space dimension of the IQI
         */
        std::pair< uint, uint > get_matrix_dim() override;
    };
}    // namespace moris::fem

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_CYLINDRICAL_STRESS_HPP_ */
