/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Stress.hpp
 *
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_STRESS_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_STRESS_HPP_

#include <map>

#include "moris_typedefs.hpp"                     //MRS/COR/src
#include "cl_Vector.hpp"                          //MRS/CNT/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_IQI.hpp"                   //FEM/INT/src

namespace moris::fem
{
    //------------------------------------------------------------------------------

    /**
     * pointwise stress measures (von Mises, principal, normal, shear, or the full stress vector) of an elasticity CM;
     * derived IQIs (e.g. IQI_Max_Stress) reuse the scalar measures and their derivatives wrt the CM flux.
     * The measures themselves are implemented once in fn_FEM_Stress_Measures.
     */
    class IQI_Stress : public IQI
    {

        //------------------------------------------------------------------------------

      protected:
        //------------------------------------------------------------------------------

        // stress type to evaluate
        enum Stress_Type mStressType;

        // flux type to evaluate
        enum CM_Function_Type mFluxType;

        enum class IQI_Constitutive_Type
        {
            ELAST_LIN_ISO,
            MAX_ENUM
        };

        //------------------------------------------------------------------------------
        /*
         * default constructor
         */
        IQI_Stress()
        {
            MORIS_ERROR( false, " IQI_Stress: Default constructor unavailable. Use enum to construct IQI. " );
        };

        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------
        /*
         * constructor
         * @param[ in ] aStressType stress type to evaluate stress in IQI
         * @param[ in ] aFluxType   flux type to evaluate for the CM to evaluate stress in IQI
         */
        IQI_Stress(
                enum Stress_Type      aStressType,
                enum CM_Function_Type aFluxType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IQI_Stress() override{};

        //------------------------------------------------------------------------------

      protected:
        //------------------------------------------------------------------------------
        /**
         * gets the stress vector from the constitutive model and sorts it into a format standardized
         * (independent of 2D-plane strain, 2D-plain stress, and 3D)
         * @param[ out ] aStressVector vector with 6 entries ( s_xx, s_yy, s_zz, s_yz, s_xz, s_xy )
         */
        void get_stress_vector( Matrix< DDRMat > &aStressVector );

        //------------------------------------------------------------------------------
        /**
         * evaluate the selected scalar stress measure (von Mises, principal, normal or shear)
         * @return the value of the stress measure
         */
        real eval_stress_value();

        //------------------------------------------------------------------------------
        /**
         * derivative of the selected scalar stress measure wrt the CM flux components (von Mises stress only)
         * @return derivative in the layout of the CM flux ( 3, 4 or 6 x 1 )
         */
        Matrix< DDRMat > eval_dStressdFlux();

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
        void compute_QI( Matrix< DDRMat > &aQI ) override;

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of the quantity of interest wrt dof types
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void compute_dQIdu( real aWStar ) override
        {
            MORIS_ERROR( false, "IQI_Stress::compute_dQIdu - not implemented." );
        }

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of the quantity of interest wrt dof types
         * @param[ in ] aDofType group of dof types wrt which derivatives are evaluated
         * @param[ in ] adQIdu   derivative of quantity of interest matrix to fill
         */
        void compute_dQIdu(
                Vector< MSI::Dof_Type > &aDofType,
                Matrix< DDRMat >        &adQIdu ) override
        {
            MORIS_ERROR( false, "IQI_Stress::compute_dQIdu() - not implemented." );
        }

        //------------------------------------------------------------------------------
        /**
         * compute matrix dimension of the IQI
         * @param[ out ] space dimension of the IQI
         */
        std::pair< uint, uint > get_matrix_dim() override;
    };
}    // namespace moris::fem

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_STRESS_HPP_ */
