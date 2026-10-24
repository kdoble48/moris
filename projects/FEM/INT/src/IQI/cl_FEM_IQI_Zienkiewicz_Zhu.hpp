/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Zienkiewicz_Zhu.hpp
 *
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_ZIENKIEWICZ_ZHU_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_ZIENKIEWICZ_ZHU_HPP_

#include <map>

#include "moris_typedefs.hpp"                     //MRS/COR/src
#include "cl_Vector.hpp"                          //MRS/CNT/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_IQI_Stress.hpp"            //FEM/INT/src

namespace moris::fem
{
    //------------------------------------------------------------------------------
    /**
     * Zienkiewicz-Zhu error indicator of a scalar stress measure s (von Mises, a principal, normal or shear stress;
     * Zienkiewicz & Zhu, IJNME 24:337-357, 1987).
     *
     * Integrand at an evaluation point:
     *
     *     e^2 = ( s* - s( sigma( u_h ) ) )^2
     *
     * - s( sigma( u_h ) ): the stress measure of the discrete solution, from the leader CM "ElastLinIso"
     *   ( sigma = C eps( u_h ), discontinuous across elements );
     * - s*: the recovered (smoothed) stress measure, the value of the dof field named in "dof_quantity"
     *   (e.g. STRESS_DOF). It is computed by a separate IWG, typically IWG_Struc_Stress (STRUC_VON_MISES_STRESS, ...),
     *   as the L2 projection of s onto the s* interpolation space:
     *
     *         int_Omega v ( s* - s( sigma( u_h ) ) ) dOmega = 0   for all test functions v
     *
     * Integrated, the IQI gives the squared estimate of the stress error in the L2 norm:
     *
     *     eta^2   = int_Omega ( s* - s( sigma( u_h ) ) )^2 dOmega   (GLOBAL)
     *     eta_K^2 = int_K ( s* - s( sigma( u_h ) ) )^2 dK per IG cell (ELEMENTAL_INT; ELEMENTAL_AVG divides by |K|)
     *
     * Take the square root for the error itself. No dQI/du (not usable as an optimization criterion).
     * The stress measure s is evaluated by IQI_Stress (see fn_FEM_Stress_Measures).
     */
    class IQI_Zienkiewicz_Zhu : public IQI_Stress
    {

        //------------------------------------------------------------------------------

      protected:
        //------------------------------------------------------------------------------
        /*
         * default constructor
         */
        IQI_Zienkiewicz_Zhu()
        {
            MORIS_ERROR( false, " IQI_Zienkiewicz_Zhu: Default constructor unavailable. Use enum to construct IQI. " );
        };

        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------
        /*
         * constructor
         * @param[ in ] aStressType scalar stress measure s (STRESS_VECTOR is not allowed)
         */
        IQI_Zienkiewicz_Zhu( enum Stress_Type aStressType );

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IQI_Zienkiewicz_Zhu() override{};

        //------------------------------------------------------------------------------

      private:
        //------------------------------------------------------------------------------
        /**
         * pointwise error e = s* - s( sigma( u_h ) ) of the recovered stress measure
         */
        real eval_error();

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
            MORIS_ERROR( false, "IQI_Zienkiewicz_Zhu::compute_dQIdu - not implemented." );
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
            MORIS_ERROR( false, "IQI_Zienkiewicz_Zhu::compute_dQIdu() - not implemented." );
        }

        //------------------------------------------------------------------------------
    };
}    // namespace moris::fem

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_ZIENKIEWICZ_ZHU_HPP_ */

