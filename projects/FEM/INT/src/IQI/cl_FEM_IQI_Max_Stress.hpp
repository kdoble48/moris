/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Max_Stress.hpp
 *
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_MAX_STRESS_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_MAX_STRESS_HPP_

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
     * integrand of a P-norm or Kreisselmeier-Steinhauser aggregate of a scalar stress measure;
     * the stress measures and their derivatives wrt the CM flux are inherited from IQI_Stress
     */
    class IQI_Max_Stress : public IQI_Stress
    {

        //------------------------------------------------------------------------------

      private:
        //------------------------------------------------------------------------------

        enum class IQI_Property_Type
        {
            REFERENCE_VALUE,
            EXPONENT,
            SHIFT,
            KS_PARAMETER,
            MAX_ENUM
        };

        //------------------------------------------------------------------------------

      protected:
        //------------------------------------------------------------------------------
        /*
         * default constructor
         */
        IQI_Max_Stress()
        {
            MORIS_ERROR( false, " IQI_Max_Stress: Default constructor unavailable. Use enum to construct IQI. " );
        }

        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------
        /*
         * constructor
         * @param[ in ] aStressType scalar stress measure to aggregate (STRESS_VECTOR is not allowed)
         */
        IQI_Max_Stress( enum Stress_Type aStressType );

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IQI_Max_Stress() override{};

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
         * compute the derivative of the quantity of interest wrt dof types (von Mises stress only)
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void compute_dQIdu( real aWStar ) override;

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
            MORIS_ERROR( false, "IQI_Max_Stress::compute_dQIdu() - not implemented." );
        }

        //------------------------------------------------------------------------------
        /**
         * read the aggregation properties: ReferenceValue (required), Exponent (required for the P-norm form),
         * Shift (default 1.0) and KSParameter (optional; > 0 selects the Kreisselmeier-Steinhauser form)
         */
        void get_aggregation_parameters( real& aRefValue, real& aExponent, real& aShift, real& aKSParameter );

        //------------------------------------------------------------------------------
        /**
         * pointwise integrand of the aggregate:
         *   P-norm (default):        ( sigma / sigma_ref - shift )^p
         *   KS (KSParameter = rho):  exp( rho ( sigma / sigma_ref - shift ) )
         * The global P-norm / KS value is assembled from the integral in the input deck.
         */
        real eval_integrand( real aStressValue );

        //------------------------------------------------------------------------------
        /**
         * derivative of the integrand wrt the stress measure
         */
        real eval_dintegrand_dstress( real aStressValue );

        //------------------------------------------------------------------------------
        /**
         * analytical geometric sensitivities are available for the von Mises stress (P-norm and KS integrands)
         * if the elasticity CM is a linear structural CM with exact spatial strain derivatives
         * (see CM_Struc_Linear::has_spatial_strain_derivatives)
         */
        bool supports_analytical_dQIdp_geometry() override;

        //------------------------------------------------------------------------------
        /**
         * derivatives of the QI contribution wrt the evaluation point location, at fixed weight
         * (see IQI::compute_dQIdx_geometry_integrand); the properties are held fixed in space
         * @param[ in ]  aWStar weight associated to the evaluation point
         * @param[ out ] adQIdx derivatives wrt the evaluation point location ( space dim x 1 )
         * @param[ out ] adQIdn derivatives wrt the normal (empty, bulk IQI)
         */
        void compute_dQIdx_geometry_integrand(
                real              aWStar,
                Matrix< DDRMat >& adQIdx,
                Matrix< DDRMat >& adQIdn ) override;

        //------------------------------------------------------------------------------
    };
}    // namespace moris::fem

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_MAX_STRESS_HPP_ */
