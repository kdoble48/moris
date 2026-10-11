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

#include "cl_FEM_IQI.hpp"                   //FEM/INT/src

namespace moris::fem
{
    //------------------------------------------------------------------------------

    class IQI_Max_Stress : public IQI
    {

        //------------------------------------------------------------------------------

      private:
        //------------------------------------------------------------------------------

        // stress type to evaluate
        enum Stress_Type mStressType;

        enum class IQI_Property_Type
        {
            REFERENCE_VALUE,
            EXPONENT,
            SHIFT,
            KS_PARAMETER,
            MAX_ENUM
        };

        enum class IQI_Constitutive_Type
        {
            ELAST_LIN_ISO,
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
            MORIS_ERROR( false, "IQI_Max_Stress::compute_dQIdu() - not implemented for a drag/lift coefficient IQI." );
        }

        //------------------------------------------------------------------------------
        /**
         * analytical geometric sensitivities are available for the von Mises stress if the elasticity CM is a
         * linear structural CM with exact spatial strain derivatives (see CM_Struc_Linear::has_spatial_strain_derivatives)
         */
        bool supports_analytical_dQIdp_geometry() override;

        //------------------------------------------------------------------------------
        /**
         * derivatives of the QI contribution wrt the evaluation point location, at fixed weight
         * (see IQI::compute_dQIdx_geometry_integrand); spatial variations of the properties are not differentiated
         * @param[ in ]  aWStar weight associated to the evaluation point
         * @param[ out ] adQIdx derivatives wrt the evaluation point location ( space dim x 1 )
         * @param[ out ] adQIdn derivatives wrt the normal (empty, bulk IQI)
         */
        void compute_dQIdx_geometry_integrand(
                real              aWStar,
                Matrix< DDRMat > &adQIdx,
                Matrix< DDRMat > &adQIdn ) override;

        //------------------------------------------------------------------------------
        /**
         * gets the stress vector from the constitutive model and sorts it into a format standardized
         * (independent of 2D-plane strain, 2D-plain stress, and 3D)
         * @param[ out ] tStressVector vector with 6 entries containing the normal and shear stress values
         */
        Matrix< DDRMat > get_stress_vector();

        //------------------------------------------------------------------------------
        /**
         * sorts a vector in the layout of the constitutive model flux into the standardized 6-component layout
         * of get_stress_vector() (independent of 2D-plane strain, 2D-plain stress, and 3D)
         * @param[ in ] aCMStress vector in the layout of the constitutive model flux
         * @return vector with 6 entries containing the normal and shear components
         */
        static Matrix< DDRMat > standardize_stress_vector( const Matrix< DDRMat > &aCMStress );

        //------------------------------------------------------------------------------
        /**
         * evaluate the selected stress measure (von Mises, principal, normal or shear)
         */
        real eval_stress_value();

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
         * derivative of the von Mises stress wrt the CM flux components (CM flux ordering)
         * @param[ in ] aVonMises the von Mises stress at the evaluation point
         */
        Matrix< DDRMat > eval_dVonMisesdFlux( real aVonMises );

        //------------------------------------------------------------------------------
        /**
         * evaluate the Von-Mises stress using the stress vector provided by the constitutive model
         * @param[ out ] tStressValue the value of the Von-Mises stress
         */
        real eval_Von_Mises_stress();

        //------------------------------------------------------------------------------
        /**
         * evaluate the i-th principal stress using the stress vector provided by the constitutive model
         * @param[ in ] aPrincipalStressIndex number i of the principal stress to get,
         *               i.e. 1,2,3 for the first, second and third principal stress, respectively
         * @param[ out ] tStressValue the value of the requested principal stress
         */
        real eval_principal_stress( uint aPrincipalStressIndex );

        //------------------------------------------------------------------------------
        /**
         * evaluate the i-th normal stress using the stress vector provided by the constitutive model
         * @param[ in ] aStressIndex number of the principal stress to get,
         *               i.e. 1,2,3 for the x-, y-, and z-normal stresses, respectively
         * @param[ out ] tStressValue the value of the requested principal stress
         */
        real eval_normal_stress( uint aStressIndex );

        //------------------------------------------------------------------------------
        /**
         * evaluate the i-th shear stress using the stress vector provided by the constitutive model
         * @param[ in ] aStressIndex number of the principal stress to get,
         *               i.e. 1,2,3 for the yz-, xz-, and xy-shear stresses, respectively
         *               in 2D the input is ignored, as there's only one shear stress
         * @param[ out ] tStressValue the value of the requested principal stress
         */
        real eval_shear_stress( uint aStressIndex );

        //------------------------------------------------------------------------------
    };
}    // namespace moris::fem

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_MAX_STRESS_HPP_ */

