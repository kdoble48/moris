/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Strain_Energy.hpp
 *
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_STRAIN_ENERGY_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_STRAIN_ENERGY_HPP_

#include <map>

#include "moris_typedefs.hpp"                     //MRS/COR/src
#include "cl_Vector.hpp"                          //MRS/CNT/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_IQI.hpp"                   //FEM/INT/src

namespace moris::fem
{
    //------------------------------------------------------------------------------

    class IQI_Strain_Energy : public IQI
    {
        //------------------------------------------------------------------------------

        enum class IWG_Property_Type
        {
            BEDDING,
            THICKNESS,
            MAX_ENUM
        };

        enum class IQI_Constitutive_Type
        {
            ELAST,
            MAX_ENUM
        };

      public:
        //------------------------------------------------------------------------------
        /*
         * constructor
         */
        IQI_Strain_Energy();

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IQI_Strain_Energy() override{};

        //------------------------------------------------------------------------------

      private:
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
        void compute_dQIdu( real aWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of the quantity of interest wrt dof types
         * @param[ in ] aDofType group of dof types wrt which derivatives are evaluated
         * @param[ in ] adQIdu   derivative of quantity of interest matrix to fill
         */
        void compute_dQIdu(
                Vector< MSI::Dof_Type > &aDofType,
                Matrix< DDRMat >        &adQIdu ) override;

        //------------------------------------------------------------------------------
        /**
         * analytical geometric sensitivities are available if the elasticity CM is a linear structural CM
         * with exact spatial strain derivatives (see CM_Struc_Linear::has_spatial_strain_derivatives) and known
         * spatial derivatives of its constitutive matrix (CM_Struc_Linear::has_spatial_constitutive_derivatives),
         * and the bedding and thickness properties are constant or have space derivative functions
         * (see classify_property_spatial_variation); otherwise the FD geometry path is used
         */
        bool supports_analytical_dQIdp_geometry() override;

        //------------------------------------------------------------------------------
        /**
         * the thickness property scales the whole QI contribution
         */
        std::shared_ptr< Property > get_geometric_sensitivity_thickness_property() override;

        //------------------------------------------------------------------------------
        /**
         * derivatives of the QI contribution wrt the evaluation point location, at fixed weight
         * (see IQI::compute_dQIdx_geometry_integrand), including the spatial variation of the constitutive matrix
         * and the bedding (the thickness is added by IQI::compute_dQIdp_geometry_analytical)
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
         * @brief checks if the gauss point location lies with the box bounds specified in input(mParameters)
         *
         */
        bool is_within_box_bounds();
    };
}    // namespace moris::fem

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_STRAIN_ENERGY_HPP_ */

