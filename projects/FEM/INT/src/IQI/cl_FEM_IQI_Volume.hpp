/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Volume.hpp
 *
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_VOLUME_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_VOLUME_HPP_

#include "moris_typedefs.hpp"                     //MRS/COR/src
#include "cl_Vector.hpp"                          //MRS/CNT/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_IQI.hpp"                   //FEM/INT/src

namespace moris::fem
{
    //------------------------------------------------------------------------------

    class IQI_Volume : public IQI
    {
        enum class IQI_Property_Type
        {
            DENSITY,
            MAX_ENUM
        };

        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------
        /*
         *  constructor
         */
        IQI_Volume();

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IQI_Volume() override{};

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
        void compute_QI( Matrix< DDRMat >& aQI ) override;

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
                Vector< MSI::Dof_Type >& aDofType,
                Matrix< DDRMat >&        adQIdu ) override;

        //------------------------------------------------------------------------------
        /**
         * analytical geometric sensitivities are available if the density is constant or has space derivative
         * functions (see classify_property_spatial_variation); otherwise the FD geometry path is used
         */
        bool supports_analytical_dQIdp_geometry() override;

        //------------------------------------------------------------------------------
        /**
         * derivatives of the QI contribution wrt the evaluation point location, at fixed weight
         * (see IQI::compute_dQIdx_geometry_integrand), the spatial variation of the density (box bounds held fixed)
         * @param[ in ]  aWStar weight associated to the evaluation point
         * @param[ out ] adQIdx derivatives wrt the evaluation point location ( space dim x 1 )
         * @param[ out ] adQIdn derivatives wrt the normal (empty, bulk IQI)
         */
        void compute_dQIdx_geometry_integrand(
                real              aWStar,
                Matrix< DDRMat >& adQIdx,
                Matrix< DDRMat >& adQIdn ) override;

        //------------------------------------------------------------------------------
        /**
         * @brief checks if the gauss point location lies with the box bounds specified in input(mParameters)
         *
         */
        bool is_within_box_bounds();
    };
}    // namespace moris::fem

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_VOLUME_HPP_ */

