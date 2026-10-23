/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Linear_Bulk.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_BULK_HPP_
#define SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_BULK_HPP_

#include <map>

#include "moris_typedefs.hpp"                     //MRS/COR/src
#include "cl_Vector.hpp"                          //MRS/CNT/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_IWG.hpp"                   //FEM/INT/src

namespace moris::fem
{
    //------------------------------------------------------------------------------

    class IWG_Isotropic_Struc_Linear_Bulk : public IWG
    {

        //------------------------------------------------------------------------------

      public:
        enum class IWG_Property_Type
        {
            LOAD,
            BEDDING,
            THICKNESS,
            MAX_ENUM
        };

        enum class IWG_Constitutive_Type
        {
            ELAST_LIN_ISO,
            ELAST_LIN_ISO_PRESSURE,
            MAX_ENUM
        };

        //------------------------------------------------------------------------------
        /*
         * constructor
         */
        IWG_Isotropic_Struc_Linear_Bulk();

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IWG_Isotropic_Struc_Linear_Bulk() override{};

        //------------------------------------------------------------------------------
        /**
         * compute the residual
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void compute_residual( real aWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * compute the jacobian
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void compute_jacobian( real aWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * compute the residual and the jacobian
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void compute_jacobian_and_residual( real aWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of the residual wrt design variables
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void compute_dRdp( real aWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * analytical geometric sensitivities are available if the elasticity CM is a linear structural CM
         * with exact spatial strain derivatives (see CM_Struc_Linear::has_spatial_strain_derivatives) and known
         * spatial derivatives of its constitutive matrix (CM_Struc_Linear::has_spatial_constitutive_derivatives),
         * and the load, bedding and thickness properties are constant or have space derivative functions
         * (see classify_property_spatial_variation); otherwise the FD geometry path is used
         */
        bool supports_analytical_dRdp_geometry() override;

        //------------------------------------------------------------------------------
        /**
         * the thickness property scales the whole residual contribution
         */
        std::shared_ptr< Property > get_geometric_sensitivity_thickness_property() override;

        //------------------------------------------------------------------------------
        /**
         * derivatives of the residual contribution wrt the evaluation point location, at fixed weight
         * (see IWG::compute_dRdx_geometry_integrand), including the spatial variation of the constitutive matrix,
         * the load and the bedding (the thickness is added by IWG::compute_dRdp_geometry_analytical)
         * @param[ in ]  aWStar weight associated to the evaluation point
         * @param[ out ] adRdx  derivatives wrt the evaluation point location, one column per space direction
         * @param[ out ] adRdn  derivatives wrt the normal components (empty for bulk)
         */
        void compute_dRdx_geometry_integrand(
                real                        aWStar,
                Vector< Matrix< DDRMat > >& adRdx,
                Vector< Matrix< DDRMat > >& adRdn ) override;

        //------------------------------------------------------------------------------
    };
    //------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_BULK_HPP_ */

