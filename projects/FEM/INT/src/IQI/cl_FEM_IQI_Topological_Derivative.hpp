/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Topological_Derivative.hpp
 *
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_TOPOLOGICAL_DERIVATIVE_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_TOPOLOGICAL_DERIVATIVE_HPP_

#include "moris_typedefs.hpp"    // MRS/COR/src
#include "cl_Vector.hpp"         // MRS/CNT/src

#include "cl_Matrix.hpp"          // LINALG/src
#include "linalg_typedefs.hpp"    // LINALG/src

#include "cl_FEM_IQI.hpp"    // FEM/INT/src

namespace moris::fem
{
    //------------------------------------------------------------------------------

    /**
     * @brief Topological derivative of the compliance (2D linear elasticity, circular inclusion) in the
     * local-stress form of Amstutz and Andrae (2006), eq. 9, with the adjoint v = -u of the compliance:
     *
     *   dJ = -(r-1)/(chi r+1) (chi+1)/2 [ 2 sigma:eps + (r-1)(chi-2)/(chi+2r-1) tr(sigma) tr(eps) ] + delta l
     *
     * sigma is the phase-local stress (the CM flux, which already carries the phase's Young's modulus) and eps the
     * strain of the leader CM "Elast"; both traces are in-plane (sigma_xx + sigma_yy, eps_xx + eps_yy).
     * r is the contrast (stiffness of the inclusion over the stiffness at the point), chi the polarisation constant
     * (lambda + 3 mu)/(lambda + mu), delta = -1 for a hole in the material and +1 for material added in the void.
     * At chi = 2 the trace term vanishes and dJ = -3(r-1)/(2r+1) sigma:eps + delta l (Teixeira et al. 2025, eq. 15).
     *
     * The IQI returns the generalised topological derivative g = s delta dJ, i.e. g = -dJ on the material set
     * (delta = -1) and g = +dJ on the void set (delta = +1) for s = +1 (Teixeira et al. 2025, eq. 16).
     *
     * Parameters ("function_parameters", separated by "/"):
     *   0: r      contrast ratio, r >= 0 (material set: alpha_out/alpha_in; void set: alpha_in/alpha_out)
     *   1: delta  -1 (material set, hole) or +1 (void set, inclusion)
     *   2: l      area weight (fixed Lagrange multiplier), default 0
     *   3: chi    polarisation constant override; 0 (default) takes chi from the CM: 3 - 4 nu for PLANE_STRAIN,
     *             (3 - nu)/(1 + nu) for PLANE_STRESS, nu from the CM property "PoissonRatio"; any positive value
     *             is used as is (2 reproduces Teixeira et al. 2025, eq. 15)
     *   4: s      output sign, +1 (default) or -1 (flipped level-set convention)
     *
     * 2D only (PLANE_STRESS or PLANE_STRAIN): no 3D or axisymmetric formula is implemented.
     * The topological derivative is not differentiated: compute_dQIdu errors (use the finite-difference dQIdu if a
     * set ever needs it).
     */
    class IQI_Topological_Derivative : public IQI
    {
        //------------------------------------------------------------------------------

        enum class IQI_Constitutive_Type
        {
            ELAST,
            MAX_ENUM
        };

        // parameters, set by set_parameters
        bool mParametersSet = false;
        real mContrast      = 0.0;
        real mDelta         = -1.0;
        real mAreaWeight    = 0.0;
        real mChiOverride   = 0.0;
        real mOutputSign    = 1.0;

      public:
        //------------------------------------------------------------------------------
        /*
         * constructor
         */
        IQI_Topological_Derivative();

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IQI_Topological_Derivative() override {};

        //------------------------------------------------------------------------------
        /**
         * set and check the parameters ( r / delta / l / chi / s )
         * @param[ in ] aParameters a list of parameters
         */
        void set_parameters( const Vector< Matrix< DDRMat > >& aParameters ) override;

        //------------------------------------------------------------------------------

      private:
        //------------------------------------------------------------------------------
        /**
         * evaluate the generalised topological derivative at the evaluation point
         */
        real evaluate_topological_derivative();

        //------------------------------------------------------------------------------
        /**
         * polarisation constant chi from the override or from the CM model type and Poisson ratio
         * @param[ in ] aModelType plane model type of the elasticity CM
         */
        real get_chi( Model_Type aModelType );

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
         * not implemented: the topological derivative is not differentiated
         */
        void compute_dQIdu( real aWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * not implemented: the topological derivative is not differentiated
         */
        void compute_dQIdu(
                Vector< MSI::Dof_Type >& aDofType,
                Matrix< DDRMat >&        adQIdu ) override;

        //------------------------------------------------------------------------------
    };
}    // namespace moris::fem

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_TOPOLOGICAL_DERIVATIVE_HPP_ */
