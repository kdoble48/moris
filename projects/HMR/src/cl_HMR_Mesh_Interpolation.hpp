/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Mesh_Interpolation.hpp
 *
 */

#pragma once

#include <utility>

#include "cl_HMR_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"

namespace moris::hmr
{
    // alignas: mtk::Mesh is 16-byte aligned (armadillo Matrix members), but as a virtual base it does not
    // align this class's non-virtual part; GCC 13 nevertheless emits 16-byte aligned stores in the
    // base-object constructor, which faults when the class sits at an 8-mod-16 offset in a derived class.
    class alignas( alignof( Mesh ) ) Interpolation_Mesh_HMR : public virtual Mesh
            , public mtk::Interpolation_Mesh
    {
        //-------------------------------------------------------------------------------
      public:
        //-------------------------------------------------------------------------------

        Interpolation_Mesh_HMR(
                std::shared_ptr< Database > aDatabase,
                uint                        aLagrangeMeshIndex )
                : Mesh(
                          std::move( aDatabase ),
                          aLagrangeMeshIndex )
        {
        }

        //-------------------------------------------------------------------------------

        Interpolation_Mesh_HMR(
                std::shared_ptr< Database > aDatabase,
                uint                        aLagrangeOrder,
                uint                        aLagrangePattern )
                : Mesh(
                          std::move( aDatabase ),
                          aLagrangeOrder,
                          aLagrangePattern )
        {
        }

        //-------------------------------------------------------------------------------

        Interpolation_Mesh_HMR(
                std::shared_ptr< Database > aDatabase,
                uint                        aOrder,
                uint                        aLagrangePattern,
                BSpline_Mesh_Base*          aDummyBSplineMesh )
                : Mesh(
                        std::move( aDatabase ),
                        aOrder,
                        aLagrangePattern,
                        aDummyBSplineMesh )
        {
        }

        //-------------------------------------------------------------------------------

        Interpolation_Mesh_HMR(
                std::shared_ptr< Database > aDatabase,
                uint                        aLagrangeOrder,
                uint                        aLagrangePattern,
                uint                        aBSplineOrder,
                uint                        aBsplinePattern )
                : Mesh(
                          std::move( aDatabase ),
                          aLagrangeOrder,
                          aLagrangePattern,
                          aBSplineOrder,
                          aBsplinePattern )
        {
        }

        //-------------------------------------------------------------------------------

        ~Interpolation_Mesh_HMR() override
        {
        }

        //-------------------------------------------------------------------------------

    };    // end class: Interpolation_Mesh_HMR

    //-------------------------------------------------------------------------------

}    // namespace moris::hmr
