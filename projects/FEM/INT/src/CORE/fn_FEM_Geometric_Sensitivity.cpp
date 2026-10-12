/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_FEM_Geometric_Sensitivity.cpp
 *
 */

#include "fn_FEM_Geometric_Sensitivity.hpp"

#include "cl_FEM_Geometry_Interpolator.hpp"
#include "cl_FEM_Field_Interpolator.hpp"

#include "fn_inv.hpp"
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_dot.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    namespace
    {
        /**
         * cross product of rows aRowA and aRowB of the matrices aA and aB (3D), as a column vector
         */
        Matrix< DDRMat >
        cross_rows(
                const Matrix< DDRMat >& aA,
                uint                    aRowA,
                const Matrix< DDRMat >& aB,
                uint                    aRowB )
        {
            return { { aA( aRowA, 1 ) * aB( aRowB, 2 ) - aA( aRowA, 2 ) * aB( aRowB, 1 ) },
                { aA( aRowA, 2 ) * aB( aRowB, 0 ) - aA( aRowA, 0 ) * aB( aRowB, 2 ) },
                { aA( aRowA, 0 ) * aB( aRowB, 1 ) - aA( aRowA, 1 ) * aB( aRowB, 0 ) } };
        }

        //------------------------------------------------------------------------------

        /**
         * physical space derivatives of the bulk IG shape functions ( number of space dimensions x number of vertices )
         * from the IG space Jacobian J = dN/dxi * X
         */
        void
        eval_bulk_ig_dNdx(
                const Matrix< DDRMat >& adNdXi,
                const Matrix< DDRMat >& aJ,
                Matrix< DDRMat >&       adNdx )
        {
            const uint tNumParam = aJ.n_rows();
            const uint tNumSpace = aJ.n_cols();

            if ( tNumParam == tNumSpace )
            {
                // standard parametrization: dN/dx = J^-1 dN/dxi
                adNdx = inv( aJ ) * adNdXi;
            }
            else
            {
                MORIS_ERROR( tNumParam == tNumSpace + 1,
                        "fn_FEM_Geometric_Sensitivity - unsupported bulk IG parametrization ( %d param vs %d space dimensions )",
                        tNumParam,
                        tNumSpace );

                // simplex area/volume coordinates with sum( xi ) = 1:
                // [ 1 ; x ] = A xi with A = [ 1 ... 1 ; J^T ], hence dxi/dx = columns 1..d of A^-1
                Matrix< DDRMat > tA( tNumParam, tNumParam, 1.0 );
                for ( uint iSpace = 0; iSpace < tNumSpace; iSpace++ )
                {
                    for ( uint iParam = 0; iParam < tNumParam; iParam++ )
                    {
                        tA( iSpace + 1, iParam ) = aJ( iParam, iSpace );
                    }
                }

                Matrix< DDRMat > tInvA = inv( tA );

                // dN/dx = ( dxi/dx )^T dN/dxi
                Matrix< DDRMat > tdXidx = tInvA( { 0, tNumParam - 1 }, { 1, tNumParam - 1 } );
                adNdx                   = trans( tdXidx ) * adNdXi;
            }
        }

        //------------------------------------------------------------------------------

        /**
         * unnormalized normal of a side and its derivative wrt the side Jacobian
         * 2D side (line): a = rotated tangent;  3D side: a = t1 x t2
         */
        void
        eval_side_area_vector(
                const Matrix< DDRMat >& aJ,
                Matrix< DDRMat >&       aA )
        {
            const uint tNumSpace = aJ.n_cols();

            if ( tNumSpace == 2 )
            {
                aA = { { aJ( 0, 1 ) }, { -aJ( 0, 0 ) } };
            }
            else
            {
                MORIS_ERROR( aJ.n_rows() == 2 && tNumSpace == 3,
                        "fn_FEM_Geometric_Sensitivity - unsupported side IG parametrization" );

                aA = cross_rows( aJ, 0, aJ, 1 );
            }
        }

        //------------------------------------------------------------------------------

        /**
         * derivative of the unnormalized side normal a wrt X( v, i )
         */
        void
        eval_side_area_vector_derivative(
                const Matrix< DDRMat >& adNdXi,
                const Matrix< DDRMat >& aJ,
                uint                    aVertex,
                uint                    aDirection,
                Matrix< DDRMat >&       adA )
        {
            const uint tNumSpace = aJ.n_cols();

            // derivative of the side Jacobian: row r gets dN_v/dxi_r in column i
            Matrix< DDRMat > tdJ( aJ.n_rows(), tNumSpace, 0.0 );
            for ( uint iRow = 0; iRow < aJ.n_rows(); iRow++ )
            {
                tdJ( iRow, aDirection ) = adNdXi( iRow, aVertex );
            }

            if ( tNumSpace == 2 )
            {
                adA = { { tdJ( 0, 1 ) }, { -tdJ( 0, 0 ) } };
            }
            else
            {
                adA = cross_rows( tdJ, 0, aJ, 1 ) + cross_rows( aJ, 0, tdJ, 1 );
            }
        }
    }    // namespace

    //------------------------------------------------------------------------------

    void
    eval_ig_measure_derivatives(
            Geometry_Interpolator* aIGGI,
            bool                   aIsSide,
            Matrix< DDRMat >&      adMeasure )
    {
        // shape function derivatives in the IG parametric space and IG space Jacobian
        const Matrix< DDRMat >& tdNdXi = aIGGI->dNdXi();
        const Matrix< DDRMat >  tJ     = tdNdXi * aIGGI->get_space_coeff();

        const uint tNumVertices = tdNdXi.n_cols();
        const uint tNumSpace    = tJ.n_cols();

        adMeasure.set_size( tNumVertices, tNumSpace, 0.0 );

        if ( !aIsSide )
        {
            // bulk: d detJ / dX( v, i ) = detJ * dN_v / dx_i
            Matrix< DDRMat > tdNdx;
            eval_bulk_ig_dNdx( tdNdXi, tJ, tdNdx );

            adMeasure = trans( tdNdx );
        }
        else
        {
            // side: detJ = c |a|, d detJ / detJ = ( a . da ) / |a|^2
            Matrix< DDRMat > tA;
            eval_side_area_vector( tJ, tA );
            const real tA2 = dot( tA, tA );

            Matrix< DDRMat > tdA;
            for ( uint iVertex = 0; iVertex < tNumVertices; iVertex++ )
            {
                for ( uint iDir = 0; iDir < tNumSpace; iDir++ )
                {
                    eval_side_area_vector_derivative( tdNdXi, tJ, iVertex, iDir, tdA );
                    adMeasure( iVertex, iDir ) = dot( tA, tdA ) / tA2;
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    eval_ig_normal_derivatives(
            Geometry_Interpolator*      aIGGI,
            const Matrix< DDRMat >&     aNormal,
            Vector< Matrix< DDRMat > >& adNormal )
    {
        const Matrix< DDRMat >& tdNdXi = aIGGI->dNdXi();
        const Matrix< DDRMat >  tJ     = tdNdXi * aIGGI->get_space_coeff();

        const uint tNumVertices = tdNdXi.n_cols();
        const uint tNumSpace    = tJ.n_cols();

        // unnormalized normal, oriented as the normal used by the IWG/IQI
        Matrix< DDRMat > tA;
        eval_side_area_vector( tJ, tA );
        const real tSign  = dot( tA, aNormal ) < 0.0 ? -1.0 : 1.0;
        const real tNormA = norm( tA );

        // unit normal consistent with the side geometry and the projector I - n n^T
        Matrix< DDRMat > tN = tSign * tA / tNormA;
        Matrix< DDRMat > tP = -1.0 * tN * trans( tN );
        for ( uint iDim = 0; iDim < tNumSpace; iDim++ )
        {
            tP( iDim, iDim ) += 1.0;
        }

        adNormal.resize( tNumSpace );

        Matrix< DDRMat > tdA;
        for ( uint iDir = 0; iDir < tNumSpace; iDir++ )
        {
            adNormal( iDir ).set_size( tNumSpace, tNumVertices, 0.0 );

            for ( uint iVertex = 0; iVertex < tNumVertices; iVertex++ )
            {
                eval_side_area_vector_derivative( tdNdXi, tJ, iVertex, iDir, tdA );

                // dn = ( I - n n^T ) da / |a|
                Matrix< DDRMat > tdN = ( tSign / tNormA ) * tP * tdA;
                for ( uint iDim = 0; iDim < tNumSpace; iDim++ )
                {
                    adNormal( iDir )( iDim, iVertex ) = tdN( iDim );
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    uint
    second_derivative_row(
            uint aJ,
            uint aK,
            uint aSpaceDim )
    {
        if ( aJ == aK )
        {
            return aJ;
        }

        if ( aSpaceDim == 2 )
        {
            return 2;
        }

        // 3D: xx, yy, zz, yz, xz, xy
        const uint tSum = aJ + aK;
        return tSum == 3 ? 3 : ( tSum == 2 ? 4 : 5 );
    }

    //------------------------------------------------------------------------------

    void
    eval_dNdx_spatial_derivative(
            Field_Interpolator* aFI,
            uint                aDirection,
            Matrix< DDRMat >&   adNdxdxk )
    {
        const uint tSpaceDim = aFI->get_space_dim();

        const Matrix< DDRMat >& td2Ndx2 = aFI->dnNdxn( 2 );

        adNdxdxk.set_size( tSpaceDim, td2Ndx2.n_cols() );

        for ( uint iDim = 0; iDim < tSpaceDim; iDim++ )
        {
            const uint tRow = second_derivative_row( iDim, aDirection, tSpaceDim );

            for ( uint iBase = 0; iBase < td2Ndx2.n_cols(); iBase++ )
            {
                adNdxdxk( iDim, iBase ) = td2Ndx2( tRow, iBase );
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    build_vector_shape_function_matrix(
            const Matrix< DDRMat >& aNRow,
            uint                    aNumFields,
            Matrix< DDRMat >&       aN )
    {
        const uint tNumBases = aNRow.numel();

        aN.set_size( aNumFields, aNumFields * tNumBases, 0.0 );

        for ( uint iField = 0; iField < aNumFields; iField++ )
        {
            for ( uint iBase = 0; iBase < tNumBases; iBase++ )
            {
                aN( iField, iField * tNumBases + iBase ) = aNRow( iBase );
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    vectorize_coefficients(
            const Matrix< DDRMat >& aCoeff,
            Matrix< DDRMat >&       aVec )
    {
        const uint tNumBases  = aCoeff.n_rows();
        const uint tNumFields = aCoeff.n_cols();

        aVec.set_size( tNumBases * tNumFields, 1 );

        for ( uint iField = 0; iField < tNumFields; iField++ )
        {
            for ( uint iBase = 0; iBase < tNumBases; iBase++ )
            {
                aVec( iField * tNumBases + iBase ) = aCoeff( iBase, iField );
            }
        }
    }

    //------------------------------------------------------------------------------

    bool
    has_active_geometry_pdv( const Matrix< DDSMat >& aGeoLocalAssembly )
    {
        for ( uint iEntry = 0; iEntry < aGeoLocalAssembly.numel(); iEntry++ )
        {
            if ( aGeoLocalAssembly( iEntry ) != -1 )
            {
                return true;
            }
        }
        return false;
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem
