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
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Property.hpp"

#include "fn_inv.hpp"
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_dot.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    namespace
    {
        // size of the translation probe of classify_property_spatial_variation, relative to the IG cell size
        // (an arbitrary irrational-looking fraction, so that the probe does not hit a period or a symmetry of the function)
        const real tTranslationProbeFraction = 0.1234567;

        //------------------------------------------------------------------------------

        /**
         * true if the two matrices have the same size and bitwise equal entries
         */
        bool
        is_identical(
                const Matrix< DDRMat >& aA,
                const Matrix< DDRMat >& aB )
        {
            if ( aA.n_rows() != aB.n_rows() || aA.n_cols() != aB.n_cols() )
            {
                return false;
            }

            for ( uint iEntry = 0; iEntry < aA.numel(); iEntry++ )
            {
                if ( aA( iEntry ) != aB( iEntry ) )
                {
                    return false;
                }
            }
            return true;
        }

        //------------------------------------------------------------------------------

        /**
         * true if the value function of the property returns the same value when the evaluation point is translated by
         * +/- tTranslationProbeFraction * ( IG cell size ) in each space direction. The IP and IG geometry coefficients are
         * translated together (the parametric point, the Jacobians and the field interpolators are unchanged) and restored.
         * A value function can only see the evaluation point through the field interpolator manager, so a function of its
         * parameters only passes, and a space dependent function fails unless it takes bitwise identical values at all
         * three points of every direction.
         */
        bool
        is_invariant_under_translation( Property& aProperty )
        {
            // get the field interpolator manager and the geometry interpolators
            Field_Interpolator_Manager* tFIManager = aProperty.get_field_interpolator_manager();

            if ( tFIManager == nullptr )
            {
                return false;
            }

            Geometry_Interpolator* tIGGI = tFIManager->get_IG_geometry_interpolator();
            Geometry_Interpolator* tIPGI = tFIManager->get_IP_geometry_interpolator();

            if ( tIGGI == nullptr || tIPGI == nullptr )
            {
                return false;
            }

            // reference value at the evaluation point
            const Matrix< DDRMat > tValue = aProperty.val();

            // unperturbed geometry coefficients ( number of bases x number of space dimensions )
            const Matrix< DDRMat > tXHatIG = tIGGI->get_space_coeff();
            const Matrix< DDRMat > tXHatIP = tIPGI->get_space_coeff();

            const uint tSpaceDim = tXHatIG.n_cols();

            if ( tXHatIP.n_cols() != tSpaceDim )
            {
                return false;
            }

            // size of the IG cell: largest extent of its vertex coordinates
            real tCellSize = 0.0;
            for ( uint iDim = 0; iDim < tSpaceDim; iDim++ )
            {
                real tMin = tXHatIG( 0, iDim );
                real tMax = tXHatIG( 0, iDim );
                for ( uint iVertex = 1; iVertex < tXHatIG.n_rows(); iVertex++ )
                {
                    tMin = std::min( tMin, tXHatIG( iVertex, iDim ) );
                    tMax = std::max( tMax, tXHatIG( iVertex, iDim ) );
                }
                tCellSize = std::max( tCellSize, tMax - tMin );
            }

            if ( !( tCellSize > 0.0 ) )
            {
                return false;
            }

            const real tShift = tTranslationProbeFraction * tCellSize;

            // value function and a copy of its parameters (the function takes them by non-const reference)
            const PropertyFunc&        tValFunction = aProperty.get_val_function();
            Vector< Matrix< DDRMat > > tParameters  = aProperty.get_parameters();

            bool             tIsInvariant = true;
            Matrix< DDRMat > tProbeValue;

            for ( uint iDim = 0; iDim < tSpaceDim && tIsInvariant; iDim++ )
            {
                for ( real tSign : { -1.0, 1.0 } )
                {
                    // translate the IG and IP geometry coefficients
                    Matrix< DDRMat > tXHatIGShifted = tXHatIG;
                    Matrix< DDRMat > tXHatIPShifted = tXHatIP;

                    for ( uint iBase = 0; iBase < tXHatIGShifted.n_rows(); iBase++ )
                    {
                        tXHatIGShifted( iBase, iDim ) += tSign * tShift;
                    }
                    for ( uint iBase = 0; iBase < tXHatIPShifted.n_rows(); iBase++ )
                    {
                        tXHatIPShifted( iBase, iDim ) += tSign * tShift;
                    }

                    tIGGI->set_space_coeff( tXHatIGShifted );
                    tIPGI->set_space_coeff( tXHatIPShifted );

                    // evaluate the value function at the translated point
                    tValFunction( tProbeValue, tParameters, tFIManager );

                    if ( !is_identical( tProbeValue, tValue ) )
                    {
                        tIsInvariant = false;
                        break;
                    }
                }
            }

            // restore the geometry coefficients
            tIGGI->set_space_coeff( tXHatIG );
            tIPGI->set_space_coeff( tXHatIP );

            return tIsInvariant;
        }

        //------------------------------------------------------------------------------

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
    Property_Spatial_Variation
    classify_property_spatial_variation( const std::shared_ptr< Property >& aProperty )
    {
        // no property: nothing varies
        if ( aProperty == nullptr )
        {
            return Property_Spatial_Variation::CONSTANT;
        }

        // a dof, dv or field dependent property varies in space with its field; its total spatial derivative is not provided
        if ( aProperty->get_dof_type_list().size() > 0 ||      //
                aProperty->get_dv_type_list().size() > 0 ||    //
                aProperty->get_field_type_list().size() > 0 )
        {
            return Property_Spatial_Variation::UNKNOWN;
        }

        // the space derivative is provided
        if ( aProperty->check_space_dependency() )
        {
            return Property_Spatial_Variation::DIFFERENTIABLE;
        }

        // default constant value function: the value is the first parameter
        if ( !aProperty->has_val_function() )
        {
            return Property_Spatial_Variation::CONSTANT;
        }

        // user value function: constant if it does not see a translation of the evaluation point
        if ( is_invariant_under_translation( *aProperty ) )
        {
            return Property_Spatial_Variation::CONSTANT;
        }

        return Property_Spatial_Variation::UNKNOWN;
    }

    //------------------------------------------------------------------------------

    bool
    has_known_spatial_variation( const std::shared_ptr< Property >& aProperty )
    {
        return classify_property_spatial_variation( aProperty ) != Property_Spatial_Variation::UNKNOWN;
    }

    //------------------------------------------------------------------------------

    bool
    is_spatially_constant( const std::shared_ptr< Property >& aProperty )
    {
        return classify_property_spatial_variation( aProperty ) == Property_Spatial_Variation::CONSTANT;
    }

    //------------------------------------------------------------------------------

    bool
    has_known_spatial_variation_nonzero( const std::shared_ptr< Property >& aProperty )
    {
        const Property_Spatial_Variation tVariation = classify_property_spatial_variation( aProperty );

        if ( tVariation == Property_Spatial_Variation::DIFFERENTIABLE )
        {
            return aProperty->val()( 0 ) != 0.0;
        }

        return tVariation == Property_Spatial_Variation::CONSTANT;
    }

    //------------------------------------------------------------------------------

    void
    eval_property_space_derivative(
            const std::shared_ptr< Property >& aProperty,
            uint                               aDirection,
            Matrix< DDRMat >&                  adPropdxk )
    {
        // no property
        if ( aProperty == nullptr )
        {
            adPropdxk.set_size( 1, 1, 0.0 );
            return;
        }

        // number of entries of the property value
        const uint tNumEntries = aProperty->val().numel();

        adPropdxk.set_size( tNumEntries, 1, 0.0 );

        // constant property (the supports_analytical_* hooks exclude properties of unknown spatial variation)
        if ( !aProperty->check_space_dependency() )
        {
            return;
        }

        // space derivative, layout of Field_Interpolator::gradx: ( number of space dimensions x number of entries )
        const Matrix< DDRMat >& tdPropdx = aProperty->dPropdx( 1 );

        if ( tdPropdx.n_cols() == tNumEntries && aDirection < tdPropdx.n_rows() )
        {
            for ( uint iEntry = 0; iEntry < tNumEntries; iEntry++ )
            {
                adPropdxk( iEntry ) = tdPropdx( aDirection, iEntry );
            }
        }
        else if ( tNumEntries == 1 && aDirection < tdPropdx.numel() )
        {
            // scalar property given as a row vector
            adPropdxk( 0 ) = tdPropdx( aDirection );
        }
        else
        {
            MORIS_ERROR( false,
                    "eval_property_space_derivative - property %s: the space derivative ( %zu x %zu ) does not match "
                    "the layout ( number of space dimensions x %u ).",
                    aProperty->get_name().c_str(),
                    static_cast< size_t >( tdPropdx.n_rows() ),
                    static_cast< size_t >( tdPropdx.n_cols() ),
                    tNumEntries );
        }
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem
