/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_Geometric_Sensitivities_Struc_Linear.cpp
 *
 * Analytical geometric sensitivities (dR/dp and dQI/dp wrt the IG vertex coordinates) of the linear
 * structural IWGs and IQIs, checked against the finite difference geometry path
 * (IWG::compute_dRdp_FD_geometry, IQI::compute_dQIdp_FD_geometry) on a cut configuration:
 * an IG simplex (bulk) or an IG facet (side set) inside a rectangular IP element, in 2D and 3D,
 * for linear and quadratic Lagrange interpolation.
 */

#include <string>
#include <catch.hpp>
#include <memory>
#include "assert.hpp"

#define protected public
#define private public
// FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG.hpp"
#include "cl_FEM_IQI.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_CM_Struc_Linear.hpp"
#undef protected
#undef private
// MTK/src
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Integrator.hpp"
#include "cl_MTK_Integration_Rule.hpp"
// LINALG/src
#include "fn_norm.hpp"
#include "fn_max.hpp"
#include "fn_trans.hpp"
// FEM/INT/src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Geometry_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_IQI_Factory.hpp"

using namespace moris;
using namespace fem;

namespace
{
    //------------------------------------------------------------------------------

    // displacement dof type group (named: a single-element nested brace list is ambiguous)
    const Vector< Vector< MSI::Dof_Type > > tGeoSensDispDofTypes = { { MSI::Dof_Type::UX } };

    //------------------------------------------------------------------------------

    // relative FD perturbation of the IG vertex coordinates
    const real tGeoSensPerturbation = 1.0e-6;

    // required agreement between analytical and FD geometric sensitivities
    const real tGeoSensTolerance = 1.0e-6;

    //------------------------------------------------------------------------------

    void
    tConstValFunc_GeoSens(
            Matrix< DDRMat >&           aPropMatrix,
            Vector< Matrix< DDRMat > >& aParameters,
            Field_Interpolator_Manager* aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    //------------------------------------------------------------------------------

    // rectangular IP element [0,2] x [0,1.5] ( x [0,1] )
    real
    ip_element_length( uint aDirection )
    {
        const real tLengths[ 3 ] = { 2.0, 1.5, 1.0 };
        return tLengths[ aDirection ];
    }

    //------------------------------------------------------------------------------

    // IP element nodal coordinates from the reference coordinates of the interpolation rule
    Matrix< DDRMat >
    rectangular_ip_coordinates(
            Geometry_Interpolator& aIPGI,
            uint                   aSpaceDim )
    {
        aIPGI.set_param_coeff();
        Matrix< DDRMat > tXiHat = aIPGI.get_space_param_coeff();

        Matrix< DDRMat > tXHat( tXiHat.n_rows(), aSpaceDim, 0.0 );
        for ( uint iNode = 0; iNode < tXiHat.n_rows(); iNode++ )
        {
            for ( uint iDim = 0; iDim < aSpaceDim; iDim++ )
            {
                tXHat( iNode, iDim ) = 0.5 * ( tXiHat( iNode, iDim ) + 1.0 ) * ip_element_length( iDim );
            }
        }
        return tXHat;
    }

    //------------------------------------------------------------------------------

    // IP parametric coordinates of physical points inside the rectangular IP element
    Matrix< DDRMat >
    ip_param_coordinates( const Matrix< DDRMat >& aX )
    {
        Matrix< DDRMat > tXi( aX.n_rows(), aX.n_cols(), 0.0 );
        for ( uint iRow = 0; iRow < aX.n_rows(); iRow++ )
        {
            for ( uint iDim = 0; iDim < aX.n_cols(); iDim++ )
            {
                tXi( iRow, iDim ) = 2.0 * aX( iRow, iDim ) / ip_element_length( iDim ) - 1.0;
            }
        }
        return tXi;
    }

    //------------------------------------------------------------------------------

    // deterministic, non-polynomial displacement coefficients
    Matrix< DDRMat >
    displacement_coefficients(
            uint aNumBases,
            uint aNumFields )
    {
        Matrix< DDRMat > tUHat( aNumBases, aNumFields, 0.0 );
        for ( uint iBase = 0; iBase < aNumBases; iBase++ )
        {
            for ( uint iField = 0; iField < aNumFields; iField++ )
            {
                tUHat( iBase, iField ) = 0.1 * std::sin( 1.3 * iBase + 0.7 * iField + 0.2 ) + 0.03 * ( ( iBase + 2 * iField ) % 5 );
            }
        }
        return tUHat;
    }

    //------------------------------------------------------------------------------

    // IG cell vertex coordinates: a triangle/tetrahedron (bulk) or a segment/triangle (side) inside the IP element
    Matrix< DDRMat >
    ig_cell_coordinates(
            uint aSpaceDim,
            bool aIsSide )
    {
        if ( aSpaceDim == 2 )
        {
            if ( aIsSide )
            {
                return { { 0.3, 0.2 }, { 1.6, 0.9 } };
            }
            return { { 0.3, 0.2 }, { 1.6, 0.4 }, { 0.7, 1.2 } };
        }

        if ( aIsSide )
        {
            return { { 0.3, 0.2, 0.1 }, { 1.5, 0.3, 0.4 }, { 0.6, 1.2, 0.7 } };
        }
        return { { 0.3, 0.2, 0.1 }, { 1.5, 0.3, 0.2 }, { 0.6, 1.2, 0.3 }, { 0.7, 0.5, 0.8 } };
    }

    //------------------------------------------------------------------------------

    // max abs difference relative to max abs reference value
    real
    relative_difference(
            const Matrix< DDRMat >& aValue,
            const Matrix< DDRMat >& aReference )
    {
        real tMaxDiff = 0.0;
        real tMaxRef  = 0.0;
        for ( uint iEntry = 0; iEntry < aReference.numel(); iEntry++ )
        {
            tMaxDiff = std::max( tMaxDiff, std::abs( aValue( iEntry ) - aReference( iEntry ) ) );
            tMaxRef  = std::max( tMaxRef, std::abs( aReference( iEntry ) ) );
        }
        return tMaxDiff / std::max( tMaxRef, 1.0e-30 );
    }

    //------------------------------------------------------------------------------

    /**
     * The geometric setup shared by the IWG and IQI checks: IP/IG geometry interpolators, the displacement
     * field interpolator, the integration rule and the pdv assembly map of the IG cell.
     */
    class Geo_Sens_Setup
    {
      public:
        uint mSpaceDim;
        uint mOrder;
        bool mIsSide;

        Geometry_Interpolator* mIPGI = nullptr;
        Geometry_Interpolator* mIGGI = nullptr;

        Vector< Vector< MSI::Dof_Type > > mDofTypes;

        Matrix< DDRMat > mIntegPoints;
        Matrix< DDRMat > mIntegWeights;
        Matrix< DDSMat > mGeoLocalAssembly;

        uint mNumDofs = 0;
        uint mNumPdvs = 0;

        Geo_Sens_Setup(
                uint aSpaceDim,
                uint aOrder,
                bool aIsSide )
                : mSpaceDim( aSpaceDim )
                , mOrder( aOrder )
                , mIsSide( aIsSide )
        {
            // IP element: rectangular quad / hex
            mtk::Geometry_Type tIPGeometryType = aSpaceDim == 2 ? mtk::Geometry_Type::QUAD : mtk::Geometry_Type::HEX;

            mtk::Interpolation_Order tOrder = aOrder == 1 ? mtk::Interpolation_Order::LINEAR : mtk::Interpolation_Order::QUADRATIC;

            mtk::Interpolation_Rule tIPRule(
                    tIPGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            mIPGI = new Geometry_Interpolator( tIPRule, mtk::CellShape::RECTANGULAR, false, false );

            Matrix< DDRMat > tTHat = { { 0.0 }, { 1.0 } };
            Matrix< DDRMat > tXHat = rectangular_ip_coordinates( *mIPGI, aSpaceDim );
            mIPGI->set_coeff( tXHat, tTHat );

            // IG cell: linear simplex (bulk) or facet (side)
            mtk::Geometry_Type tIGGeometryType;
            if ( aSpaceDim == 2 )
            {
                tIGGeometryType = aIsSide ? mtk::Geometry_Type::LINE : mtk::Geometry_Type::TRI;
            }
            else
            {
                tIGGeometryType = aIsSide ? mtk::Geometry_Type::TRI : mtk::Geometry_Type::TET;
            }

            mtk::Interpolation_Rule tIGRule(
                    tIGGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            mIGGI = new Geometry_Interpolator( tIGRule, tIPRule, mtk::CellShape::GENERAL, aIsSide, false );

            Matrix< DDRMat > tXIG = ig_cell_coordinates( aSpaceDim, aIsSide );

            mIGGI->set_space_coeff( tXIG );
            mIGGI->set_time_coeff( tTHat );
            mIGGI->set_space_param_coeff( ip_param_coordinates( tXIG ) );
            mIGGI->set_time_param_coeff( { { -1.0 }, { 1.0 } } );

            // integration rule on the IG cell
            mtk::Integration_Order tIntegOrder;
            if ( aSpaceDim == 2 )
            {
                tIntegOrder = aIsSide ? mtk::Integration_Order::BAR_3 : mtk::Integration_Order::TRI_3;
            }
            else
            {
                tIntegOrder = aIsSide ? mtk::Integration_Order::TRI_3 : mtk::Integration_Order::TET_4;
            }

            mtk::Integration_Rule tIntegRule(
                    tIGGeometryType,
                    mtk::Integration_Type::GAUSS,
                    tIntegOrder,
                    mtk::Geometry_Type::LINE,
                    mtk::Integration_Type::GAUSS,
                    mtk::Integration_Order::BAR_1 );

            mtk::Integrator tIntegrator( tIntegRule );
            tIntegrator.get_points( mIntegPoints );
            tIntegrator.get_weights( mIntegWeights );

            // all IG vertex coordinates are active pdvs
            uint tNumVertices = tXIG.n_rows();
            mGeoLocalAssembly.set_size( tNumVertices, aSpaceDim, -1 );
            for ( uint iVertex = 0; iVertex < tNumVertices; iVertex++ )
            {
                for ( uint iDim = 0; iDim < aSpaceDim; iDim++ )
                {
                    mGeoLocalAssembly( iVertex, iDim ) = iVertex * aSpaceDim + iDim;
                }
            }
            mNumPdvs = tNumVertices * aSpaceDim;

            // displacement dofs
            if ( aSpaceDim == 2 )
            {
                mDofTypes = { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } };
            }
            else
            {
                mDofTypes = { { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } };
            }
        }

        ~Geo_Sens_Setup()
        {
            delete mIPGI;
            delete mIGGI;
        }

        // create the displacement field interpolator on the IP element
        Field_Interpolator*
        create_displacement_FI()
        {
            mtk::Geometry_Type tIPGeometryType = mSpaceDim == 2 ? mtk::Geometry_Type::QUAD : mtk::Geometry_Type::HEX;

            mtk::Interpolation_Order tOrder = mOrder == 1 ? mtk::Interpolation_Order::LINEAR : mtk::Interpolation_Order::QUADRATIC;

            mtk::Interpolation_Rule tFIRule(
                    tIPGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            Field_Interpolator* tFI = new Field_Interpolator( mSpaceDim, tFIRule, mIPGI, mDofTypes( 0 ) );

            uint tNumBases = tFI->get_number_of_space_time_bases();
            tFI->set_coeff( displacement_coefficients( tNumBases, mSpaceDim ) );

            mNumDofs = tNumBases * mSpaceDim;

            return tFI;
        }

        // set the set maps for the displacement dof type
        void
        populate_set_dof_maps( Set* aSet )
        {
            aSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

            aSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            aSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

            aSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            aSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

            aSet->mResDofAssemblyMap.resize( 1 );
            aSet->mResDofAssemblyMap( 0 ) = { { 0, static_cast< sint >( mNumDofs ) - 1 } };

            aSet->mJacDofAssemblyMap.resize( 1 );
            aSet->mJacDofAssemblyMap( 0 ) = { { 0, static_cast< sint >( mNumDofs ) - 1 } };

            aSet->mResidual.resize( 1 );
            aSet->mResidual( 0 ).set_size( mNumDofs, 1, 0.0 );
            aSet->mJacobian.set_size( mNumDofs, mNumDofs, 0.0 );
        }

        // set the evaluation point and return the integration weight
        real
        set_evaluation_point(
                Field_Interpolator_Manager& aFIManager,
                uint                        aGP )
        {
            Matrix< DDRMat > tLocalIntegPoint = mIntegPoints.get_column( aGP );
            aFIManager.set_space_time_from_local_IG_point( tLocalIntegPoint );

            return mIntegWeights( aGP ) * mIGGI->det_J();
        }
    };

    //------------------------------------------------------------------------------

    // elastic CM, plane stress in 2D
    std::shared_ptr< Constitutive_Model >
    create_elastic_CM(
            uint       aSpaceDim,
            Model_Type aModelType = Model_Type::PLANE_STRESS )
    {
        std::shared_ptr< Property > tPropEMod = std::make_shared< Property >();
        tPropEMod->set_parameters( { { { 2.0 } } } );
        tPropEMod->set_val_function( tConstValFunc_GeoSens );

        std::shared_ptr< Property > tPropNu = std::make_shared< Property >();
        tPropNu->set_parameters( { { { 0.3 } } } );
        tPropNu->set_val_function( tConstValFunc_GeoSens );

        CM_Factory tCMFactory;

        std::shared_ptr< Constitutive_Model > tCM = tCMFactory.create_CM( Constitutive_Type::STRUC_LIN_ISO );

        tCM->set_dof_type_list( tGeoSensDispDofTypes );
        tCM->set_property( tPropEMod, "YoungsModulus" );
        tCM->set_property( tPropNu, "PoissonRatio" );
        tCM->set_local_properties();
        tCM->set_model_type( aModelType );
        tCM->set_space_dim( aSpaceDim );

        return tCM;
    }

    //------------------------------------------------------------------------------

    // constant property
    std::shared_ptr< Property >
    create_const_property( const Matrix< DDRMat >& aValue )
    {
        std::shared_ptr< Property > tProp = std::make_shared< Property >();
        tProp->set_parameters( { aValue } );
        tProp->set_val_function( tConstValFunc_GeoSens );
        return tProp;
    }

    //------------------------------------------------------------------------------

    /**
     * compare analytical and FD dR/dp_geo of an IWG at all integration points of the IG cell
     * @return max relative difference over the integration points
     */
    real
    check_IWG_geometric_sensitivities(
            const std::shared_ptr< IWG >& aIWG,
            uint                          aSpaceDim,
            uint                          aOrder,
            bool                          aIsSide )
    {
        Geo_Sens_Setup tSetup( aSpaceDim, aOrder, aIsSide );

        // set
        Set* tSet = new Set();
        tSet->set_set_type( aIsSide ? Element_Type::SIDESET : Element_Type::BULK );
        aIWG->set_set_pointer( tSet );

        // field interpolators and set maps
        Vector< Field_Interpolator* > tFIs( 1, tSetup.create_displacement_FI() );
        tSetup.populate_set_dof_maps( tSet );

        aIWG->get_global_dof_type_list();
        aIWG->mRequestedLeaderGlobalDofTypes = tGeoSensDispDofTypes;

        // field interpolator manager
        Vector< Vector< gen::PDV_Type > >   tDummyDv;
        Vector< Vector< mtk::Field_Type > > tDummyField;
        Field_Interpolator_Manager          tFIManager( tGeoSensDispDofTypes, tDummyDv, tDummyField, tSet );

        tFIManager.mFI                     = tFIs;
        tFIManager.mIPGeometryInterpolator = tSetup.mIPGI;
        tFIManager.mIGGeometryInterpolator = tSetup.mIGGI;

        tSet->mLeaderFIManager = &tFIManager;
        aIWG->set_field_interpolator_manager( &tFIManager );

        // geometry sensitivities
        tSet->mdRdp.resize( 2 );
        tSet->mdRdp( 1 ).set_size( tSetup.mNumDofs, tSetup.mNumPdvs, 0.0 );

        // the cluster measure contribution is shared by the FD and the analytical path
        aIWG->mActiveCMEAFlag = false;

        // the IWG must report analytical geometric sensitivities
        REQUIRE( aIWG->has_analytical_dRdp_geometry() );

        Vector< Matrix< IndexMat > > tVertexIndices( 0 );

        real tMaxError = 0.0;
        for ( uint iGP = 0; iGP < tSetup.mIntegPoints.n_cols(); iGP++ )
        {
            real tWStar = tSetup.set_evaluation_point( tFIManager, iGP );

            if ( aIsSide )
            {
                Matrix< DDRMat > tNormal;
                tSetup.mIGGI->get_normal( tNormal );
                aIWG->set_normal( tNormal );
            }

            // FD reference
            aIWG->reset_eval_flags();
            tSet->mdRdp( 1 ).fill( 0.0 );
            aIWG->compute_dRdp_FD_geometry(
                    tWStar,
                    tGeoSensPerturbation,
                    FDScheme_Type::POINT_3_CENTRAL,
                    tSetup.mGeoLocalAssembly,
                    tVertexIndices );
            Matrix< DDRMat > tdRdpFD = tSet->mdRdp( 1 );

            // analytical
            aIWG->reset_eval_flags();
            tSet->mdRdp( 1 ).fill( 0.0 );
            aIWG->compute_dRdp_geometry_analytical(
                    tWStar,
                    tGeoSensPerturbation,
                    FDScheme_Type::POINT_3_CENTRAL,
                    tSetup.mGeoLocalAssembly );
            Matrix< DDRMat > tdRdp = tSet->mdRdp( 1 );

            REQUIRE( norm( tdRdpFD ) > 0.0 );

            tMaxError = std::max( tMaxError, relative_difference( tdRdp, tdRdpFD ) );
        }

        // the FI manager deletes the field interpolators
        tFIManager.mIPGeometryInterpolator = nullptr;
        tFIManager.mIGGeometryInterpolator = nullptr;
        tSet->mLeaderFIManager             = nullptr;
        delete tSet;

        return tMaxError;
    }

    //------------------------------------------------------------------------------

    /**
     * compare analytical and FD dQI/dp_geo of an IQI at all integration points of the IG cell
     * @return max relative difference over the integration points
     */
    real
    check_IQI_geometric_sensitivities(
            const std::shared_ptr< IQI >& aIQI,
            uint                          aSpaceDim,
            uint                          aOrder )
    {
        Geo_Sens_Setup tSetup( aSpaceDim, aOrder, false );

        // set
        Set* tSet = new Set();
        tSet->set_set_type( Element_Type::BULK );
        aIQI->set_name( "QI" );
        aIQI->set_set_pointer( tSet );

        // field interpolators and set maps
        Vector< Field_Interpolator* > tFIs( 1, tSetup.create_displacement_FI() );
        tSetup.populate_set_dof_maps( tSet );

        aIQI->mRequestedLeaderGlobalDofTypes = tGeoSensDispDofTypes;

        // QI storage
        tSet->mRequestedIQINamesAssemblyMap[ "QI" ] = 0;
        tSet->mQI.resize( 1 );
        tSet->mQI( 0 ).set_size( 1, 1, 0.0 );

        // field interpolator manager
        Vector< Vector< gen::PDV_Type > >   tDummyDv;
        Vector< Vector< mtk::Field_Type > > tDummyField;
        Field_Interpolator_Manager          tFIManager( tGeoSensDispDofTypes, tDummyDv, tDummyField, tSet );

        tFIManager.mFI                     = tFIs;
        tFIManager.mIPGeometryInterpolator = tSetup.mIPGI;
        tFIManager.mIGGeometryInterpolator = tSetup.mIGGI;

        tSet->mLeaderFIManager = &tFIManager;
        aIQI->set_field_interpolator_manager( &tFIManager );

        // geometry sensitivities
        tSet->mdQIdp.resize( 2 );
        tSet->mdQIdp( 1 ).resize( 1 );
        tSet->mdQIdp( 1 )( 0 ).set_size( 1, tSetup.mNumPdvs, 0.0 );

        // the cluster measure contribution is shared by the FD and the analytical path
        aIQI->mActiveCMEAFlag = false;

        // the IQI must report analytical geometric sensitivities
        REQUIRE( aIQI->has_analytical_dQIdp_geometry() );

        Vector< Matrix< IndexMat > > tVertexIndices( 0 );

        real tMaxError = 0.0;
        for ( uint iGP = 0; iGP < tSetup.mIntegPoints.n_cols(); iGP++ )
        {
            real tWStar = tSetup.set_evaluation_point( tFIManager, iGP );

            // FD reference
            aIQI->reset_eval_flags();
            tSet->mdQIdp( 1 )( 0 ).fill( 0.0 );
            aIQI->compute_dQIdp_FD_geometry(
                    tWStar,
                    tGeoSensPerturbation,
                    FDScheme_Type::POINT_3_CENTRAL,
                    tSetup.mGeoLocalAssembly,
                    tVertexIndices );
            Matrix< DDRMat > tdQIdpFD = tSet->mdQIdp( 1 )( 0 );

            // analytical
            aIQI->reset_eval_flags();
            tSet->mdQIdp( 1 )( 0 ).fill( 0.0 );
            aIQI->compute_dQIdp_geometry_analytical(
                    tWStar,
                    tGeoSensPerturbation,
                    FDScheme_Type::POINT_3_CENTRAL,
                    tSetup.mGeoLocalAssembly );
            Matrix< DDRMat > tdQIdp = tSet->mdQIdp( 1 )( 0 );

            REQUIRE( norm( tdQIdpFD ) > 0.0 );

            tMaxError = std::max( tMaxError, relative_difference( tdQIdp, tdQIdpFD ) );
        }

        tFIManager.mIPGeometryInterpolator = nullptr;
        tFIManager.mIGGeometryInterpolator = nullptr;
        tSet->mLeaderFIManager             = nullptr;
        delete tSet;

        return tMaxError;
    }

    //------------------------------------------------------------------------------
}    // namespace

//------------------------------------------------------------------------------

TEST_CASE( "Analytic_Geo_Sens_Test_Strain_Layout", "[fem],[Analytic_Geo_Sens],[Analytic_Geo_Sens_Test_Strain_Layout]" )
{
    // the strain-displacement matrix built from dN/dx must reproduce CM_Struc_Linear::testStrain()
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        Vector< Model_Type > tModelTypes = { Model_Type::PLANE_STRESS };
        if ( iSpaceDim == 2 )
        {
            tModelTypes.push_back( Model_Type::PLANE_STRAIN );
        }

        for ( Model_Type tModelType : tModelTypes )
        {
            Geo_Sens_Setup tSetup( iSpaceDim, 2, false );

            Set tSet;
            tSet.set_set_type( Element_Type::BULK );

            Vector< Field_Interpolator* > tFIs( 1, tSetup.create_displacement_FI() );
            tSetup.populate_set_dof_maps( &tSet );

            Vector< Vector< gen::PDV_Type > >   tDummyDv;
            Vector< Vector< mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager          tFIManager( tGeoSensDispDofTypes, tDummyDv, tDummyField, &tSet );
            tFIManager.mFI                     = tFIs;
            tFIManager.mIPGeometryInterpolator = tSetup.mIPGI;
            tFIManager.mIGGeometryInterpolator = tSetup.mIGGI;

            std::shared_ptr< Constitutive_Model > tCM = create_elastic_CM( iSpaceDim, tModelType );
            tCM->set_field_interpolator_manager( &tFIManager );

            tSetup.set_evaluation_point( tFIManager, 0 );
            tCM->reset_eval_flags();

            CM_Struc_Linear* tCMStrucLinear = dynamic_cast< CM_Struc_Linear* >( tCM.get() );
            REQUIRE( tCMStrucLinear != nullptr );
            REQUIRE( tCMStrucLinear->has_spatial_strain_derivatives() );

            Matrix< DDRMat > tB;
            tCMStrucLinear->build_test_strain_from_gradient( tFIs( 0 )->dnNdxn( 1 ), tB );

            CHECK( relative_difference( tB, tCM->testStrain() ) < 1.0e-14 );

            tFIManager.mIPGeometryInterpolator = nullptr;
            tFIManager.mIGGeometryInterpolator = nullptr;
        }
    }
}

//------------------------------------------------------------------------------

TEST_CASE( "Analytic_Geo_Sens_IWG_Struc_Linear_Bulk", "[fem],[Analytic_Geo_Sens],[Analytic_Geo_Sens_IWG_Bulk]" )
{
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        for ( uint iOrder = 1; iOrder < 3; iOrder++ )
        {
            IWG_Factory tIWGFactory;

            std::shared_ptr< IWG > tIWG = tIWGFactory.create_IWG( IWG_Type::STRUC_LINEAR_BULK );
            tIWG->set_residual_dof_type( tGeoSensDispDofTypes );
            tIWG->set_dof_type_list( tGeoSensDispDofTypes, mtk::Leader_Follower::LEADER );
            tIWG->set_constitutive_model( create_elastic_CM( iSpaceDim ), "ElastLinIso" );

            // bedding and body load
            tIWG->set_property( create_const_property( { { 0.7 } } ), "Bedding" );
            tIWG->set_property( create_const_property( Matrix< DDRMat >( iSpaceDim, 1, 0.4 ) ), "Load" );

            real tError = check_IWG_geometric_sensitivities( tIWG, iSpaceDim, iOrder, false );

            std::cout << "Analytic_Geo_Sens IWG_Struc_Linear_Bulk dim " << iSpaceDim << " order " << iOrder
                      << " max rel. error vs FD " << tError << '\n';

            CHECK( tError < tGeoSensTolerance );
        }
    }
}

//------------------------------------------------------------------------------

TEST_CASE( "Analytic_Geo_Sens_IWG_Struc_Linear_Bulk_Plane_Strain", "[fem],[Analytic_Geo_Sens],[Analytic_Geo_Sens_IWG_Bulk]" )
{
    IWG_Factory tIWGFactory;

    std::shared_ptr< IWG > tIWG = tIWGFactory.create_IWG( IWG_Type::STRUC_LINEAR_BULK );
    tIWG->set_residual_dof_type( tGeoSensDispDofTypes );
    tIWG->set_dof_type_list( tGeoSensDispDofTypes, mtk::Leader_Follower::LEADER );
    tIWG->set_constitutive_model( create_elastic_CM( 2, Model_Type::PLANE_STRAIN ), "ElastLinIso" );

    real tError = check_IWG_geometric_sensitivities( tIWG, 2, 2, false );

    std::cout << "Analytic_Geo_Sens IWG_Struc_Linear_Bulk plane strain max rel. error vs FD " << tError << '\n';

    CHECK( tError < tGeoSensTolerance );
}

//------------------------------------------------------------------------------

TEST_CASE( "Analytic_Geo_Sens_IWG_Struc_Linear_Neumann", "[fem],[Analytic_Geo_Sens],[Analytic_Geo_Sens_IWG_Neumann]" )
{
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        for ( uint iOrder = 1; iOrder < 3; iOrder++ )
        {
            IWG_Factory tIWGFactory;

            std::shared_ptr< IWG > tIWG = tIWGFactory.create_IWG( IWG_Type::STRUC_LINEAR_NEUMANN );
            tIWG->set_residual_dof_type( tGeoSensDispDofTypes );
            tIWG->set_dof_type_list( tGeoSensDispDofTypes, mtk::Leader_Follower::LEADER );

            // traction and pressure
            Matrix< DDRMat > tTraction( iSpaceDim, 1, 0.3 );
            tTraction( 0 ) = -0.6;
            tIWG->set_property( create_const_property( tTraction ), "Traction" );
            tIWG->set_property( create_const_property( { { 0.8 } } ), "Pressure" );

            real tError = check_IWG_geometric_sensitivities( tIWG, iSpaceDim, iOrder, true );

            std::cout << "Analytic_Geo_Sens IWG_Struc_Linear_Neumann dim " << iSpaceDim << " order " << iOrder
                      << " max rel. error vs FD " << tError << '\n';

            CHECK( tError < tGeoSensTolerance );
        }
    }
}

//------------------------------------------------------------------------------

TEST_CASE( "Analytic_Geo_Sens_IWG_Struc_Linear_Dirichlet", "[fem],[Analytic_Geo_Sens],[Analytic_Geo_Sens_IWG_Dirichlet]" )
{
    Vector< IWG_Type > tIWGTypes = {
        IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE,
        IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE
    };

    for ( IWG_Type tIWGType : tIWGTypes )
    {
        for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
        {
            for ( uint iOrder = 1; iOrder < 3; iOrder++ )
            {
                std::shared_ptr< Constitutive_Model > tCM = create_elastic_CM( iSpaceDim );

                // Nitsche parameter, constant through its (empty) cluster
                std::shared_ptr< Property > tPropEMod = create_const_property( { { 2.0 } } );

                SP_Factory tSPFactory;

                std::shared_ptr< Stabilization_Parameter > tSPNitsche =
                        tSPFactory.create_SP( Stabilization_Type::DIRICHLET_NITSCHE );
                tSPNitsche->set_parameters( { { { 100.0 } } } );
                tSPNitsche->set_property( tPropEMod, "Material", mtk::Leader_Follower::LEADER );

                Cluster tCluster;
                tSPNitsche->set_cluster( &tCluster );

                IWG_Factory tIWGFactory;

                std::shared_ptr< IWG > tIWG = tIWGFactory.create_IWG( tIWGType );
                tIWG->set_residual_dof_type( tGeoSensDispDofTypes );
                tIWG->set_dof_type_list( tGeoSensDispDofTypes, mtk::Leader_Follower::LEADER );
                tIWG->set_stabilization_parameter( tSPNitsche, "DirichletNitsche" );
                tIWG->set_constitutive_model( tCM, "ElastLinIso" );

                Matrix< DDRMat > tImposedDisp( iSpaceDim, 1, 0.05 );
                tIWG->set_property( create_const_property( tImposedDisp ), "Dirichlet" );

                real tError = check_IWG_geometric_sensitivities( tIWG, iSpaceDim, iOrder, true );

                std::cout << "Analytic_Geo_Sens IWG_Struc_Linear_Dirichlet type " << static_cast< uint >( tIWGType )
                          << " dim " << iSpaceDim << " order " << iOrder << " max rel. error vs FD " << tError << '\n';

                CHECK( tError < tGeoSensTolerance );
            }
        }
    }
}

//------------------------------------------------------------------------------

TEST_CASE( "Analytic_Geo_Sens_IQI_Volume", "[fem],[Analytic_Geo_Sens],[Analytic_Geo_Sens_IQI_Volume]" )
{
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        IQI_Factory tIQIFactory;

        std::shared_ptr< IQI > tIQI = tIQIFactory.create_IQI( IQI_Type::VOLUME );
        tIQI->set_property( create_const_property( { { 1.3 } } ), "Density", mtk::Leader_Follower::LEADER );

        real tError = check_IQI_geometric_sensitivities( tIQI, iSpaceDim, 1 );

        std::cout << "Analytic_Geo_Sens IQI_Volume dim " << iSpaceDim << " max rel. error vs FD " << tError << '\n';

        CHECK( tError < tGeoSensTolerance );
    }
}

//------------------------------------------------------------------------------

TEST_CASE( "Analytic_Geo_Sens_IQI_Strain_Energy", "[fem],[Analytic_Geo_Sens],[Analytic_Geo_Sens_IQI_Strain_Energy]" )
{
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        for ( uint iOrder = 1; iOrder < 3; iOrder++ )
        {
            IQI_Factory tIQIFactory;

            std::shared_ptr< IQI > tIQI = tIQIFactory.create_IQI( IQI_Type::STRAIN_ENERGY );
            tIQI->set_dof_type_list( tGeoSensDispDofTypes, mtk::Leader_Follower::LEADER );
            tIQI->set_constitutive_model( create_elastic_CM( iSpaceDim ), "Elast", mtk::Leader_Follower::LEADER );
            tIQI->set_property( create_const_property( { { 0.7 } } ), "Bedding", mtk::Leader_Follower::LEADER );

            real tError = check_IQI_geometric_sensitivities( tIQI, iSpaceDim, iOrder );

            std::cout << "Analytic_Geo_Sens IQI_Strain_Energy dim " << iSpaceDim << " order " << iOrder
                      << " max rel. error vs FD " << tError << '\n';

            CHECK( tError < tGeoSensTolerance );
        }
    }
}

//------------------------------------------------------------------------------

TEST_CASE( "Analytic_Geo_Sens_IQI_Max_Stress", "[fem],[Analytic_Geo_Sens],[Analytic_Geo_Sens_IQI_Max_Stress]" )
{
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        for ( uint iOrder = 1; iOrder < 3; iOrder++ )
        {
            IQI_Factory tIQIFactory;

            std::shared_ptr< IQI > tIQI = tIQIFactory.create_IQI( IQI_Type::MAX_VON_MISES_STRESS );
            tIQI->set_dof_type_list( tGeoSensDispDofTypes, mtk::Leader_Follower::LEADER );
            tIQI->set_constitutive_model( create_elastic_CM( iSpaceDim ), "ElastLinIso", mtk::Leader_Follower::LEADER );
            tIQI->set_property( create_const_property( { { 0.5 } } ), "ReferenceValue", mtk::Leader_Follower::LEADER );
            tIQI->set_property( create_const_property( { { 8.0 } } ), "Exponent", mtk::Leader_Follower::LEADER );
            tIQI->set_property( create_const_property( { { 0.0 } } ), "Shift", mtk::Leader_Follower::LEADER );

            real tError = check_IQI_geometric_sensitivities( tIQI, iSpaceDim, iOrder );

            std::cout << "Analytic_Geo_Sens IQI_Max_Stress dim " << iSpaceDim << " order " << iOrder
                      << " max rel. error vs FD " << tError << '\n';

            CHECK( tError < tGeoSensTolerance );
        }
    }
}

//------------------------------------------------------------------------------
