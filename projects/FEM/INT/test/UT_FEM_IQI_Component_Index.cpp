/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IQI_Component_Index.cpp
 *
 * Component selection by vectorial_field_index: a set index is used as given and
 * never reset by a default (MB-054, MB-055, MB-060), and IQI_Dof resolves scalar
 * fields and single members of a vector group (IQI audit B-1, B-2, B-3).
 * Exact values on one QUAD4 with u = b + A x, T = c + g.x.
 *
 */

#include <string>
#include <catch.hpp>
#include <memory>

#include "assert.hpp"

#define protected public
#define private public
// FEM/INT/src
#include "cl_FEM_Model.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Stabilization_Parameter.hpp"
#undef protected
#undef private
// FEM/INT/src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_IQI_Factory.hpp"

using namespace moris;
using namespace fem;

namespace
{
    // a stabilization parameter with two constant values, like (tauM, tauC)
    class SP_Two_Values_Test : public Stabilization_Parameter
    {
      public:
        void
        eval_SP() override
        {
            mPPVal = { { 2.0 }, { 5.0 } };
        }
    };

    // one QUAD4 (unit square), linear in time, with u = b + A x and T = c + g.x
    struct Component_Index_Fixture
    {
        Vector< Vector< MSI::Dof_Type > > mDofTypes = { { MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::TEMP } };

        Geometry_Interpolator*        mGI = nullptr;
        Vector< Field_Interpolator* > mFIs;
        fem::Set*                     mSet       = nullptr;
        Field_Interpolator_Manager*   mFIManager = nullptr;

        Component_Index_Fixture()
        {
            mtk::Interpolation_Rule tRule( mtk::Geometry_Type::QUAD,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            mGI = new Geometry_Interpolator( tRule );

            Matrix< DDRMat > tXHat = { { 0.0, 0.0 }, { 1.0, 0.0 }, { 1.0, 1.0 }, { 0.0, 1.0 } };
            Matrix< DDRMat > tTHat = { { 0.0 }, { 1.0 } };
            mGI->set_coeff( tXHat, tTHat );

            // nodal values, repeated for both time levels
            Matrix< DDRMat > tUHat( 8, 2, 0.0 );
            Matrix< DDRMat > tTempHat( 8, 1, 0.0 );
            for ( uint iNode = 0; iNode < 8; iNode++ )
            {
                real tX = tXHat( iNode % 4, 0 );
                real tY = tXHat( iNode % 4, 1 );

                tUHat( iNode, 0 )    = 0.1 + 1.0 * tX + 2.0 * tY;
                tUHat( iNode, 1 )    = 0.2 + 3.0 * tX + 4.0 * tY;
                tTempHat( iNode, 0 ) = 5.0 + 6.0 * tX + 7.0 * tY;
            }

            mFIs.resize( 2 );
            mFIs( 0 ) = new Field_Interpolator( 2, tRule, mGI, mDofTypes( 0 ) );
            mFIs( 0 )->set_coeff( tUHat );
            mFIs( 1 ) = new Field_Interpolator( 1, tRule, mGI, mDofTypes( 1 ) );
            mFIs( 1 )->set_coeff( tTempHat );

            // set with VX, VY in FI 0 and TEMP in FI 1
            mSet = new fem::Set();
            mSet->set_set_type( fem::Element_Type::BULK );
            mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );
            mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            for ( MSI::Dof_Type tType : { MSI::Dof_Type::VX, MSI::Dof_Type::VY } )
            {
                mSet->mUniqueDofTypeMap( static_cast< int >( tType ) ) = 0;
                mSet->mLeaderDofTypeMap( static_cast< int >( tType ) ) = 0;
            }
            mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;
            mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;

            // residual: 16 velocity dofs, then 8 temperature dofs
            mSet->mResDofAssemblyMap.resize( 2 );
            mSet->mResDofAssemblyMap( 0 ) = { { 0, 15 } };
            mSet->mResDofAssemblyMap( 1 ) = { { 16, 23 } };
            mSet->mResidual.resize( 1 );
            mSet->mResidual( 0 ).set_size( 24, 1, 0.0 );

            moris::map< std::string, moris_index > tRequestedIQINamesAssemblyMap;
            tRequestedIQINamesAssemblyMap[ "IQI" ] = 0;
            mSet->mRequestedIQINamesAssemblyMap    = tRequestedIQINamesAssemblyMap;
            mSet->mQI.resize( 1 );
            mSet->mQI( 0 ).set_size( 1, 1, 0.0 );

            Vector< Vector< gen::PDV_Type > >   tDummyDv;
            Vector< Vector< mtk::Field_Type > > tDummyField;
            mFIManager                          = new Field_Interpolator_Manager( mDofTypes, tDummyDv, tDummyField, mSet );
            mFIManager->mFI                     = mFIs;
            mFIManager->mIPGeometryInterpolator = mGI;
            mFIManager->mIGGeometryInterpolator = mGI;
            mSet->mLeaderFIManager              = mFIManager;

            // evaluate at the element centre, x = ( 0.5, 0.5 )
            Matrix< DDRMat > tParamPoint = { { 0.0 }, { 0.0 }, { 0.0 } };
            mFIManager->set_space_time( tParamPoint );
        }

        ~Component_Index_Fixture()
        {
            // the set owns the FI manager, which owns the FIs
            delete mSet;
            delete mGI;
        }

        std::shared_ptr< IQI >
        make_dof_iqi(
                const Vector< MSI::Dof_Type >& aQuantity,
                sint                           aIndex )
        {
            fem::IQI_Factory       tFactory;
            std::shared_ptr< IQI > tIQI = tFactory.create_IQI( fem::IQI_Type::DOF );
            tIQI->set_name( "IQI" );
            tIQI->set_dof_type_list( mDofTypes, mtk::Leader_Follower::LEADER );
            tIQI->set_quantity_dof_type( aQuantity );
            tIQI->set_output_type_index( aIndex );
            tIQI->set_set_pointer( mSet );
            tIQI->get_global_dof_type_list();
            tIQI->mRequestedLeaderGlobalDofTypes = mDofTypes;
            tIQI->set_field_interpolator_manager( mFIManager );
            return tIQI;
        }

        real
        global_qi( const std::shared_ptr< IQI >& aIQI )
        {
            mSet->mQI( 0 ).fill( 0.0 );
            aIQI->compute_QI( 1.0 );
            return mSet->mQI( 0 )( 0 );
        }
    };

    // sum of the dQIdu entries in rows [ aStart, aStop ]
    real
    block_sum(
            const Matrix< DDRMat >& adQIdu,
            uint                    aStart,
            uint                    aStop )
    {
        real tSum = 0.0;
        for ( uint iRow = aStart; iRow <= aStop; iRow++ )
        {
            tSum += adQIdu( iRow );
        }
        return tSum;
    }
}    // namespace

//------------------------------------------------------------------------------

TEST_CASE( "IQI_Component_Index_Helper", "[moris],[fem],[IQI_Component_Index]" )
{
    fem::IQI_Factory       tFactory;
    std::shared_ptr< IQI > tIQI = tFactory.create_IQI( fem::IQI_Type::DOF );

    // unset: default, or an error if there is no default
    tIQI->set_output_type_index( -1 );
    CHECK( tIQI->get_component_index( 3 ) == 0 );
    CHECK( tIQI->get_component_index( 3, 2 ) == 2 );
    CHECK_THROWS( tIQI->get_component_index( 3, -1 ) );

    // set: used as given, never replaced by the default
    tIQI->set_output_type_index( 1 );
    CHECK( tIQI->get_component_index( 3 ) == 1 );
    CHECK( tIQI->get_component_index( 3, 0 ) == 1 );

    // out of range
    tIQI->set_output_type_index( 3 );
    CHECK_THROWS( tIQI->get_component_index( 3 ) );
}

//------------------------------------------------------------------------------

TEST_CASE( "IQI_Dof_Component_Index", "[moris],[fem],[IQI_Component_Index]" )
{
    Component_Index_Fixture tFix;

    real tUX   = 0.1 + 0.5 * 1.0 + 0.5 * 2.0;    // 1.6
    real tUY   = 0.2 + 0.5 * 3.0 + 0.5 * 4.0;    // 3.7
    real tTemp = 5.0 + 0.5 * 6.0 + 0.5 * 7.0;    // 11.5

    SECTION( "vector group with an index (regression)" )
    {
        std::shared_ptr< IQI > tIQI = tFix.make_dof_iqi( { MSI::Dof_Type::VX, MSI::Dof_Type::VY }, 1 );
        CHECK( tFix.global_qi( tIQI ) == Approx( tUY ) );

        Matrix< DDRMat > tQI;
        tIQI->compute_QI( tQI );
        CHECK( tQI( 0 ) == Approx( tUY ) );
    }

    SECTION( "vector group without an index: whole vector for output, error as a scalar QI" )
    {
        std::shared_ptr< IQI > tIQI = tFix.make_dof_iqi( { MSI::Dof_Type::VX, MSI::Dof_Type::VY }, -1 );

        Matrix< DDRMat > tQI;
        tIQI->compute_QI( tQI );
        REQUIRE( tQI.numel() == 2 );
        CHECK( tQI( 0 ) == Approx( tUX ) );
        CHECK( tQI( 1 ) == Approx( tUY ) );

        CHECK_THROWS( tFix.global_qi( tIQI ) );
    }

    SECTION( "B-2: a single member of a vector group is that component" )
    {
        std::shared_ptr< IQI > tIQI = tFix.make_dof_iqi( { MSI::Dof_Type::VY }, -1 );
        CHECK( tFix.global_qi( tIQI ) == Approx( tUY ) );

        Matrix< DDRMat > tQI;
        tIQI->compute_QI( tQI );
        REQUIRE( tQI.numel() == 1 );
        CHECK( tQI( 0 ) == Approx( tUY ) );

        // dQIdu: the VY shape functions (sum 1), nothing on VX
        Matrix< DDRMat >        tdQIdu;
        Vector< MSI::Dof_Type > tVel = { MSI::Dof_Type::VX, MSI::Dof_Type::VY };
        tIQI->compute_dQIdu( tVel, tdQIdu );
        REQUIRE( tdQIdu.numel() == 16 );
        CHECK( block_sum( tdQIdu, 0, 7 ) == Approx( 0.0 ).margin( 1e-12 ) );
        CHECK( block_sum( tdQIdu, 8, 15 ) == Approx( 1.0 ) );

        // same through the set residual
        tFix.mSet->mResidual( 0 ).fill( 0.0 );
        tIQI->compute_dQIdu( 1.0 );
        CHECK( block_sum( tFix.mSet->mResidual( 0 ), 0, 7 ) == Approx( 0.0 ).margin( 1e-12 ) );
        CHECK( block_sum( tFix.mSet->mResidual( 0 ), 8, 15 ) == Approx( 1.0 ) );
    }

    SECTION( "a set index is not reset for a single dof_quantity" )
    {
        std::shared_ptr< IQI > tIQI = tFix.make_dof_iqi( { MSI::Dof_Type::VX }, 1 );
        CHECK( tFix.global_qi( tIQI ) == Approx( tUY ) );
        CHECK( tFix.global_qi( tIQI ) == Approx( tUY ) );
    }

    SECTION( "B-1: scalar field, analytic dQIdu" )
    {
        std::shared_ptr< IQI > tIQI = tFix.make_dof_iqi( { MSI::Dof_Type::TEMP }, -1 );
        CHECK( tFix.global_qi( tIQI ) == Approx( tTemp ) );

        tFix.mSet->mResidual( 0 ).fill( 0.0 );
        tIQI->compute_dQIdu( 1.0 );
        CHECK( block_sum( tFix.mSet->mResidual( 0 ), 0, 15 ) == Approx( 0.0 ).margin( 1e-12 ) );
        CHECK( block_sum( tFix.mSet->mResidual( 0 ), 16, 23 ) == Approx( 1.0 ) );
    }

    SECTION( "B-3: scalar field, spatial derivative without an index (output path)" )
    {
        std::shared_ptr< IQI > tIQI = tFix.make_dof_iqi( { MSI::Dof_Type::TEMP }, -1 );
        tIQI->set_parameters( { { { 1, 1 } } } );    // d/dy, first order

        Matrix< DDRMat > tQI;
        tIQI->compute_QI( tQI );
        REQUIRE( tQI.numel() == 1 );
        CHECK( tQI( 0 ) == Approx( 7.0 ) );
    }
}

//------------------------------------------------------------------------------

TEST_CASE( "IQI_Eigen_Value_Index", "[moris],[fem],[IQI_Component_Index],[MB-054]" )
{
    Component_Index_Fixture tFix;

    fem::FEM_Model tModel;
    tModel.set_eigen_values( std::make_shared< Vector< real > >( Vector< real >{ 1.0, 2.0, 3.0 } ) );
    tFix.mSet->mFemModel = &tModel;

    fem::IQI_Factory       tFactory;
    std::shared_ptr< IQI > tIQI = tFactory.create_IQI( fem::IQI_Type::EIGEN_VALUE );
    tIQI->set_name( "IQI" );
    tIQI->set_set_pointer( tFix.mSet );
    tIQI->set_parameters( { { { 0 } } } );

    // vectorial_field_index selects the eigenvalue, no dof_quantity, criterion path
    tIQI->set_output_type_index( 2 );
    CHECK( tFix.global_qi( tIQI ) == Approx( 3.0 ) );

    // the criterion path must not reset the index used by the output path
    Matrix< DDRMat > tQI;
    tIQI->compute_QI( tQI );
    CHECK( tQI( 0 ) == Approx( 3.0 ) );

    // unset: the eigen vector index from function_parameters
    tIQI->set_parameters( { { { 1 } } } );
    tIQI->set_output_type_index( -1 );
    tIQI->compute_QI( tQI );
    CHECK( tQI( 0 ) == Approx( 2.0 ) );

    // out of range
    tIQI->set_output_type_index( 3 );
    CHECK_THROWS( tIQI->compute_QI( tQI ) );

    tFix.mSet->mFemModel = nullptr;
}

//------------------------------------------------------------------------------

TEST_CASE( "IQI_Stabilization_Index", "[moris],[fem],[IQI_Component_Index],[MB-055]" )
{
    Component_Index_Fixture tFix;

    std::shared_ptr< Stabilization_Parameter > tSP = std::make_shared< SP_Two_Values_Test >();

    fem::IQI_Factory       tFactory;
    std::shared_ptr< IQI > tIQI = tFactory.create_IQI( fem::IQI_Type::STABILIZATION );
    tIQI->set_name( "IQI" );
    tIQI->set_set_pointer( tFix.mSet );
    tIQI->set_stabilization_parameter( tSP, "Stabilization" );

    // second SP value, with no dof_quantity (the EXA use: tauC of an NS SP)
    tIQI->set_output_type_index( 1 );
    CHECK( tFix.global_qi( tIQI ) == Approx( 5.0 ) );

    Matrix< DDRMat > tQI;
    tIQI->compute_QI( tQI );
    CHECK( tQI( 0 ) == Approx( 5.0 ) );

    // unset: the first value
    tIQI->set_output_type_index( -1 );
    tIQI->compute_QI( tQI );
    CHECK( tQI( 0 ) == Approx( 2.0 ) );

    // out of range
    tIQI->set_output_type_index( 2 );
    CHECK_THROWS( tIQI->compute_QI( tQI ) );
}
