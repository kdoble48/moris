/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_WRK_TD_Design_Operator.cpp
 *
 */

#include "cl_WRK_TD_Design_Operator.hpp"

#include <algorithm>
#include <cmath>
#include <deque>
#include <unordered_map>

#include "cl_Communication_Tools.hpp"
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"
#include "fn_Parsing_Tools.hpp"
#include "HDF5_Tools.hpp"

#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Vertex_Interpolation.hpp"
#include "cl_MTK_Field.hpp"
#include "cl_MTK_Interpolation_Function_Factory.hpp"
#include "cl_MTK_Interpolation_Function_Base.hpp"
#include "cl_MTK_Integration_Rule.hpp"
#include "cl_MTK_Integration_Coeffs_Base.hpp"

#include "fn_trans.hpp"
#include "op_times.hpp"
#include "op_plus.hpp"

#include "cl_FEM_Field.hpp"
#include "cl_GEN_Geometry_Engine.hpp"

namespace moris::wrk
{
    namespace
    {
        //------------------------------------------------------------------------------

        /** Gauss rule with p_L + 1 points per direction: exact for the mass (degree 2 p_L per direction) */
        mtk::Integration_Order
        integration_order_for_cell(
                mtk::Geometry_Type       aGeometryType,
                mtk::Interpolation_Order aInterpolationOrder )
        {
            uint tOrder = 0;
            switch ( aInterpolationOrder )
            {
                case mtk::Interpolation_Order::LINEAR:
                    tOrder = 1;
                    break;
                case mtk::Interpolation_Order::QUADRATIC:
                    tOrder = 2;
                    break;
                case mtk::Interpolation_Order::CUBIC:
                    tOrder = 3;
                    break;
                default:
                    MORIS_ERROR( false, "TD_Design_Operator - unsupported Lagrange interpolation order of the background cells." );
            }

            if ( aGeometryType == mtk::Geometry_Type::QUAD )
            {
                const mtk::Integration_Order tQuad[ 3 ] = { mtk::Integration_Order::QUAD_2x2, mtk::Integration_Order::QUAD_3x3, mtk::Integration_Order::QUAD_4x4 };
                return tQuad[ tOrder - 1 ];
            }
            if ( aGeometryType == mtk::Geometry_Type::HEX )
            {
                const mtk::Integration_Order tHex[ 3 ] = { mtk::Integration_Order::HEX_2x2x2, mtk::Integration_Order::HEX_3x3x3, mtk::Integration_Order::HEX_4x4x4 };
                return tHex[ tOrder - 1 ];
            }

            MORIS_ERROR( false, "TD_Design_Operator - background cells must be QUAD or HEX (HMR tensor meshes)." );
            return mtk::Integration_Order::UNDEFINED;
        }

        //------------------------------------------------------------------------------

        /** determinant and inverse of a 2 x 2 or 3 x 3 Jacobian */
        real
        invert_jacobian( const Matrix< DDRMat >& aJ, Matrix< DDRMat >& aInverse )
        {
            uint tDim = aJ.n_rows();
            aInverse.set_size( tDim, tDim );

            if ( tDim == 2 )
            {
                real tDet        = aJ( 0, 0 ) * aJ( 1, 1 ) - aJ( 0, 1 ) * aJ( 1, 0 );
                aInverse( 0, 0 ) = aJ( 1, 1 ) / tDet;
                aInverse( 0, 1 ) = -aJ( 0, 1 ) / tDet;
                aInverse( 1, 0 ) = -aJ( 1, 0 ) / tDet;
                aInverse( 1, 1 ) = aJ( 0, 0 ) / tDet;
                return tDet;
            }

            real tDet = aJ( 0, 0 ) * ( aJ( 1, 1 ) * aJ( 2, 2 ) - aJ( 1, 2 ) * aJ( 2, 1 ) )
                      - aJ( 0, 1 ) * ( aJ( 1, 0 ) * aJ( 2, 2 ) - aJ( 1, 2 ) * aJ( 2, 0 ) )
                      + aJ( 0, 2 ) * ( aJ( 1, 0 ) * aJ( 2, 1 ) - aJ( 1, 1 ) * aJ( 2, 0 ) );

            aInverse( 0, 0 ) = ( aJ( 1, 1 ) * aJ( 2, 2 ) - aJ( 1, 2 ) * aJ( 2, 1 ) ) / tDet;
            aInverse( 0, 1 ) = ( aJ( 0, 2 ) * aJ( 2, 1 ) - aJ( 0, 1 ) * aJ( 2, 2 ) ) / tDet;
            aInverse( 0, 2 ) = ( aJ( 0, 1 ) * aJ( 1, 2 ) - aJ( 0, 2 ) * aJ( 1, 1 ) ) / tDet;
            aInverse( 1, 0 ) = ( aJ( 1, 2 ) * aJ( 2, 0 ) - aJ( 1, 0 ) * aJ( 2, 2 ) ) / tDet;
            aInverse( 1, 1 ) = ( aJ( 0, 0 ) * aJ( 2, 2 ) - aJ( 0, 2 ) * aJ( 2, 0 ) ) / tDet;
            aInverse( 1, 2 ) = ( aJ( 0, 2 ) * aJ( 1, 0 ) - aJ( 0, 0 ) * aJ( 1, 2 ) ) / tDet;
            aInverse( 2, 0 ) = ( aJ( 1, 0 ) * aJ( 2, 1 ) - aJ( 1, 1 ) * aJ( 2, 0 ) ) / tDet;
            aInverse( 2, 1 ) = ( aJ( 0, 1 ) * aJ( 2, 0 ) - aJ( 0, 0 ) * aJ( 2, 1 ) ) / tDet;
            aInverse( 2, 2 ) = ( aJ( 0, 0 ) * aJ( 1, 1 ) - aJ( 0, 1 ) * aJ( 1, 0 ) ) / tDet;
            return tDet;
        }

        //------------------------------------------------------------------------------

        /** finds the field with the given label among the FEM fields; errors if missing or not an ELEMENTAL ip_cell field */
        mtk::Field*
        find_elemental_field(
                const Vector< std::shared_ptr< mtk::Field > >& aFields,
                const std::string&                             aName )
        {
            for ( const std::shared_ptr< mtk::Field >& iField : aFields )
            {
                if ( iField->get_label() != aName )
                {
                    continue;
                }

                // the NODAL channel overwrites cut-cell vertices with the last write (IMPL_PLAN risk R5): ELEMENTAL only
                MORIS_ERROR( iField->get_field_entity_type() == mtk::Field_Entity_Type::ELEMENTAL,
                        "TD_Design_Operator - FEM field '%s' must have field_entity_type ELEMENTAL.",
                        aName.c_str() );

                // the per-phase rows of a cut cell add up to g_e only as integrals per interpolation cell
                fem::Field* tFEMField = dynamic_cast< fem::Field* >( iField.get() );
                MORIS_ERROR( tFEMField != nullptr and tFEMField->get_elemental_normalization() == fem::Elemental_Normalization::IP_CELL,
                        "TD_Design_Operator - FEM field '%s' must set field_elemental_normalization = ip_cell.",
                        aName.c_str() );

                return iField.get();
            }

            MORIS_ERROR( false, "TD_Design_Operator - no FEM field named '%s' (FEM FIELDS field_name).", aName.c_str() );
            return nullptr;
        }

        //------------------------------------------------------------------------------

        /** adds the set rows of an ELEMENTAL field to the base cells of their interpolation cells; returns the set rows */
        uint
        add_field_to_base_cells(
                mtk::Field*     aField,
                Vector< real >& aCellValues )
        {
            mtk::Interpolation_Mesh* tFieldMesh = aField->get_mesh_pair().get_interpolation_mesh();
            const Matrix< DDRMat >&  tValues    = aField->get_values();

            MORIS_ERROR( tValues.n_rows() == tFieldMesh->get_num_elems(),
                    "TD_Design_Operator - ELEMENTAL field '%s' has %u rows for %u interpolation cells.",
                    aField->get_label().c_str(),
                    (uint)tValues.n_rows(),
                    tFieldMesh->get_num_elems() );

            uint tNumSet = 0;
            for ( uint iRow = 0; iRow < tValues.n_rows(); iRow++ )
            {
                // rows of cells outside the IQI's sets keep the initial value MORIS_REAL_MIN
                if ( tValues( iRow ) == MORIS_REAL_MIN )
                {
                    continue;
                }

                moris_index tBaseCell = tFieldMesh->get_mtk_cell( iRow ).get_base_cell()->get_index();

                MORIS_ERROR( tBaseCell >= 0 and ( uint ) tBaseCell < aCellValues.size(),
                        "TD_Design_Operator - base cell index %d of field '%s' outside the background mesh (%u cells).",
                        (int)tBaseCell,
                        aField->get_label().c_str(),
                        (uint)aCellValues.size() );

                aCellValues( tBaseCell ) += tValues( iRow );
                tNumSet++;
            }

            return tNumSet;
        }

        //------------------------------------------------------------------------------
    }    // namespace

    //------------------------------------------------------------------------------

    TD_Design_Operator::TD_Design_Operator( const Parameter_List& aParameterList )
            : mADVFieldName( aParameterList.get< std::string >( "adv_field" ) )
            , mCutRatioFieldName( aParameterList.get< std::string >( "cut_ratio_field" ) )
            , mVoidTDValue( aParameterList.get< real >( "void_td_value" ) )
            , mFilterGamma( aParameterList.get< real >( "filter_gamma" ) )
            , mOutputFile( aParameterList.get< std::string >( "output_file" ) )
    {
        MORIS_ERROR( par_size() == 1,
                "TD_Design_Operator - the TD_OPERATOR sublist runs in serial only (the operator lives on one rank)." );

        MORIS_ERROR( not mADVFieldName.empty(), "TD_Design_Operator - adv_field (the level-set GEN design) is not set." );

        std::string tFieldNames = aParameterList.get< std::string >( "td_field_names" );
        MORIS_ERROR( not tFieldNames.empty(), "TD_Design_Operator - td_field_names is not set." );
        mTDFieldNames = split_string( tFieldNames, "," );

        MORIS_ERROR( mFilterGamma >= 0.0, "TD_Design_Operator - filter_gamma must be non-negative." );
    }

    //------------------------------------------------------------------------------

    void
    TD_Design_Operator::build_operators( gen::Geometry_Engine& aGeometryEngine, uint aNumADVs )
    {
        // ADV positions of the design's coefficients; empty if no design has this name
        Vector< sint > tPositions = aGeometryEngine.get_design_adv_positions( mADVFieldName );
        MORIS_ERROR( tPositions.size() > 0,
                "TD_Design_Operator - GEN has no B-spline design named '%s' (adv_field).",
                mADVFieldName.c_str() );

        // with an ADV filter the ADVs are not the coefficients of the level set
        MORIS_ERROR( not aGeometryEngine.has_adv_filter(),
                "TD_Design_Operator - not supported with adv_filter_radius > 0 (the ADVs are not the B-spline coefficients)." );
        MORIS_ERROR( not aGeometryEngine.has_additive_discretized_designs(),
                "TD_Design_Operator - not supported with discretization_mode = \"additive\"." );

        // mesh and discretization index of the design
        Vector< std::shared_ptr< mtk::Field > > tGENFields = aGeometryEngine.get_mtk_fields();

        std::shared_ptr< mtk::Field > tDesignField = nullptr;
        for ( const std::shared_ptr< mtk::Field >& iField : tGENFields )
        {
            if ( iField->get_label() == mADVFieldName )
            {
                tDesignField = iField;
            }
        }
        MORIS_ERROR( tDesignField != nullptr, "TD_Design_Operator - GEN has no field labelled '%s' (adv_field).", mADVFieldName.c_str() );

        sint tDiscretizationIndex = tDesignField->get_discretization_mesh_index();
        MORIS_ERROR( tDiscretizationIndex >= 0, "TD_Design_Operator - design '%s' is not discretized on a B-spline mesh.", mADVFieldName.c_str() );

        this->build_operators(
                tDesignField->get_mesh_pair().get_interpolation_mesh(),
                (uint)tDiscretizationIndex,
                tPositions,
                aNumADVs );
    }

    //------------------------------------------------------------------------------

    void
    TD_Design_Operator::build_operators(
            mtk::Interpolation_Mesh* aMesh,
            uint                     aDiscretizationIndex,
            const Vector< sint >&    aCoefficientADVPositions,
            uint                     aNumADVs )
    {
        Tracer tTracer( "WRK", "TD Design Operator", "Build operators" );

        MORIS_ERROR( par_size() == 1, "TD_Design_Operator::build_operators - serial only." );

        mMesh                = aMesh;
        mDiscretizationIndex = aDiscretizationIndex;
        mNumADVs             = aNumADVs;
        mNumCells            = aMesh->get_num_elems();

        uint tNumCoefficients = aCoefficientADVPositions.size();
        uint tDim             = aMesh->get_spatial_dim();

        // B_A is in the cell's Lagrange space only if the Lagrange order is at least the B-spline order
        MORIS_ERROR( aMesh->get_discretization_order( aDiscretizationIndex ) <= aMesh->get_order(),
                "TD_Design_Operator - design B-spline order %u exceeds the Lagrange order %u of the mesh.",
                aMesh->get_discretization_order( aDiscretizationIndex ),
                aMesh->get_order() );

        //------------------------------------------------------------------------------
        // 1. coefficients present in a T-matrix of an owned cell -> unknowns

        Vector< sint > tCoefficientToUnknown( tNumCoefficients, -1 );
        for ( uint iCell = 0; iCell < mNumCells; iCell++ )
        {
            const mtk::Cell& tCell = aMesh->get_mtk_cell( iCell );
            if ( tCell.get_owner() != par_rank() )
            {
                continue;
            }
            for ( mtk::Vertex* iVertex : tCell.get_vertex_pointers() )
            {
                Matrix< IndexMat > tIndices = iVertex->get_interpolation( aDiscretizationIndex )->get_indices();
                for ( uint iC = 0; iC < tIndices.numel(); iC++ )
                {
                    MORIS_ERROR( tIndices( iC ) >= 0 and ( uint ) tIndices( iC ) < tNumCoefficients,
                            "TD_Design_Operator - T-matrix coefficient index %d outside the design's %u coefficients.",
                            (int)tIndices( iC ),
                            tNumCoefficients );
                    tCoefficientToUnknown( tIndices( iC ) ) = 0;
                }
            }
        }

        mUnknownToADV.clear();
        mNumDroppedCoefficients = 0;
        for ( uint iCoefficient = 0; iCoefficient < tNumCoefficients; iCoefficient++ )
        {
            if ( tCoefficientToUnknown( iCoefficient ) == 0 )
            {
                MORIS_ERROR( aCoefficientADVPositions( iCoefficient ) >= 0 and ( uint ) aCoefficientADVPositions( iCoefficient ) < aNumADVs,
                        "TD_Design_Operator - coefficient %u of the domain is not an ADV.",
                        iCoefficient );

                tCoefficientToUnknown( iCoefficient ) = mUnknownToADV.size();
                mUnknownToADV.push_back( (uint)aCoefficientADVPositions( iCoefficient ) );
            }
            else if ( aCoefficientADVPositions( iCoefficient ) >= 0 )
            {
                // zero T column: the function does not reach the domain; keeping it would make gamma K + M singular
                mNumDroppedCoefficients++;
            }
        }
        uint tNumUnknowns = mUnknownToADV.size();
        MORIS_ERROR( tNumUnknowns > 0, "TD_Design_Operator - no design coefficient reaches the domain." );

        //------------------------------------------------------------------------------
        // 2. cell matrices: B = N T_e at the Gauss points; pattern and values assembled in CSR

        mtk::Interpolation_Function_Factory tFactory;

        Vector< Vector< uint > >   tRowColumns( tNumUnknowns );
        Vector< Vector< uint > >   tCellUnknowns( mNumCells );
        Vector< Matrix< DDRMat > > tCellMass( mNumCells );
        Vector< Matrix< DDRMat > > tCellStiffness( mNumCells );
        Vector< Matrix< DDRMat > > tCellIntegrals( mNumCells );
        mCellMeasures.resize( mNumCells );
        std::fill( mCellMeasures.begin(), mCellMeasures.end(), 0.0 );

        for ( uint iCell = 0; iCell < mNumCells; iCell++ )
        {
            const mtk::Cell& tCell = aMesh->get_mtk_cell( iCell );
            if ( tCell.get_owner() != par_rank() )
            {
                continue;
            }
            uint tCellIndex = tCell.get_index();
            MORIS_ERROR( tCellIndex < mNumCells, "TD_Design_Operator - cell index outside the mesh." );

            Vector< mtk::Vertex* > tVertices = tCell.get_vertex_pointers();
            uint                   tNumNodes = tVertices.size();

            // local unknowns of the cell and the cell T-matrix ( nodes x local unknowns )
            std::unordered_map< uint, uint > tLocal;
            Vector< uint >&                  tUnknowns = tCellUnknowns( tCellIndex );
            for ( mtk::Vertex* iVertex : tVertices )
            {
                Matrix< IndexMat > tIndices = iVertex->get_interpolation( aDiscretizationIndex )->get_indices();
                for ( uint iC = 0; iC < tIndices.numel(); iC++ )
                {
                    uint tUnknown = (uint)tCoefficientToUnknown( tIndices( iC ) );
                    if ( tLocal.find( tUnknown ) == tLocal.end() )
                    {
                        tLocal[ tUnknown ] = tUnknowns.size();
                        tUnknowns.push_back( tUnknown );
                    }
                }
            }
            uint             tNumLocal = tUnknowns.size();
            Matrix< DDRMat > tT( tNumNodes, tNumLocal, 0.0 );
            Matrix< DDRMat > tX( tNumNodes, tDim, 0.0 );
            for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
            {
                const mtk::Vertex_Interpolation* tInterpolation = tVertices( iNode )->get_interpolation( aDiscretizationIndex );
                Matrix< IndexMat >               tIndices       = tInterpolation->get_indices();
                const Matrix< DDRMat >&          tWeights       = *tInterpolation->get_weights();
                for ( uint iC = 0; iC < tIndices.numel(); iC++ )
                {
                    tT( iNode, tLocal[ (uint)tCoefficientToUnknown( tIndices( iC ) ) ] ) += tWeights( iC );
                }
                Matrix< DDRMat > tCoords = tVertices( iNode )->get_coords();
                for ( uint iDim = 0; iDim < tDim; iDim++ )
                {
                    tX( iNode, iDim ) = tCoords( iDim );
                }
            }

            // Lagrange shape functions of the cell and the Gauss rule
            mtk::Interpolation_Function_Base* tShape = tFactory.create_interpolation_function(
                    tCell.get_geometry_type(),
                    mtk::Interpolation_Type::LAGRANGE,
                    tCell.get_interpolation_order() );
            MORIS_ERROR( tShape->get_number_of_bases() == tNumNodes, "TD_Design_Operator - cell vertices do not match its Lagrange basis." );

            mtk::Integration_Rule tRule(
                    tCell.get_geometry_type(),
                    mtk::Integration_Type::GAUSS,
                    integration_order_for_cell( tCell.get_geometry_type(), tCell.get_interpolation_order() ),
                    mtk::Integration_Type::GAUSS,
                    mtk::Integration_Order::BAR_1 );
            std::unique_ptr< mtk::Integration_Coeffs_Base > tCoeffs = tRule.create_space_coeffs();
            Matrix< DDRMat >                                tPoints;
            Matrix< DDRMat >                                tWeights;
            tCoeffs->get_points( tPoints );
            tCoeffs->get_weights( tWeights );

            Matrix< DDRMat > tMe( tNumLocal, tNumLocal, 0.0 );
            Matrix< DDRMat > tKe( tNumLocal, tNumLocal, 0.0 );
            Matrix< DDRMat > tWe( tNumLocal, 1, 0.0 );
            Matrix< DDRMat > tN;
            Matrix< DDRMat > tdNdXi;
            Matrix< DDRMat > tJInverse;
            for ( uint iPoint = 0; iPoint < tWeights.numel(); iPoint++ )
            {
                Matrix< DDRMat > tXi = tPoints.get_column( iPoint );
                tShape->eval_N( tXi, tN );
                tShape->eval_dNdXi( tXi, tdNdXi );

                Matrix< DDRMat > tJ   = tdNdXi * tX;    // ( dim x dim ), J_ij = dx_j / dxi_i
                real             tDet = invert_jacobian( tJ, tJInverse );
                MORIS_ERROR( tDet > 0.0, "TD_Design_Operator - non-positive Jacobian on cell %u.", tCellIndex );

                real             tW    = tWeights( iPoint ) * tDet;
                Matrix< DDRMat > tB    = tN * tT;                        // ( 1 x local )
                Matrix< DDRMat > tdBdx = tJInverse * ( tdNdXi * tT );    // ( dim x local )

                tMe = tMe + tW * ( trans( tB ) * tB );
                tKe = tKe + tW * ( trans( tdBdx ) * tdBdx );
                tWe = tWe + tW * trans( tB );
                mCellMeasures( tCellIndex ) += tW;
            }
            delete tShape;

            tCellMass( tCellIndex )      = tMe;
            tCellStiffness( tCellIndex ) = tKe;
            tCellIntegrals( tCellIndex ) = tWe;

            for ( uint iA = 0; iA < tNumLocal; iA++ )
            {
                for ( uint iB = 0; iB < tNumLocal; iB++ )
                {
                    tRowColumns( tUnknowns( iA ) ).push_back( tUnknowns( iB ) );
                }
            }
        }

        // CSR pattern
        mRowOffsets.resize( tNumUnknowns + 1 );
        mRowOffsets( 0 ) = 0;
        for ( uint iRow = 0; iRow < tNumUnknowns; iRow++ )
        {
            Vector< uint >& tColumns = tRowColumns( iRow );
            std::sort( tColumns.begin(), tColumns.end() );
            tColumns.data().erase( std::unique( tColumns.begin(), tColumns.end() ), tColumns.end() );
            mRowOffsets( iRow + 1 ) = mRowOffsets( iRow ) + tColumns.size();
        }
        mColumns.resize( mRowOffsets( tNumUnknowns ) );
        for ( uint iRow = 0; iRow < tNumUnknowns; iRow++ )
        {
            std::copy( tRowColumns( iRow ).begin(), tRowColumns( iRow ).end(), mColumns.begin() + mRowOffsets( iRow ) );
        }
        tRowColumns.clear();

        // values, and the per-cell integrals in CSR by cell
        mMassValues.resize( mColumns.size() );
        mStiffnessValues.resize( mColumns.size() );
        std::fill( mMassValues.begin(), mMassValues.end(), 0.0 );
        std::fill( mStiffnessValues.begin(), mStiffnessValues.end(), 0.0 );
        mCellOffsets.resize( mNumCells + 1 );
        mCellOffsets( 0 ) = 0;
        mCellUnknowns.clear();
        mCellIntegrals.clear();
        mDomainMeasure = 0.0;
        for ( uint iCell = 0; iCell < mNumCells; iCell++ )
        {
            const Vector< uint >& tUnknowns = tCellUnknowns( iCell );
            for ( uint iA = 0; iA < tUnknowns.size(); iA++ )
            {
                uint tRow = tUnknowns( iA );
                for ( uint iB = 0; iB < tUnknowns.size(); iB++ )
                {
                    const uint* tFirst = mColumns.memptr() + mRowOffsets( tRow );
                    const uint* tLast  = mColumns.memptr() + mRowOffsets( tRow + 1 );
                    uint        tPos   = std::lower_bound( tFirst, tLast, tUnknowns( iB ) ) - mColumns.memptr();
                    mMassValues( tPos ) += tCellMass( iCell )( iA, iB );
                    mStiffnessValues( tPos ) += tCellStiffness( iCell )( iA, iB );
                }
                mCellUnknowns.push_back( tRow );
                mCellIntegrals.push_back( tCellIntegrals( iCell )( iA ) );
            }
            mCellOffsets( iCell + 1 ) = mCellUnknowns.size();
            mDomainMeasure += mCellMeasures( iCell );
        }

        //------------------------------------------------------------------------------
        // 3. sanity checks: partition of unity, consistent row sums, constants in the kernel of K, no zero rows

        Vector< real > tIntegrals( tNumUnknowns, 0.0 );
        for ( uint iEntry = 0; iEntry < mCellUnknowns.size(); iEntry++ )
        {
            tIntegrals( mCellUnknowns( iEntry ) ) += mCellIntegrals( iEntry );
        }

        mMassSum          = 0.0;
        mMinDiagonal      = MORIS_REAL_MAX;
        real tMaxRowError = 0.0;
        real tMaxKRowSum  = 0.0;
        real tMaxKDiag    = 0.0;
        real tMaxIntegral = 0.0;
        for ( uint iRow = 0; iRow < tNumUnknowns; iRow++ )
        {
            real tMassRowSum      = 0.0;
            real tStiffnessRowSum = 0.0;
            for ( uint iEntry = mRowOffsets( iRow ); iEntry < mRowOffsets( iRow + 1 ); iEntry++ )
            {
                tMassRowSum += mMassValues( iEntry );
                tStiffnessRowSum += mStiffnessValues( iEntry );
                if ( mColumns( iEntry ) == iRow )
                {
                    mMinDiagonal = std::min( mMinDiagonal, mFilterGamma * mStiffnessValues( iEntry ) + mMassValues( iEntry ) );
                    tMaxKDiag    = std::max( tMaxKDiag, std::abs( mStiffnessValues( iEntry ) ) );
                }
            }
            mMassSum += tMassRowSum;
            tMaxRowError = std::max( tMaxRowError, std::abs( tMassRowSum - tIntegrals( iRow ) ) );
            tMaxKRowSum  = std::max( tMaxKRowSum, std::abs( tStiffnessRowSum ) );
            tMaxIntegral = std::max( tMaxIntegral, std::abs( tIntegrals( iRow ) ) );
        }

        MORIS_ERROR( std::abs( mMassSum - mDomainMeasure ) <= 1.0e-10 * mDomainMeasure,
                "TD_Design_Operator - sum of the mass matrix %.15e differs from the domain measure %.15e: the B-splines do not "
                "form a partition of unity on the domain (padding or missing coefficients).",
                mMassSum,
                mDomainMeasure );
        MORIS_ERROR( tMaxRowError <= 1.0e-10 * tMaxIntegral,
                "TD_Design_Operator - row sums of M differ from int B_A by %.3e.",
                tMaxRowError );
        MORIS_ERROR( tMaxKRowSum <= 1.0e-10 * tMaxKDiag,
                "TD_Design_Operator - K 1 = %.3e, constants are not in the kernel of the Laplacian.",
                tMaxKRowSum );
        MORIS_ERROR( mMinDiagonal > 0.0,
                "TD_Design_Operator - gamma K + M has a zero row (min diagonal %.3e).",
                mMinDiagonal );

        //------------------------------------------------------------------------------
        // 4. factor gamma K + M (Cholesky: errors unless SPD)

        this->factorize();

        mNumBuilds++;

        MORIS_LOG_INFO(
                "TD_Design_Operator: design '%s' (mesh index %u): %u unknowns of %u coefficients (%u dropped, zero T column), "
                "%u ADVs, %u cells, nnz %u, envelope %u",
                mADVFieldName.c_str(),
                aDiscretizationIndex,
                tNumUnknowns,
                tNumCoefficients,
                mNumDroppedCoefficients,
                aNumADVs,
                mNumCells,
                (uint)mColumns.size(),
                (uint)mFactor.size() );
        MORIS_LOG_INFO(
                "TD_Design_Operator: sum M = %.15e, |D| = %.15e, max |K 1| = %.3e, min diag( gamma K + M ) = %.3e, "
                "min Cholesky pivot = %.3e, gamma = %.3e",
                mMassSum,
                mDomainMeasure,
                tMaxKRowSum,
                mMinDiagonal,
                mMinPivot,
                mFilterGamma );
    }

    //------------------------------------------------------------------------------

    void
    TD_Design_Operator::factorize()
    {
        uint tNumUnknowns = mUnknownToADV.size();

        // reverse Cuthill-McKee ordering of the matrix graph (keeps the envelope small)
        Vector< uint > tDegree( tNumUnknowns );
        for ( uint iRow = 0; iRow < tNumUnknowns; iRow++ )
        {
            tDegree( iRow ) = mRowOffsets( iRow + 1 ) - mRowOffsets( iRow );
        }

        Vector< uint > tOrder;
        Vector< bool > tVisited( tNumUnknowns, false );
        Vector< uint > tNeighbors;
        for ( uint iStart = 0; iStart < tNumUnknowns; iStart++ )
        {
            if ( tVisited( iStart ) )
            {
                continue;
            }

            // start the component at its unvisited node of minimum degree
            uint tRoot = iStart;
            {
                std::deque< uint > tQueue( 1, iStart );
                Vector< bool >     tSeen( tNumUnknowns, false );
                tSeen( iStart ) = true;
                while ( not tQueue.empty() )
                {
                    uint tNode = tQueue.front();
                    tQueue.pop_front();
                    if ( tDegree( tNode ) < tDegree( tRoot ) )
                    {
                        tRoot = tNode;
                    }
                    for ( uint iEntry = mRowOffsets( tNode ); iEntry < mRowOffsets( tNode + 1 ); iEntry++ )
                    {
                        if ( not tSeen( mColumns( iEntry ) ) and not tVisited( mColumns( iEntry ) ) )
                        {
                            tSeen( mColumns( iEntry ) ) = true;
                            tQueue.push_back( mColumns( iEntry ) );
                        }
                    }
                }
            }

            std::deque< uint > tQueue( 1, tRoot );
            tVisited( tRoot ) = true;
            while ( not tQueue.empty() )
            {
                uint tNode = tQueue.front();
                tQueue.pop_front();
                tOrder.push_back( tNode );

                tNeighbors.clear();
                for ( uint iEntry = mRowOffsets( tNode ); iEntry < mRowOffsets( tNode + 1 ); iEntry++ )
                {
                    if ( not tVisited( mColumns( iEntry ) ) )
                    {
                        tVisited( mColumns( iEntry ) ) = true;
                        tNeighbors.push_back( mColumns( iEntry ) );
                    }
                }
                std::stable_sort( tNeighbors.begin(), tNeighbors.end(), [ &tDegree ]( uint aFirst, uint aSecond ) { return tDegree( aFirst ) < tDegree( aSecond ); } );
                for ( uint iNeighbor : tNeighbors )
                {
                    tQueue.push_back( iNeighbor );
                }
            }
        }
        std::reverse( tOrder.begin(), tOrder.end() );

        mPermutation = tOrder;
        Vector< uint > tInverse( tNumUnknowns );
        for ( uint iNew = 0; iNew < tNumUnknowns; iNew++ )
        {
            tInverse( mPermutation( iNew ) ) = iNew;
        }

        // envelope of the permuted lower triangle
        mEnvelopeFirst.resize( tNumUnknowns );
        mEnvelopeOffsets.resize( tNumUnknowns + 1 );
        mEnvelopeOffsets( 0 ) = 0;
        for ( uint iNew = 0; iNew < tNumUnknowns; iNew++ )
        {
            uint tOld   = mPermutation( iNew );
            uint tFirst = iNew;
            for ( uint iEntry = mRowOffsets( tOld ); iEntry < mRowOffsets( tOld + 1 ); iEntry++ )
            {
                tFirst = std::min( tFirst, tInverse( mColumns( iEntry ) ) );
            }
            mEnvelopeFirst( iNew )       = tFirst;
            mEnvelopeOffsets( iNew + 1 ) = mEnvelopeOffsets( iNew ) + ( iNew - tFirst + 1 );
        }

        // scatter gamma K + M into the envelope
        mFactor.resize( mEnvelopeOffsets( tNumUnknowns ) );
        std::fill( mFactor.begin(), mFactor.end(), 0.0 );
        for ( uint iNew = 0; iNew < tNumUnknowns; iNew++ )
        {
            uint tOld = mPermutation( iNew );
            for ( uint iEntry = mRowOffsets( tOld ); iEntry < mRowOffsets( tOld + 1 ); iEntry++ )
            {
                uint tColumn = tInverse( mColumns( iEntry ) );
                if ( tColumn <= iNew )
                {
                    mFactor( mEnvelopeOffsets( iNew ) + tColumn - mEnvelopeFirst( iNew ) ) =
                            mFilterGamma * mStiffnessValues( iEntry ) + mMassValues( iEntry );
                }
            }
        }

        // row-wise envelope Cholesky A = L L^T (fill stays inside the envelope)
        mMinPivot = MORIS_REAL_MAX;
        for ( uint iRow = 0; iRow < tNumUnknowns; iRow++ )
        {
            uint  tFirstI = mEnvelopeFirst( iRow );
            real* tRowI   = mFactor.memptr() + mEnvelopeOffsets( iRow );    // tRowI[ k - tFirstI ] = L( iRow, k )
            for ( uint iCol = tFirstI; iCol < iRow; iCol++ )
            {
                uint        tFirstJ = mEnvelopeFirst( iCol );
                const real* tRowJ   = mFactor.memptr() + mEnvelopeOffsets( iCol );
                real        tSum    = tRowI[ iCol - tFirstI ];
                for ( uint k = std::max( tFirstI, tFirstJ ); k < iCol; k++ )
                {
                    tSum -= tRowI[ k - tFirstI ] * tRowJ[ k - tFirstJ ];
                }
                tRowI[ iCol - tFirstI ] = tSum / tRowJ[ iCol - tFirstJ ];
            }
            real tPivot = tRowI[ iRow - tFirstI ];
            for ( uint k = tFirstI; k < iRow; k++ )
            {
                tPivot -= tRowI[ k - tFirstI ] * tRowI[ k - tFirstI ];
            }
            MORIS_ERROR( tPivot > 0.0,
                    "TD_Design_Operator - gamma K + M is not positive definite (Cholesky pivot %.3e at row %u).",
                    tPivot,
                    iRow );
            mMinPivot               = std::min( mMinPivot, tPivot );
            tRowI[ iRow - tFirstI ] = std::sqrt( tPivot );
        }
    }

    //------------------------------------------------------------------------------

    Vector< real >
    TD_Design_Operator::to_unknowns( const Vector< real >& aADVVector ) const
    {
        MORIS_ERROR( aADVVector.size() == mNumADVs,
                "TD_Design_Operator - vector of length %u, expected the %u ADVs.",
                (uint)aADVVector.size(),
                mNumADVs );

        Vector< real > tUnknowns( mUnknownToADV.size() );
        for ( uint iUnknown = 0; iUnknown < mUnknownToADV.size(); iUnknown++ )
        {
            tUnknowns( iUnknown ) = aADVVector( mUnknownToADV( iUnknown ) );
        }
        return tUnknowns;
    }

    //------------------------------------------------------------------------------

    Vector< real >
    TD_Design_Operator::to_adv_vector( const Vector< real >& aUnknownVector ) const
    {
        Vector< real > tADVVector( mNumADVs, 0.0 );
        for ( uint iUnknown = 0; iUnknown < mUnknownToADV.size(); iUnknown++ )
        {
            tADVVector( mUnknownToADV( iUnknown ) ) = aUnknownVector( iUnknown );
        }
        return tADVVector;
    }

    //------------------------------------------------------------------------------

    Vector< real >
    TD_Design_Operator::csr_product( const Vector< real >& aValues, const Vector< real >& aUnknownVector ) const
    {
        Vector< real > tResult( mUnknownToADV.size(), 0.0 );
        for ( uint iRow = 0; iRow < mUnknownToADV.size(); iRow++ )
        {
            real tSum = 0.0;
            for ( uint iEntry = mRowOffsets( iRow ); iEntry < mRowOffsets( iRow + 1 ); iEntry++ )
            {
                tSum += aValues( iEntry ) * aUnknownVector( mColumns( iEntry ) );
            }
            tResult( iRow ) = tSum;
        }
        return tResult;
    }

    //------------------------------------------------------------------------------

    Vector< real >
    TD_Design_Operator::assemble_load( const Vector< real >& aCellValues ) const
    {
        MORIS_ERROR( aCellValues.size() == mNumCells,
                "TD_Design_Operator::assemble_load - %u cell values for %u cells.",
                (uint)aCellValues.size(),
                mNumCells );

        Vector< real > tLoad( mUnknownToADV.size(), 0.0 );
        for ( uint iCell = 0; iCell < mNumCells; iCell++ )
        {
            for ( uint iEntry = mCellOffsets( iCell ); iEntry < mCellOffsets( iCell + 1 ); iEntry++ )
            {
                tLoad( mCellUnknowns( iEntry ) ) += aCellValues( iCell ) * mCellIntegrals( iEntry );
            }
        }
        return this->to_adv_vector( tLoad );
    }

    //------------------------------------------------------------------------------

    Vector< real >
    TD_Design_Operator::filter( const Vector< real >& aLoad ) const
    {
        Vector< real > tRight = this->to_unknowns( aLoad );
        uint           tNum   = tRight.size();

        // permuted right-hand side
        Vector< real > tY( tNum );
        for ( uint iNew = 0; iNew < tNum; iNew++ )
        {
            tY( iNew ) = tRight( mPermutation( iNew ) );
        }

        // L y = b
        for ( uint iRow = 0; iRow < tNum; iRow++ )
        {
            uint        tFirst = mEnvelopeFirst( iRow );
            const real* tRowI  = mFactor.memptr() + mEnvelopeOffsets( iRow );
            real        tSum   = tY( iRow );
            for ( uint k = tFirst; k < iRow; k++ )
            {
                tSum -= tRowI[ k - tFirst ] * tY( k );
            }
            tY( iRow ) = tSum / tRowI[ iRow - tFirst ];
        }

        // L^T x = y
        for ( uint iRow = tNum; iRow-- > 0; )
        {
            uint        tFirst = mEnvelopeFirst( iRow );
            const real* tRowI  = mFactor.memptr() + mEnvelopeOffsets( iRow );
            tY( iRow ) /= tRowI[ iRow - tFirst ];
            for ( uint k = tFirst; k < iRow; k++ )
            {
                tY( k ) -= tRowI[ k - tFirst ] * tY( iRow );
            }
        }

        Vector< real > tSolution( tNum );
        for ( uint iNew = 0; iNew < tNum; iNew++ )
        {
            tSolution( mPermutation( iNew ) ) = tY( iNew );
        }
        return this->to_adv_vector( tSolution );
    }

    //------------------------------------------------------------------------------

    real
    TD_Design_Operator::mass_product( const Vector< real >& aFirst, const Vector< real >& aSecond ) const
    {
        Vector< real > tFirst  = this->to_unknowns( aFirst );
        Vector< real > tSecond = this->csr_product( mMassValues, this->to_unknowns( aSecond ) );

        real tProduct = 0.0;
        for ( uint iUnknown = 0; iUnknown < tFirst.size(); iUnknown++ )
        {
            tProduct += tFirst( iUnknown ) * tSecond( iUnknown );
        }
        return tProduct;
    }

    //------------------------------------------------------------------------------

    Vector< real >
    TD_Design_Operator::mass_times( const Vector< real >& aVector ) const
    {
        return this->to_adv_vector( this->csr_product( mMassValues, this->to_unknowns( aVector ) ) );
    }

    //------------------------------------------------------------------------------

    Vector< real >
    TD_Design_Operator::stiffness_times( const Vector< real >& aVector ) const
    {
        return this->to_adv_vector( this->csr_product( mStiffnessValues, this->to_unknowns( aVector ) ) );
    }

    //------------------------------------------------------------------------------

    Vector< real >
    TD_Design_Operator::get_coefficient_integrals() const
    {
        Vector< real > tIntegrals( mUnknownToADV.size(), 0.0 );
        for ( uint iEntry = 0; iEntry < mCellUnknowns.size(); iEntry++ )
        {
            tIntegrals( mCellUnknowns( iEntry ) ) += mCellIntegrals( iEntry );
        }
        return this->to_adv_vector( tIntegrals );
    }

    //------------------------------------------------------------------------------

    Vector< real >
    TD_Design_Operator::collect_cell_values( const Vector< std::shared_ptr< mtk::Field > >& aFields ) const
    {
        Vector< real > tCellValues( mNumCells, 0.0 );

        uint tNumSetTotal = 0;
        for ( const std::string& iName : mTDFieldNames )
        {
            mtk::Field* tField  = find_elemental_field( aFields, iName );
            uint        tNumSet = add_field_to_base_cells( tField, tCellValues );
            tNumSetTotal += tNumSet;

            // a phase that is absent from the design (e.g. no void at a full-solid start) leaves its field empty
            if ( tNumSet == 0 )
            {
                MORIS_LOG_INFO( "TD_Design_Operator: FEM field '%s' has no values in this evaluation (phase absent?).", iName.c_str() );
            }
        }

        // fields are populated with the solution output (output criterion of the time solver)
        MORIS_ERROR( tNumSetTotal > 0,
                "TD_Design_Operator - the td_field_names fields have no values: FEM fields are populated at the solution output "
                "(VIS and the time solver's output criterion)." );

        if ( not mCutRatioFieldName.empty() )
        {
            Vector< real > tCutRatio( mNumCells, 0.0 );
            add_field_to_base_cells( find_elemental_field( aFields, mCutRatioFieldName ), tCutRatio );

            for ( uint iCell = 0; iCell < mNumCells; iCell++ )
            {
                tCellValues( iCell ) += ( 1.0 - tCutRatio( iCell ) ) * mVoidTDValue;
            }
        }

        return tCellValues;
    }

    //------------------------------------------------------------------------------

    void
    TD_Design_Operator::perform(
            gen::Geometry_Engine&                          aGeometryEngine,
            const Vector< std::shared_ptr< mtk::Field > >& aFields,
            const Vector< real >&                          aADVs )
    {
        Tracer tTracer( "WRK", "TD Design Operator", "Perform" );

        // (re)build on the first evaluation and whenever the design mesh or the ADV vector changed
        bool tRebuild = mNumBuilds == 0 or aADVs.size() != mNumADVs;
        if ( not tRebuild )
        {
            Vector< std::shared_ptr< mtk::Field > > tGENFields = aGeometryEngine.get_mtk_fields();
            for ( const std::shared_ptr< mtk::Field >& iField : tGENFields )
            {
                if ( iField->get_label() == mADVFieldName )
                {
                    tRebuild = iField->get_mesh_pair().get_interpolation_mesh() != mMesh;
                }
            }
        }
        if ( tRebuild )
        {
            this->build_operators( aGeometryEngine, aADVs.size() );
        }

        mCellValues    = this->collect_cell_values( aFields );
        mLoad          = this->assemble_load( mCellValues );
        mDirection     = this->filter( mLoad );
        mEvaluatedADVs = aADVs;
        mNumEvaluations++;

        real tIntegral = 0.0;
        for ( uint iCell = 0; iCell < mNumCells; iCell++ )
        {
            tIntegral += mCellValues( iCell ) * mCellMeasures( iCell );
        }
        real tLoadSum = 0.0;
        for ( real iLoad : mLoad )
        {
            tLoadSum += iLoad;
        }
        MORIS_LOG_INFO( "TD_Design_Operator: evaluation %u, int g = %.15e, sum b = %.15e, |g~|_M = %.15e",
                mNumEvaluations,
                tIntegral,
                tLoadSum,
                std::sqrt( std::max( 0.0, this->mass_product( mDirection, mDirection ) ) ) );

        if ( not mOutputFile.empty() )
        {
            this->write_output();
        }
    }

    //------------------------------------------------------------------------------

    void
    TD_Design_Operator::write_output() const
    {
        std::string tBase = mOutputFile;
        if ( tBase.size() > 5 and tBase.substr( tBase.size() - 5 ) == ".hdf5" )
        {
            tBase = tBase.substr( 0, tBase.size() - 5 );
        }
        std::string tFileName = tBase + "_" + std::to_string( mNumEvaluations ) + ".hdf5";

        uint             tDim = mMesh->get_spatial_dim();
        Matrix< DDRMat > tCentroids( mNumCells, tDim, 0.0 );
        for ( uint iCell = 0; iCell < mNumCells; iCell++ )
        {
            const mtk::Cell&       tCell     = mMesh->get_mtk_cell( iCell );
            Vector< mtk::Vertex* > tVertices = tCell.get_vertex_pointers();
            Matrix< DDRMat >       tCenter( 1, tDim, 0.0 );
            uint                   tNumCorners = tDim == 2 ? 4 : 8;
            for ( uint iVertex = 0; iVertex < tNumCorners; iVertex++ )
            {
                tCenter = tCenter + tVertices( iVertex )->get_coords();
            }
            for ( uint iDim = 0; iDim < tDim; iDim++ )
            {
                tCentroids( tCell.get_index(), iDim ) = tCenter( iDim ) / tNumCorners;
            }
        }

        hid_t  tFile   = create_hdf5_file( tFileName, false );
        herr_t tStatus = 0;
        save_vector_to_hdf5_file( tFile, "CellValues", mCellValues.data(), tStatus );
        save_vector_to_hdf5_file( tFile, "CellMeasures", mCellMeasures.data(), tStatus );
        save_matrix_to_hdf5_file( tFile, "CellCentroids", tCentroids, tStatus );
        save_vector_to_hdf5_file( tFile, "Load", mLoad.data(), tStatus );
        save_vector_to_hdf5_file( tFile, "Direction", mDirection.data(), tStatus );
        save_vector_to_hdf5_file( tFile, "ADVs", mEvaluatedADVs.data(), tStatus );
        close_hdf5_file( tFile );

        MORIS_LOG_INFO( "TD_Design_Operator: wrote %s", tFileName.c_str() );
    }

    //------------------------------------------------------------------------------

}    // namespace moris::wrk
