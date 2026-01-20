/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_WRK_Reinitialize_Performer.cpp
 *
 */

#include "cl_WRK_Reinitialize_Performer.hpp"

#include <memory>
#include "cl_Matrix.hpp"
#include "cl_Tracer.hpp"
#include "cl_Logger.hpp"
#include "cl_MTK_Mesh_Pair.hpp"

#include "cl_MTK_Field.hpp"
#include "cl_MTK_Field_Discrete.hpp"
#include "cl_MTK_Mapper.hpp"
#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_MTK_Reader_Exodus.hpp"

#include "cl_MSI_Dof_Type_Enums.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_HMR_File.hpp"
#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR_Mesh_Integration.hpp"
#include "HMR_Globals.hpp"

#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_MDL_Model.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_MSI_Solver_Interface.hpp"

#include "cl_SOL_Dist_Vector.hpp"
#include "cl_SOL_Dist_Map.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Warehouse.hpp"
#include "cl_MSI_Equation_Model.hpp"
#include "fn_sort.hpp"

// XTK includes for interface-based SDF
#include "cl_XTK_Cut_Integration_Mesh.hpp"

// SDF includes for distance computation
#include "cl_SDF_From_Interface.hpp"

namespace moris::wrk
{
    //------------------------------------------------------------------------------

    Reinitialize_Performer::Reinitialize_Performer( const std::shared_ptr< Library_IO >& aLibrary )
            : mLibrary( aLibrary )
    {
        // get the parameter lists
        Module_Parameter_Lists tMORISParameterList = aLibrary->get_parameters_for_module( Module_Type::MORISGENERAL );
        Module_Parameter_Lists tMSIParameterList   = aLibrary->get_parameters_for_module( Module_Type::MSI );

        mAdofMeshIndex = tMSIParameterList( 0 )( 0 ).get< moris::sint >( tMORISParameterList( 2 )( 0 ).get< std::string >( "dof_type" ) );

        // get the adv field name that wll be reinitialized
        mADVFiledName = tMORISParameterList( 2 )( 0 ).get< std::string >( "adv_field" );

        // get msi string to dof type map
        moris::map< std::string, MSI::Dof_Type > tMSIDofTypeMap =
                moris::MSI::get_msi_dof_type_map();

        // get the quantity dof type from parameter list
        string_to_vector(
                tMORISParameterList( 2 )( 0 ).get< std::string >( "dof_type" ),
                mDofTypes,
                tMSIDofTypeMap );

        mReinitializationFrequency = tMORISParameterList( 2 )( 0 ).get< sint >( "reinitialization_frequency" );

        // get the mesh output info
        mOutputMeshFile = tMORISParameterList( 2 )( 0 ).get< std::string >( "output_mesh_file" );
        mTimeOffset     = tMORISParameterList( 2 )( 0 ).get< real >( "time_offset" );
    }

    //------------------------------------------------------------------------------

    void
    Reinitialize_Performer::perform(
            Vector< std::shared_ptr< hmr::HMR > >&             aHMRPerformers,
            Vector< std::shared_ptr< gen::Geometry_Engine > >& aGENPerformer,
            Vector< std::shared_ptr< mtk::Mesh_Manager > >&    aMTKPerformer,
            Vector< std::shared_ptr< mdl::Model > >            aMDLPerformer )
    {
        // Tracer to trace the time
        Tracer tTracer( "WRK", "Reinitialize ADVs", "Perform Reinitialize" );

        // initialize and populate the fields
        Vector< std::shared_ptr< mtk::Field > > tGENFields;
        tGENFields.append( aGENPerformer( 0 )->get_mtk_fields() );

        // find the index of the desired adv field that will be reinitialized
        auto itr = std::find_if( tGENFields.begin(), tGENFields.end(), [ & ]( std::shared_ptr< mtk::Field > const & aFiled )    //
                { return aFiled->get_label() == mADVFiledName; } );

        // find the index of the adv field
        moris_index tADVFieldIndex = std::distance( tGENFields.begin(), itr );

        // get the the adv discretization mesh index
        uint tDiscretizationMeshIndex = ( *itr )->get_discretization_mesh_index();

        // get interpolation mesh from mesh pair
        moris::mtk::Mesh* tTargetMesh = ( *itr )->get_mesh_pair().get_interpolation_mesh();

        // get the solution field and get a matrix of the solutions
        // generate a cell containing the indices of the bspline coefficients
        // since indices are consecutive and they start from 0
        Vector< moris_index > tLocalCoeffIndices( tTargetMesh->get_num_entities( mtk::EntityRank::BSPLINE ) );
        std::iota( tLocalCoeffIndices.begin(), tLocalCoeffIndices.end(), 0 );

        moris::sol::Dist_Vector* tPartialSolutionVector = aMDLPerformer( 0 )->get_solver_interface()->get_solution_vector( mDofTypes, tLocalCoeffIndices );
        tPartialSolutionVector->extract_copy( mCoefficients );

        // delete the pointer as it is not needed anymore
        delete tPartialSolutionVector;

        // create field object for this mesh ,the discretization index is zero as there is only one discretization in the newly constructed IP mesh
        std::shared_ptr< mtk::Field_Discrete > tFieldSource = std::make_shared< mtk::Field_Discrete >( aMTKPerformer( 0 )->get_mesh_pair( 0 ), mAdofMeshIndex );

        // unlock fields and set the coeff
        tFieldSource->unlock_field();
        tFieldSource->set_coefficients( mCoefficients );

        // compute the nodal values based on the coeff
        tFieldSource->compute_nodal_values();

        // create field object for this mesh ,the discretization index is zero as there is only one discretization in the newly constructed IP mesh
        std::shared_ptr< mtk::Field_Discrete > tFieldTarget = std::make_shared< mtk::Field_Discrete >( ( *itr )->get_mesh_pair(), tDiscretizationMeshIndex );
        tFieldTarget->set_label( mADVFiledName );

        // set the nodal values
        tFieldTarget->unlock_field();
        tFieldTarget->set_values( tFieldSource->get_values() );

        // invoke the mapper and map to the target field
        mtk::Mapper tMapper;
        tFieldTarget->unlock_field();
        tMapper.map_input_field_to_output_field_2( tFieldTarget.get() );

        // compute the nodal value
        tFieldTarget->compute_nodal_values();

        // get the coefficents and store them
        mCoefficients = tFieldTarget->get_coefficients();

        // clip the values and
        this->impose_upper_lower_bound( aGENPerformer, tFieldTarget.get() );

        // replace the newly constructed field
        tGENFields( tADVFieldIndex ) = tFieldTarget;

        // store the fields
        mMTKFields = tGENFields;

        // output the fields if asked
        if ( mOutputMeshFile != "" )
        {
            this->output_fields( tFieldTarget.get(), tFieldSource.get(), mOutputMeshFile );
        }
    }

    //------------------------------------------------------------------------------

    moris::sint
    Reinitialize_Performer::get_reinitialization_frequency() const
    {
        return mReinitializationFrequency;
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat > const &
    Reinitialize_Performer::get_coefficients() const
    {
        return mCoefficients;
    }
    //------------------------------------------------------------------------------
    void
    Reinitialize_Performer::impose_upper_lower_bound( Vector< std::shared_ptr< gen::Geometry_Engine > >& aGENPerformer, mtk::Field* aField )
    {
        // lower bound and upper bound are defined on proc 0 and they need to be communicated to other
        // Note:  we make an assumption that all the lower bounds and upper bounds are equal

        // initialize  the upper and lower abound
        moris::real tLowerBound;
        moris::real tUpperBound;

        // assign the values on processor 0
        if ( par_rank() == 0 )
        {
            tLowerBound = aGENPerformer( 0 )->get_lower_bounds()( 0 );
            tUpperBound = aGENPerformer( 0 )->get_upper_bounds()( 0 );
        }    // Bcast the values to other processeors
        MPI_Bcast( &tLowerBound, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
        MPI_Bcast( &tUpperBound, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

        // clip the values of the adv
        for ( uint iADV = 0; iADV < mCoefficients.numel(); iADV++ )
        {
            mCoefficients( iADV ) = std::max( tLowerBound, std::min( mCoefficients( iADV ), tUpperBound ) );
        }

        // update the field based on the newly clipped coeff
        aField->unlock_field();
        aField->set_coefficients( mCoefficients );
        aField->compute_nodal_values();
    }

    //------------------------------------------------------------------------------

    Vector< std::shared_ptr< mtk::Field > >
    Reinitialize_Performer::get_mtk_fields() const
    {
        return mMTKFields;
    }

    //------------------------------------------------------------------------------

    void
    Reinitialize_Performer::output_fields( mtk::Field* aTarget, mtk::Field* aSource, std::string aExoFileName ) const
    {
        //
        Tracer tTracer( "WRK", "Reinitialize ADVs", "Outputting Fields" );
        // time shift
        real tTimeShift = 0.0;

        if ( mTimeOffset > 0 )
        {
            // get optimization iteration
            uint tOptIter = gLogger.get_opt_iteration();

            // set name
            std::string tOptIterStrg = std::to_string( tOptIter );
            aExoFileName += ".e-s." + std::string( 4 - tOptIterStrg.length(), '0' ) + tOptIterStrg;

            // determine time shift
            tTimeShift = tOptIter * mTimeOffset;
        }
        // call the lagrange mesh
        moris::mtk::Mesh* tMesh = aTarget->get_mesh_pair().get_interpolation_mesh();
        // Write mesh
        mtk::Writer_Exodus tWriter( tMesh );
        tWriter.write_mesh( "./", aExoFileName, "./", "gen_temp.exo" );

        // write time to file
        tWriter.set_time( tTimeShift );

        // Set nodal fields based on field names
        Vector< std::string > tNodalFieldNames = { "Mapped_Field", "Original_Field" };
        tWriter.set_nodal_fields( tNodalFieldNames );

        // Create field on mesh
        tWriter.write_nodal_field( tNodalFieldNames( 0 ), aTarget->get_values() );
        tWriter.write_nodal_field( tNodalFieldNames( 1 ), aSource->get_values() );

        // Finalize
        tWriter.close_file( true );
    }

    //------------------------------------------------------------------------------

    void
    Reinitialize_Performer::set_reinit_mode( ReinitMode aMode )
    {
        mReinitMode = aMode;
    }

    //------------------------------------------------------------------------------

    ReinitMode
    Reinitialize_Performer::get_reinit_mode() const
    {
        return mReinitMode;
    }

    //------------------------------------------------------------------------------

    void
    Reinitialize_Performer::set_material_phase( moris_index aPhase )
    {
        mMaterialPhase = aPhase;
    }

    //------------------------------------------------------------------------------

    void
    Reinitialize_Performer::compute_sdf_from_interface(
            xtk::Cut_Integration_Mesh&                         aCutMesh,
            mtk::Mesh*                                         aInterpMesh,
            Vector< std::shared_ptr< gen::Geometry_Engine > >& aGENPerformer,
            Vector< std::shared_ptr< mtk::Mesh_Manager > >&    aMTKPerformer,
            std::shared_ptr< mtk::Field >&                     aSDFField )
    {
        Tracer tTracer( "WRK", "Reinitialize ADVs", "Compute SDF from Interface" );

        // Get interface facets from XTK
        Vector< moris_index > const& tInterfaceFacetIndices = aCutMesh.get_interface_facets();

        if ( tInterfaceFacetIndices.size() == 0 )
        {
            MORIS_LOG_WARNING( "No interface facets found in XTK mesh, skipping SDF reinitialization" );
            return;
        }

        // Get spatial dimension
        uint tSpatialDim = aInterpMesh->get_spatial_dim();
        uint tNodesPerFacet = ( tSpatialDim == 3 ) ? 3 : 2;  // Triangles in 3D, segments in 2D

        // Count total facet nodes (may have duplicates, handled by connectivity)
        uint tNumFacets = tInterfaceFacetIndices.size();

        // Build facet connectivity and collect unique vertices
        // For simplicity, we collect facet vertex coordinates directly
        // Note: This is a simplified approach - production code may need
        // to handle the facet-to-vertex mapping more carefully

        // Get facet connectivity from XTK's facet-based connectivity
        std::shared_ptr< xtk::Facet_Based_Connectivity > tFaceConn = aCutMesh.get_face_connectivity();

        MORIS_ERROR( tFaceConn != nullptr,
                "compute_sdf_from_interface - Face connectivity not available in XTK mesh" );

        // Collect all unique facet node coordinates
        std::map< moris_index, Matrix< DDRMat > > tUniqueVertexCoords;
        Matrix< IndexMat > tFacetConn( tNumFacets, tNodesPerFacet );

        for ( uint iFacet = 0; iFacet < tNumFacets; ++iFacet )
        {
            moris_index tFacetIndex = tInterfaceFacetIndices( iFacet );

            // Get vertices of this facet directly from mFacetVertices member
            Vector< moris::mtk::Vertex* > const& tFacetVertices = tFaceConn->mFacetVertices( tFacetIndex );

            for ( uint iNode = 0; iNode < tNodesPerFacet && iNode < tFacetVertices.size(); ++iNode )
            {
                if ( tFacetVertices( iNode ) == nullptr )
                {
                    continue;  // Skip null vertices
                }

                moris_index tVertexIndex = tFacetVertices( iNode )->get_index();
                tFacetConn( iFacet, iNode ) = tVertexIndex;

                // Store vertex coordinates if not already present
                if ( tUniqueVertexCoords.find( tVertexIndex ) == tUniqueVertexCoords.end() )
                {
                    tUniqueVertexCoords[ tVertexIndex ] = tFacetVertices( iNode )->get_coords();
                }
            }
        }

        // Build facet node coordinate matrix
        uint tNumUniqueVerts = tUniqueVertexCoords.size();
        Matrix< DDRMat > tFacetNodeCoords( tNumUniqueVerts, tSpatialDim );

        // Create mapping from original indices to compact indices
        std::map< moris_index, moris_index > tVertexIndexMap;
        moris_index tCompactIdx = 0;
        for ( auto const& [tOrigIdx, tCoords] : tUniqueVertexCoords )
        {
            tVertexIndexMap[ tOrigIdx ] = tCompactIdx;
            for ( uint iDim = 0; iDim < tSpatialDim; ++iDim )
            {
                tFacetNodeCoords( tCompactIdx, iDim ) = tCoords( iDim );
            }
            ++tCompactIdx;
        }

        // Update facet connectivity to use compact indices
        for ( uint iFacet = 0; iFacet < tNumFacets; ++iFacet )
        {
            for ( uint iNode = 0; iNode < tNodesPerFacet; ++iNode )
            {
                tFacetConn( iFacet, iNode ) = tVertexIndexMap[ tFacetConn( iFacet, iNode ) ];
            }
        }

        // Get interpolation mesh node coordinates
        uint tNumNodes = aInterpMesh->get_num_nodes();
        Matrix< DDRMat > tNodeCoords( tNumNodes, tSpatialDim );

        for ( uint iNode = 0; iNode < tNumNodes; ++iNode )
        {
            Matrix< DDRMat > tCoords = aInterpMesh->get_mtk_vertex( iNode ).get_coords();
            for ( uint iDim = 0; iDim < tSpatialDim; ++iDim )
            {
                tNodeCoords( iNode, iDim ) = tCoords( iDim );
            }
        }

        // Determine bulk phase per node using GEN
        Vector< moris_index > tNodeBulkPhase( tNumNodes );
        for ( uint iNode = 0; iNode < tNumNodes; ++iNode )
        {
            Matrix< DDRMat > tCoords = aInterpMesh->get_mtk_vertex( iNode ).get_coords();
            tNodeBulkPhase( iNode ) = aGENPerformer( 0 )->get_phase_index( iNode, tCoords );
        }

        // Compute SDF using SDF_From_Interface
        Matrix< DDRMat > tSDF;
        sdf::SDF_From_Interface::compute(
                tNumNodes,
                tNodeCoords,
                tFacetNodeCoords,
                tFacetConn,
                tNodeBulkPhase,
                mMaterialPhase,
                tSDF );

        // Create MTK Field with computed SDF values
        aSDFField = std::make_shared< mtk::Field_Discrete >(
                aMTKPerformer( 0 )->get_mesh_pair( 0 ), mAdofMeshIndex );
        aSDFField->set_label( mADVFiledName );
        aSDFField->unlock_field();
        aSDFField->set_values( tSDF );

        // Map nodal values to coefficients
        mtk::Mapper tMapper;
        aSDFField->unlock_field();
        tMapper.map_input_field_to_output_field_2( aSDFField.get() );

        // Compute coefficients from nodal values
        mCoefficients = aSDFField->get_coefficients();

        // Clip values to bounds
        this->impose_upper_lower_bound( aGENPerformer, aSDFField.get() );

        // Store field for GEN update
        Vector< std::shared_ptr< mtk::Field > > tGENFields;
        tGENFields.append( aGENPerformer( 0 )->get_mtk_fields() );

        // Find and replace the ADV field
        auto itr = std::find_if( tGENFields.begin(), tGENFields.end(),
                [ & ]( std::shared_ptr< mtk::Field > const& aField ) {
                    return aField->get_label() == mADVFiledName;
                } );

        if ( itr != tGENFields.end() )
        {
            moris_index tADVFieldIndex = std::distance( tGENFields.begin(), itr );
            tGENFields( tADVFieldIndex ) = aSDFField;
        }

        mMTKFields = tGENFields;

        MORIS_LOG_INFO( "SDF reinitialization complete: %d nodes, %d interface facets",
                tNumNodes, tNumFacets );
    }

}    // namespace moris::wrk
