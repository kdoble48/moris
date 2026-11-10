/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_HMR_T_Matrix.cpp
 *
 */

#include <catch.hpp>

#include "cl_HMR_T_Matrix.hpp"                //HMR/src
#include "cl_HMR_Background_Mesh_Base.hpp"    //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp"       //HMR/src
#include "cl_HMR_Factory.hpp"                 //HMR/src
#include "cl_HMR_Lagrange_Mesh_Base.hpp"      //HMR/src
#include "cl_HMR_Parameters.hpp"              //HMR/src

#include "cl_Communication_Manager.hpp"    //COM/src
#include "cl_Communication_Tools.hpp"      //COM/src

#include "paths.hpp"
#include "HDF5_Tools.hpp"
#include "fn_check_equal.hpp"

namespace moris::hmr
{
    // Test class for T-matrix
    template< uint N >
    class T_Matrix_Test : public T_Matrix< N >
    {
      public:
        // Constructor
        using T_Matrix< N >::T_Matrix;

        // Test evaluation
        void
        evaluate_shape_function_test(
                const Matrix< DDRMat >& aXi,
                Matrix< DDRMat >&       aN )
        {
            T_Matrix< N >::evaluate_shape_function( aXi, aN );
        }
    };

    // -----------------------------------------------------------------------------------------------------------------

    TEST_CASE( "HMR T-matrix", "[moris],[mesh],[hmr],[hmr_t_matrix]" )
    {
        // these tests are only performed in serial. They have nothing to do with parallel.
        if ( par_size() == 1 )
        {
            // Epsilon for shape function checks
            real tEpsilon = 1e-12;

            // Create parameters
            auto tParameters = new Parameters;

            // Open HDF5 file with read-only access
            std::string tFileName = get_base_moris_dir() + "/projects/HMR/test/data/T-Matrix_test_values.hdf5";
            hid_t       tFileID   = open_hdf5_file( tFileName, true, true );
            herr_t      tStatus;

            // Loop over number of dimensions
            for ( uint iNumberOfDimensions = 2; iNumberOfDimensions <= 3; iNumberOfDimensions++ )
            {
                // Create one element
                tParameters->set_number_of_elements_per_dimension( Vector< uint >( iNumberOfDimensions, 1 ) );

                // Domain from -1 to 1 in each dimension
                tParameters->set_domain_dimensions( Vector< real >( iNumberOfDimensions, 2 ) );
                tParameters->set_domain_offset( Vector< real >( iNumberOfDimensions, -1 ) );

                // Loop over orders
                for ( uint iOrder = 1; iOrder <= 3; iOrder++ )
                {
                    // set buffer size
                    tParameters->set_refinement_buffer( iOrder );
                    tParameters->set_staircase_buffer( iOrder );

                    // Loop over truncation
                    for ( bool iTruncation : { false, true } )
                    {
                        // create factory
                        Factory tFactory( tParameters );

                        // create background mesh object
                        Background_Mesh_Base* tBackgroundMesh = tFactory.create_background_mesh();

                        // create B-Spline Mesh
                        BSpline_Mesh_Base* tBSplineMesh = tFactory.create_dummy_bspline_mesh(
                                tBackgroundMesh,
                                0,
                                iOrder );

                        // create container of B-Spline meshes
                        Vector< BSpline_Mesh_Base* > tBSplineMeshes( 1, tBSplineMesh );

                        // create B-Spline Mesh
                        Lagrange_Mesh_Base* tLagrangeMesh = tFactory.create_lagrange_mesh(
                                tBackgroundMesh,
                                tBSplineMeshes,
                                0,
                                iOrder );

                        // create T-Matrix object
                        T_Matrix_Base* tTMatrix;
                        if ( iNumberOfDimensions == 2 )
                        {
                            tTMatrix = new T_Matrix_Test< 2 >( tLagrangeMesh, tBSplineMesh, iTruncation );
                        }
                        else
                        {
                            tTMatrix = new T_Matrix_Test< 3 >( tLagrangeMesh, tBSplineMesh, iTruncation );
                        }

                        // ask Lagrange mesh for number of nodes per element
                        luint tNumberOfNodes = tLagrangeMesh->get_number_of_nodes_on_proc();

                        // shape function vector
                        Matrix< DDRMat > tN( tNumberOfNodes, 1 );

                        // loop over all nodes
                        for ( uint iNodeIndex = 0; iNodeIndex < tNumberOfNodes; iNodeIndex++ )
                        {
                            // get pointer to node
                            Basis* tNode = tLagrangeMesh->get_basis_by_memory_index( iNodeIndex );

                            // get node coordinate
                            auto tXYZ = tNode->get_coords();

                            // Evaluate shape function
                            if ( iNumberOfDimensions == 2 )
                            {
                                static_cast< T_Matrix_Test< 2 >* >( tTMatrix )->evaluate_shape_function_test( tXYZ, tN );
                            }
                            else
                            {
                                static_cast< T_Matrix_Test< 3 >* >( tTMatrix )->evaluate_shape_function_test( tXYZ, tN );
                            }

                            // Check shape function vector
                            for ( uint iShapeFunctionIndex = 0; iShapeFunctionIndex < tNumberOfNodes; iShapeFunctionIndex++ )
                            {
                                if ( iShapeFunctionIndex == iNodeIndex )
                                {
                                    CHECK( std::abs( tN( iShapeFunctionIndex ) - 1 ) < tEpsilon );
                                }
                                else
                                {
                                    CHECK( std::abs( tN( iShapeFunctionIndex ) ) < tEpsilon );
                                }
                            }
                        }

                        // Perform refinement
                        tBackgroundMesh->set_activation_pattern( 0 );
                        for ( uint iLevel = 0; iLevel < 2; iLevel++ )
                        {
                            // Get single element on this level
                            Background_Element_Base* tElement = tBackgroundMesh->get_element( 0 );

                            // flag element for refinement
                            tElement->put_on_refinement_queue();

                            // Refine element
                            tBackgroundMesh->perform_refinement( 0 );
                        }

                        // Update B-spline mesh after refinement
                        tBSplineMesh->update_mesh();

                        // Calculate T-matrix
                        Matrix< DDRMat > tTMatrixCalculated;
                        Vector< Basis* >   tBasis;
                        tTMatrix->calculate_t_matrix( tBSplineMesh->get_element( 0 )->get_memory_index(),
                                tTMatrixCalculated,
                                tBasis );

                        // HDF5 label: Dimensions, order, truncation
                        std::string tLabel = "D" + std::to_string( iNumberOfDimensions )
                                           + "_O" + std::to_string( iOrder )
                                           + "_T" + std::to_string( iTruncation );

                        // Read expected T-matrix from HDF5 file and check with calculated values
                        Matrix< DDRMat > tTMatrixExpected;
                        load_matrix_from_hdf5_file( tFileID, tLabel, tTMatrixExpected, tStatus );
                        CHECK_EQUAL( tTMatrixCalculated, tTMatrixExpected, );

                        // tidy up memory
                        delete tTMatrix;
                        delete tBSplineMesh;
                        delete tLagrangeMesh;
                        delete tBackgroundMesh;
                    }
                }
            }

            // close file
            close_hdf5_file( tFileID );

            // delete settings object
            delete tParameters;
        }
    }

    // -----------------------------------------------------------------------------------------------------------------

    namespace
    {
        /**
         * Finds the active element of a mesh whose node bounding box contains a point (strictly inside)
         *
         * @param aMesh Lagrange mesh, its activation pattern is selected here
         * @param aPoint point inside the requested element
         * @return pointer to the element, nullptr if no element contains the point
         */
        Element*
        find_lagrange_element_containing_point(
                Lagrange_Mesh_Base*     aMesh,
                const Matrix< DDRMat >& aPoint )
        {
            aMesh->select_activation_pattern();

            uint tNumberOfDimensions   = aPoint.numel();
            uint tNumberOfNodesPerElem = aMesh->get_number_of_bases_per_element();

            for ( luint iElement = 0; iElement < aMesh->get_number_of_elements(); iElement++ )
            {
                Element* tElement = aMesh->get_element( iElement );

                Matrix< DDRMat > tMin( tNumberOfDimensions, 1, MORIS_REAL_MAX );
                Matrix< DDRMat > tMax( tNumberOfDimensions, 1, -MORIS_REAL_MAX );

                for ( uint iNode = 0; iNode < tNumberOfNodesPerElem; iNode++ )
                {
                    Matrix< DDRMat > tCoords = tElement->get_basis( iNode )->get_coords();

                    for ( uint iDim = 0; iDim < tNumberOfDimensions; iDim++ )
                    {
                        tMin( iDim ) = std::min( tMin( iDim ), tCoords( iDim ) );
                        tMax( iDim ) = std::max( tMax( iDim ), tCoords( iDim ) );
                    }
                }

                bool tInside = true;
                for ( uint iDim = 0; iDim < tNumberOfDimensions; iDim++ )
                {
                    tInside = tInside && aPoint( iDim ) > tMin( iDim ) && aPoint( iDim ) < tMax( iDim );
                }

                if ( tInside )
                {
                    return tElement;
                }
            }

            return nullptr;
        }

        /**
         * Value of the linear tensor-product B-spline polynomial of an element at a point; the polynomial is extended
         * beyond the element. The element spans [aLower, aLower + 1] in every direction.
         *
         * @param aCorner local corner of the basis: 0 (lower) or 1 (upper) per direction
         * @param aLower lower corner of the B-spline element
         * @param aPoint evaluation point
         * @return value of the extended basis function
         */
        real
        extended_linear_bspline_value(
                const Vector< luint >&  aCorner,
                const Matrix< DDRMat >& aLower,
                const Matrix< DDRMat >& aPoint )
        {
            real tValue = 1.0;
            for ( uint iDim = 0; iDim < aCorner.size(); iDim++ )
            {
                real tLocal = aPoint( iDim ) - aLower( iDim );
                tValue *= aCorner( iDim ) == 0 ? 1.0 - tLocal : tLocal;
            }
            return tValue;
        }

        /**
         * Creates a T-matrix object for a Lagrange and a B-spline mesh of a given dimension
         */
        T_Matrix_Base*
        create_test_t_matrix(
                uint                aNumberOfDimensions,
                Lagrange_Mesh_Base* aLagrangeMesh,
                BSpline_Mesh_Base*  aBSplineMesh )
        {
            if ( aNumberOfDimensions == 2 )
            {
                return new T_Matrix< 2 >( aLagrangeMesh, aBSplineMesh, true );
            }
            return new T_Matrix< 3 >( aLagrangeMesh, aBSplineMesh, true );
        }
    }    // namespace

    // -----------------------------------------------------------------------------------------------------------------

    // MB-104: the extended T-matrix (cell agglomeration) evaluates the root B-spline element's linear basis at the
    // Lagrange nodes of an element some elements away. The node coordinates in the root element's parameter space were
    // wrong for Lagrange order > 1 (element shift scaled by 1/P) and when the Lagrange element is finer than the
    // B-spline element (ijk of different levels mixed). The weights are compared with the extended linear polynomials
    // evaluated at the physical node coordinates.
    TEST_CASE( "HMR extended T-matrix", "[moris],[mesh],[hmr],[hmr_t_matrix],[hmr_extended_t_matrix]" )
    {
        if ( par_size() == 1 )
        {
            real tEpsilon = 1e-12;

            for ( uint iNumberOfDimensions = 2; iNumberOfDimensions <= 3; iNumberOfDimensions++ )
            {
                for ( uint iLagrangeOrder = 1; iLagrangeOrder <= 3; iLagrangeOrder++ )
                {
                    // the Lagrange mesh is on the B-spline pattern, or on a pattern refined once more
                    for ( bool iLagrangeRefined : { false, true } )
                    {
                        // 4 x 4 (x 4) elements of size 1, domain [0,4]^N
                        Parameters* tParameters = new Parameters;
                        tParameters->set_number_of_elements_per_dimension( Vector< uint >( iNumberOfDimensions, 4 ) );
                        tParameters->set_domain_dimensions( Vector< real >( iNumberOfDimensions, 4.0 ) );
                        tParameters->set_domain_offset( Vector< real >( iNumberOfDimensions, 0.0 ) );
                        tParameters->set_refinement_buffer( 1 );
                        tParameters->set_staircase_buffer( 1 );

                        Factory               tFactory( tParameters );
                        Background_Mesh_Base* tBackgroundMesh = tFactory.create_background_mesh();

                        uint tBSplinePattern  = 0;
                        uint tLagrangePattern = iLagrangeRefined ? 1 : 0;

                        if ( iLagrangeRefined )
                        {
                            // refine the Lagrange pattern once everywhere
                            tBackgroundMesh->reset_pattern( tLagrangePattern );
                            tBackgroundMesh->set_activation_pattern( tLagrangePattern );
                            for ( luint iElement = 0; iElement < tBackgroundMesh->get_number_of_active_elements_on_proc(); iElement++ )
                            {
                                tBackgroundMesh->get_element( iElement )->put_on_refinement_queue();
                            }
                            tBackgroundMesh->perform_refinement( tLagrangePattern );
                        }

                        BSpline_Mesh_Base* tBSplineMesh = tFactory.create_dummy_bspline_mesh( tBackgroundMesh, tBSplinePattern, 1 );

                        Vector< BSpline_Mesh_Base* > tBSplineMeshes( 1, tBSplineMesh );
                        Lagrange_Mesh_Base*          tLagrangeMesh = tFactory.create_lagrange_mesh(
                                tBackgroundMesh,
                                tBSplineMeshes,
                                tLagrangePattern,
                                iLagrangeOrder );

                        T_Matrix_Base* tTMatrix = create_test_t_matrix( iNumberOfDimensions, tLagrangeMesh, tBSplineMesh );

                        // root B-spline element: [1,2]^N, found through the level-0 ancestor of a Lagrange element in it
                        Matrix< DDRMat > tBSplineLower( iNumberOfDimensions, 1, 1.0 );
                        Matrix< DDRMat > tPointInRoot( iNumberOfDimensions, 1, 1.4 );
                        Element*         tLagrangeInRoot = find_lagrange_element_containing_point( tLagrangeMesh, tPointInRoot );
                        REQUIRE( tLagrangeInRoot != nullptr );

                        Background_Element_Base* tRootBackgroundElement = tLagrangeInRoot->get_background_element();
                        while ( tRootBackgroundElement->get_level() > 0 )
                        {
                            tRootBackgroundElement = tRootBackgroundElement->get_parent();
                        }
                        Element* tBSplineElement = tBSplineMesh->get_element_by_memory_index( tRootBackgroundElement->get_memory_index() );

                        // corner of each basis of the root element, from the basis ijk relative to the smallest one
                        uint            tNumberOfBases = tBSplineMesh->get_number_of_bases_per_element();
                        Vector< luint > tMinIJK( iNumberOfDimensions, MORIS_LUINT_MAX );
                        for ( uint iBasis = 0; iBasis < tNumberOfBases; iBasis++ )
                        {
                            const luint* tIJK = tBSplineElement->get_basis( iBasis )->get_ijk();
                            for ( uint iDim = 0; iDim < iNumberOfDimensions; iDim++ )
                            {
                                tMinIJK( iDim ) = std::min( tMinIJK( iDim ), tIJK[ iDim ] );
                            }
                        }

                        // Lagrange elements to extend into: one and two elements away (and diagonal)
                        Vector< Vector< real > > tTargetCoordinates;
                        if ( iNumberOfDimensions == 2 )
                        {
                            tTargetCoordinates = { { 2.3, 1.4 }, { 3.8, 2.3 }, { 0.3, 0.8 } };
                        }
                        else
                        {
                            tTargetCoordinates = { { 2.3, 1.4, 1.6 }, { 3.8, 2.3, 0.2 } };
                        }

                        Vector< Matrix< DDRMat > > tTargetPoints;
                        for ( const Vector< real >& tCoordinates : tTargetCoordinates )
                        {
                            Matrix< DDRMat > tPoint( iNumberOfDimensions, 1 );
                            for ( uint iDim = 0; iDim < iNumberOfDimensions; iDim++ )
                            {
                                tPoint( iDim ) = tCoordinates( iDim );
                            }
                            tTargetPoints.push_back( tPoint );
                        }

                        for ( const Matrix< DDRMat >& tTargetPoint : tTargetPoints )
                        {
                            Element* tLagrangeElement = find_lagrange_element_containing_point( tLagrangeMesh, tTargetPoint );
                            REQUIRE( tLagrangeElement != nullptr );
                            REQUIRE( tLagrangeElement->get_level() == ( iLagrangeRefined ? 1u : 0u ) );

                            Vector< Vector< mtk::Vertex* > > tBSplineBasis;
                            Vector< Matrix< DDRMat > >       tWeights;
                            tTMatrix->evaluate_extended_t_matrix( tBSplineElement, tLagrangeElement, tBSplineBasis, tWeights );

                            uint tNumberOfNodes = tLagrangeMesh->get_number_of_bases_per_element();
                            REQUIRE( tBSplineBasis.size() == tNumberOfNodes );

                            for ( uint iNode = 0; iNode < tNumberOfNodes; iNode++ )
                            {
                                CAPTURE( iNumberOfDimensions, iLagrangeOrder, iLagrangeRefined, iNode );

                                Matrix< DDRMat > tNodeCoords = tLagrangeElement->get_basis( iNode )->get_coords();

                                real tSumOfWeights     = 0.0;
                                uint tNumberOfExpected = 0;

                                // every basis of the root element with a nonzero value at the node must be returned with that value
                                for ( uint iBasis = 0; iBasis < tNumberOfBases; iBasis++ )
                                {
                                    Basis*          tBasis = tBSplineElement->get_basis( iBasis );
                                    const luint*    tIJK   = tBasis->get_ijk();
                                    Vector< luint > tCorner( iNumberOfDimensions );
                                    for ( uint iDim = 0; iDim < iNumberOfDimensions; iDim++ )
                                    {
                                        tCorner( iDim ) = tIJK[ iDim ] - tMinIJK( iDim );
                                    }

                                    real tExpected = extended_linear_bspline_value( tCorner, tBSplineLower, tNodeCoords );

                                    real tComputed = 0.0;
                                    for ( uint iEntry = 0; iEntry < tBSplineBasis( iNode ).size(); iEntry++ )
                                    {
                                        if ( tBSplineBasis( iNode )( iEntry ) == static_cast< mtk::Vertex* >( tBasis ) )
                                        {
                                            tComputed = tWeights( iNode )( iEntry );
                                        }
                                    }

                                    CHECK( std::abs( tComputed - tExpected ) < tEpsilon );

                                    if ( std::abs( tExpected ) > tEpsilon )
                                    {
                                        tNumberOfExpected++;
                                    }
                                }

                                // no other basis is returned
                                CHECK( tBSplineBasis( iNode ).size() == tNumberOfExpected );

                                // partition of unity
                                for ( uint iEntry = 0; iEntry < tWeights( iNode ).numel(); iEntry++ )
                                {
                                    tSumOfWeights += tWeights( iNode )( iEntry );
                                }
                                CHECK( std::abs( tSumOfWeights - 1.0 ) < tEpsilon );
                            }
                        }

                        delete tTMatrix;
                        delete tLagrangeMesh;
                        delete tBSplineMesh;
                        delete tBackgroundMesh;
                        delete tParameters;
                    }
                }
            }
        }
    }
}    // namespace moris::hmr
