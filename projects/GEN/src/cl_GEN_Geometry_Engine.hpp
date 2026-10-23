/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Geometry_Engine.hpp
 *
 */

#pragma once

// GEN
#include "st_GEN_Geometry_Engine_Parameters.hpp"
#include "cl_GEN_Phase_Table.hpp"
#include "cl_GEN_ADV_Manager.hpp"
#include "cl_GEN_Node_Manager.hpp"
#include "cl_GEN_PDV_Host_Manager.hpp"

// MTK
#include "cl_MTK_Enums.hpp"
#include <unordered_map>
#include <map>
#include <array>

namespace moris::gen
{
    class Geometry;
    class Property;
}

namespace moris::MSI
{
    // Forward declare design variable interface
    class Design_Variable_Interface;
}

namespace moris::sol
{
    class Dist_Vector;
}

namespace moris::gen
{
    class Geometry_Engine
    {

      protected:

        // Protected for access in test version
        Node_Manager                        mNodeManager;
        Vector< std::shared_ptr< Geometry > > mGeometries;
        Vector< std::shared_ptr< Property > > mProperties;

      private:

        // Phase Table
        Phase_Table mPhaseTable;

        // Spatial dimensions
        uint mNumSpatialDimensions;

        // ADVs
        ADV_Manager       mADVManager;
        Vector< real >    mInitialPrimitiveADVs;
        Vector< sint >    mFullADVIds;
        Matrix< IdMat >   mFullijklIDs;
        sol::Dist_Vector* mOwnedADVs     = nullptr;
        sol::Dist_Vector* mPrimitiveADVs = nullptr;

        // Kreissl-Maute ADV filter (adv_filter_radius; Kreissl and Maute 2012, SMO 46:311-326, Eq. 12, with r = r_s h): the
        // filtered rows of F (coefficients = F * ADVs) on the full ADV vector; positions without a row are identity rows.
        // Built on processor 0 by build_adv_filter().
        struct ADV_Filter_Row
        {
            uint           mPosition;    // position of the filtered coefficient in the full ADV vector
            Vector< uint > mColumns;     // positions of the ADVs it averages
            Vector< real > mWeights;     // normalized weights max( 0, r - |x_i - x_j| ) / sum
        };

        // Neighbor search of the ADV filter: coefficients are sorted into buckets of edge length r (integer bucket
        // coordinates), and global ADV IDs are mapped to positions in the full ADV vector
        using Filter_Bucket_Key = std::array< sint, 3 >;
        using Filter_Buckets    = std::map< Filter_Bucket_Key, Vector< uint > >;
        using ADV_Position_Map  = std::unordered_map< sint, uint >;
        Vector< ADV_Filter_Row > mADVFilterRows;
        Vector< real >           mUnfilteredADVs;    // the optimizer's ADVs last passed to set_advs()

        // IQIs
        Vector< std::string > mRequestedIQIs;

        size_t      mActiveGeometryIndex = 0;
        std::string mGeometryFieldFile;
        std::string mOutputMeshFile;
        bool        mShapeSensitivities = false;
        real        mTimeOffset;

        // PDVs
        PDV_Host_Manager                     mPDVHostManager;
        Intersection_Node* mQueuedIntersectionNode = nullptr;

        // diagnostic information
        bool        mDiagnostics    = false;
        std::string mDiagnosticPath;
        std::string mDiagnosticId;

      public:

        /**
         * Constructor using cell of cell of parameter lists
         *
         * @param aParameterLists GEN parameter lists (see fn_PRM_GEN_Parameters.hpp)
         * @param aLibrary Library used for pulling user-defined functions
         * @param aMesh Mesh for discrete or mesh based geomtries
         */
        explicit Geometry_Engine(
                Module_Parameter_Lists        aParameterLists,
                const std::shared_ptr< Library_IO >& aLibrary = nullptr,
                mtk::Mesh*                           aMesh    = nullptr );

        /**
         * Constructor
         *
         * @param aMesh Mesh for getting B-spline information
         * @param aParameters Optional geometry engine parameters
         */
        explicit Geometry_Engine(
                mtk::Interpolation_Mesh*   aMesh,
                Geometry_Engine_Parameters aParameters = {} );

        /**
         * Destructor
         */
        ~Geometry_Engine();

        /**
         * Sets new ADVs for the geometry engine.
         *
         * @param aNewADVs vector of new advs to use
         */
        void set_advs( const Vector< real >& aNewADVs );

        /**
         * Gets the advs from the geometry engine
         *
         * @return vector of advs
         */
        Vector< real >& get_advs();

        /**
         * Whether any design filters its B-spline coefficients (adv_filter_radius > 0). With a filter, set_advs() and
         * get_advs() exchange the optimizer's variables psi, the fields use F * psi, and get_dcriteria_dadv() returns
         * F^T times the coefficient sensitivities.
         *
         * @return true if an ADV filter is active
         */
        bool has_adv_filter() const;

        /**
         * Applies the ADV filter to a full ADV vector: F * aADVs (identity outside the filtered designs).
         *
         * @param aADVs Full ADV vector (processor 0)
         * @return Filtered vector
         */
        Vector< real > apply_adv_filter( const Vector< real >& aADVs ) const;

        /**
         * Applies the transposed ADV filter to a full vector of sensitivities: F^T * aValues.
         *
         * @param aValues Derivatives with respect to the filtered coefficients (processor 0)
         * @return Derivatives with respect to the ADVs
         */
        Vector< real > apply_adv_filter_transpose( const Vector< real >& aValues ) const;

        /**
         * Gets the global ADV IDs owned by a named design's field, in coefficient order. Lets a
         * caller (e.g. the SDF reinitialization workflow) write new coefficients into exactly that
         * design's slice of the global ADV vector, leaving every other design's ADVs and all PDVs
         * untouched -- avoiding the full distribute_advs() rebuild that desyncs FEM-coupled PDVs.
         *
         * @param aDesignName Name of the geometry/design (matches the "name" parameter)
         * @return Determining ADV IDs of that design's field coefficients; empty if not found
         */
        Vector< sint > get_design_determining_adv_ids( const std::string& aDesignName );

        /**
         * Positions in the optimizer's full ADV vector (the vector passed to set_advs()) of a named
         * design's field coefficients, in coefficient order. The determining ADV IDs above are GLOBAL
         * IDs of the distributed ADV map, not positions: set_advs() places aNewADVs( k ) at global ID
         * mFullADVIds( k ). Writing a coefficient at its ID instead of its position shifts the design
         * whenever IDs and positions differ (MB-017: by one index in a serial MBB run).
         *
         * @param aDesignName Name of the geometry/design (matches the "name" parameter)
         * @return Position of each coefficient in the full ADV vector, -1 if its ID is not in it; empty if not found
         */
        Vector< sint > get_design_adv_positions( const std::string& aDesignName );

        /**
         * Gets the names of all designs (geometries and properties) discretized with discretization_mode = "additive".
         *
         * @return Names of additively discretized designs (empty if none)
         */
        Vector< std::string > get_additive_discretized_design_names() const;

        /**
         * Gets whether any design (geometry or property) is discretized with discretization_mode = "additive".
         * Such designs are incompatible with HMR remeshing and level-set reinitialization.
         *
         * @return true if at least one design is discretized additively
         */
        bool has_additive_discretized_designs() const;

        /**
         * Get vector with ijkl IDs. All Ids are on proc 0, all others return empty vec
         *
         * @return vector of aijkl Ids
         */
        Matrix< IdMat >& get_IjklIDs();

        /**
         * Gets the lower bounds from the geometry engine
         *
         * @return vector of lower bounds
         */
        Vector< real >& get_lower_bounds();

        /**
         * Gets the upper bounds from the geometry engine
         *
         * @return vector of upper bounds
         */
        Vector< real >& get_upper_bounds();

        /**
         * Lets MDL know about the stored requested IQIs through the PDV host manager. This has to be done after
         * the model is set, that's why it's a separate call that needs to be performed at the right time.
         */
        void communicate_requested_IQIs();

        /**
         * Import a phase function pointer
         * used in plato workflow where platoanalyze does the forward analysis and no library is present
         */
        void
        set_phase_function(
                PHASE_FUNCTION      aPhaseFunction,
                uint                aNumPhases,
                const Vector< std::string >& aPhaseNames = {} );

        /**
         * Import dcriteria/dx from file
         * used in plato workflow where platoanalyze does the forward analysis
         */
        void
        set_dQIdp(
                const Vector< Matrix< DDRMat >* >& adQIdp,
                Matrix< DDSMat >*                aMap );

        /**
         * Gets the sensitivities of the criteria with respect to the advs
         *
         * @return Vector of sensitivities
         */
        Matrix< DDRMat > get_dcriteria_dadv();

        /**
         * Gets the design variable interface from the geometry engine
         *
         * @return member pdv host manager pointer
         */
        MSI::Design_Variable_Interface* get_design_variable_interface();

        /**
         * Returns if the element consisting of provided node indices and coordinates is intersected by the given geometry.
         *
         * @param aGeometryIndex Geometry index
         * @param aNodeIndices Node indices
         * @param aNodeCoordinates Node coordinates
         * @return If the element is intersected
         */
        bool is_intersected(
                uint                      aGeometryIndex,
                const Matrix< IndexMat >& aNodeIndices,
                const Matrix< DDRMat >&   aNodeCoordinates );

        /**
         * Returns if the given geometry can change geometric region inside a background cell whose corner nodes all lie
         * in one region (an analytic level set evaluated exactly at derived nodes, MB-031). If so, cut-cell detection
         * must also test the points at which the decomposition will create vertices.
         *
         * @param aGeometryIndex Geometry index
         * @return If interior points need to be tested
         */
        bool can_change_region_inside_cell( uint aGeometryIndex );

        /**
         * Returns if the element is intersected by the given geometry, testing its nodes and additionally the given
         * interior points (the vertices the decomposition will create inside the element). The interior points are only
         * evaluated for geometries that can change region inside a cell (see can_change_region_inside_cell()).
         *
         * @param aGeometryIndex Geometry index
         * @param aNodeIndices Node indices
         * @param aNodeCoordinates Node coordinates, indexed by node index
         * @param aInteriorPoints Global coordinates of the interior points (row vectors)
         * @return If the element is intersected
         */
        bool is_intersected(
                uint                                                         aGeometryIndex,
                const Matrix< IndexMat >&                                    aNodeIndices,
                Vector< std::shared_ptr< moris::Matrix< moris::DDRMat > > >* aNodeCoordinates,
                const Vector< Matrix< DDRMat > >&                            aInteriorPoints );

        /**
         * Returns if the element consisting of provided node indices and coordinates is intersected by the given geometry.
         *
         * @param aGeometryIndex Geometry index
         * @param aNodeIndices Node indices
         * @param aNodeCoordinates Node coordinates
         * @return If the element is intersected
         */
        bool is_intersected(
                uint                                         aGeometryIndex,
                const Matrix< IndexMat >&                    aNodeIndices,
                Vector< std::shared_ptr< Matrix< DDRMat > > >* aNodeCoordinates );

        /**
         * Determines if the given edge is intersected, and queues an intersection node if it is. If an intersection
         * node has been queued, questions can be asked about the queued node:
         *
         * @param aEdgeFirstNodeIndex First node index on the intersection edge
         * @param aEdgeSecondNodeIndex Second node index on the intersection edge
         * @param aEdgeFirstNodeParametricCoordinates Local coordinates of the first node inside the background element
         * @param aEdgeSecondNodeParametricCoordinates Local coordinates of the second node inside the background element
         * @param aBackgroundElementNodeIndices Node indices of the background element
         * @param aBackgroundGeometryType Geometry type of the background element
         * @param aBackgroundInterpolationOrder Interpolation order of the background element
         *
         * @return If the edge is intersected and a node has been queued
         */
        bool queue_intersection(
                uint                     aEdgeFirstNodeIndex,
                uint                     aEdgeSecondNodeIndex,
                const Matrix< DDRMat >&  aEdgeFirstNodeParametricCoordinates,
                const Matrix< DDRMat >&  aEdgeSecondNodeParametricCoordinates,
                const Matrix< DDUMat >&  aBackgroundElementNodeIndices,
                mtk::Geometry_Type       aBackgroundGeometryType,
                mtk::Interpolation_Order aBackgroundInterpolationOrder );

        /**
         * Returns if the queued intersection has the first parent node on the active geometry interface.
         *
         * @return If the first parent node is on the interface
         */
        bool queued_intersection_first_parent_on_interface();

        /**
         * Returns if the queued intersection has the second parent node on the active geometry interface.
         *
         * @return If the second parent node is on the interface
         */
        bool queued_intersection_second_parent_on_interface();

        /**
         * Gets the local coordinate of the queued intersection node.
         *
         * @return Intersection node local coordinate (between -1 and 1)
         */
        real get_queued_intersection_local_coordinate();

        /**
         * Gets the global coordinates of the queued intersection node.
         *
         * @return Intersection node global coordinates
         */
        Matrix< DDRMat > get_queued_intersection_global_coordinates();

        /**
         * Admit the queued intersection as a unique, permanent node for sensitivity calculations.
         */
        void admit_queued_intersection();

        /**
         * Update the queued intersection node with its node ID and node owner.
         *
         * @param aNodeIndex Node Index
         * @param aNodeId Node ID
         * @param aNodeOwner Node owner
         */
        void update_intersection_node(
                uint        aNodeIndex,
                moris_id    aNodeId,
                moris_index aNodeOwner );

        /**
         * Creates and registers new derived nodes based on the given information.
         *
         * @param aVertexIndices Indices of the parent cell
         * @param aParametricCoordinates Parametric coordinates of each new derived node to create
         * @param aBackgroundGeometryType Geometry type of the background element
         * @param aBackgroundInterpolationOrder Interpolation order of the background element
         */
        void create_new_derived_nodes(
                const Vector< Matrix< IndexMat > >& aVertexIndices,
                const Vector< Matrix< DDRMat > >&   aParametricCoordinates,
                mtk::Geometry_Type                aBackgroundGeometryType,
                mtk::Interpolation_Order          aBackgroundInterpolationOrder );

        /**
         * Overloaded version of creating derived nodes. Calls other version internally.
         *
         * @param aNewNodeParentCell MTK cells
         * @param aParametricCoordinates Parametric coordinates for creating the derived node
         */
        void create_new_derived_nodes(
                Vector< mtk::Cell* >&                                aNewNodeParentCell,
                const Vector< std::shared_ptr< Matrix< DDRMat > > >& aParametricCoordinates );

        /**
         * Get the total number of phases in the phase table
         */
        size_t get_num_phases();

        /**
         * Gets the phase of a given node
         *
         * @param aNodeIndex Node index
         * @param aNodeCoordinates Node coordinates
         * @return Phase index
         */
        size_t get_phase_index(
                moris_index             aNodeIndex,
                const Matrix< DDRMat >& aNodeCoordinates );

        /**
         * Gets the geometric region of a node with respect to a given geometry.
         *
         * @param aGeometryIndex Geometry index
         * @param aNodeIndex Node index
         * @param aCoordinates Node coordinates
         * @return Geometric region
         */
        Geometric_Region get_geometric_region(
                uint                    aGeometryIndex,
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aNodeCoordinates );

        /**
         * Determines the phase of an element based on the geometric regions at each vertex.
         *
         * @param aCell Cell with vertices for geometric region evaluation
         * @return element phase index
         */
        moris_index get_element_phase_index( const mtk::Cell& aCell );

        /**
         * @brief Returns the number of geometries
         */
        size_t get_number_of_geometries();

        //-------------------------------------------------------------------------------

        /**
         * @brief Returns the number of phases
         */
        size_t get_num_bulk_phase();

        //-------------------------------------------------------------------------------

        /**
         * @brief Returns the active geometry index
         */
        size_t get_active_geometry_index();

        //-------------------------------------------------------------------------------

        /**
         * @brief Advance the active geometry index
         */
        void advance_geometry_index();

        //-------------------------------------------------------------------------------

        /**
         * Return the number of fields that can be used for refinement
         *
         * @return Number of fields for refinement
         */
        uint get_num_refinement_fields();

        //-------------------------------------------------------------------------------

        const Vector< uint >& get_num_refinements( uint aFieldIndex );

        //-------------------------------------------------------------------------------

        const Vector< uint >& get_refinement_mesh_indices( uint aFieldIndex );

        /**
         * Gets all of the MTK fields that the geometry engine is using.
         *
         * @return MTK fields
         */
        Vector< std::shared_ptr< mtk::Field > > get_mtk_fields();

        /**
         * Gets the index of an HMR user-defined refinement function for the given field index
         *
         * @param aFieldIndex Index of the field
         * @param aRefinementIndex The current refinement step being performed
         * @return User-defined function index, or -1 to use default refinement
         */
        sint get_refinement_function_index(
                uint aFieldIndex,
                uint aRefinementIndex );

        /**
         * Discretize GEN fields on the given mesh and distribute parallel ADVs based on these fields.
         *
         * @param aMeshPair Mesh for discretizing fields
         */
        void distribute_advs(
                mtk::Mesh_Pair                                 aMeshPair,
                const Vector< std::shared_ptr< mtk::Field > >& aFields        = {},
                mtk::EntityRank                                aADVEntityRank = mtk::EntityRank::BSPLINE );

        /**
         * Resets the information that the geometry engine stores about a mesh.
         *
         * @param aMesh Mesh for computing level set data
         */
        void reset_mesh_information( mtk::Interpolation_Mesh* aMesh );

        /**
         * Outputs geometry and property fields on the given mesh, and writes level set fields to a text file.
         * Uses output locations and file names stored from a parameter list or previous call to an output function.
         *
         * @param aMesh Mesh to evaluate fields on
         */
        void output_fields( mtk::Mesh* aMesh );

        /**
         * Creates geometry and property fields on the given mesh, and writes the mesh to an exodus file.
         *
         * @param aMesh Mesh to evaluate fields on
         * @param aExodusFileName Name of an Exodus file to write the mesh to
         */
        void output_fields_on_mesh(
                mtk::Mesh*  aMesh,
                std::string aExodusFileName );

        /**
         * Writes all geometry fields to separate text files with the given base file name (suffix appended on to
         * identify each individual geometry by index.
         *
         * @param aMesh Mesh to evaluate fields on
         * @param aFileName Base name of text files to write the geometry field data to
         */
        void write_geometry_fields(
                mtk::Mesh*         aMesh,
                const std::string& aBaseFileName );

        /**
         * Assign PDV hosts based on properties constructed through parameter lists
         *
         * @param aMeshManager Mesh manager
         */
        void create_pdvs( const mtk::Mesh_Pair& aMeshPair );

        //-------------------------------------------------------------------------------

        void
        print_gen_vertices(
                const std::string& aFile,
                mtk::Mesh*         aMesh );

        //-------------------------------------------------------------------------------

        void
        setup_diagnostics(
                bool                aDiagnostics,
                std::string const & aDiagnosticPath,
                std::string const & aDiagnosticLabel );

        std::string
        get_diagnostic_file_name( std::string const & aLabel ) const;

        //-------------------------------------------------------------------------------

      private:
        /**
         * Builds the ADV filter rows (mADVFilterRows) for every geometry with adv_filter_radius > 0 from the Greville
         * coordinates of its B-spline coefficients. Called at the end of distribute_advs(); serial only.
         */
        void build_adv_filter();

        /**
         * Appends the filter rows of one B-spline field to mADVFilterRows: for each coefficient i with an ADV, the
         * normalized weights max( 0, r - |x_i - x_j| ) of all coefficients j within the radius (Kreissl and Maute 2012,
         * Eq. 12).
         *
         * @param aCoordinates Greville coordinates of the field's coefficients, one row per coefficient
         * @param aADVIds Global ADV ID of each coefficient
         * @param aRadius Filter radius in physical units
         * @param aPositionOfId Position of each global ADV ID in the full ADV vector
         * @return Number of nonzero filter weights added
         */
        uint add_adv_filter_rows(
                const Matrix< DDRMat >& aCoordinates,
                const Vector< sint >&   aADVIds,
                real                    aRadius,
                const ADV_Position_Map& aPositionOfId );

        static void communicate_missing_owned_coefficients(
                mtk::Mesh_Pair&  aMeshPair,
                Matrix< IdMat >& aAllCoefIds,
                Matrix< IdMat >& aAllCoefOwners,
                Matrix< IdMat >& aAllCoefijklIds,
                Vector< uint >&    aNumCoeff,
                uint             aFieldIndex,
                uint             aDiscretizationMeshIndex,
                mtk::MeshType    aMeshType );

        /**
         * Create PDV_Type hosts with the specified PDV_Type types on the interpolation mesh
         *
         * @param aPDVTypes PDV_Type types; set->group->individual
         * @param aMeshIndex Interpolation mesh index
         */
        void create_interpolation_pdvs(
                mtk::Interpolation_Mesh*         aInterpolationMesh,
                mtk::Integration_Mesh*           aIntegrationMesh,
                Vector< Vector< Vector< PDV_Type > > > aPDVTypes );

        /**
         * Create PDV_Type hosts with PDVs for each of the spatial dimensions on the integration mesh
         *
         * @param aMeshIndex Integration mesh index
         */
        void set_integration_pdv_types( mtk::Integration_Mesh* aIntegrationMesh );

        /**
         * Initialize the PDV type list.
         */
        void initialize_pdv_type_list();

        /**
         * Decides how to construct the phase table based on the given arguments.
         *
         * @param aParameterLists GEN parameter lists (see fn_PRM_GEN_Parameters.hpp)
         * @return Phase table
         */
        static Phase_Table create_phase_table(
                const Module_Parameter_Lists& aParameterLists,
                const std::shared_ptr< Library_IO >& aLibrary );

        /**
         * Decides how to construct the phase table based on the given parameter lists
         *
         * @param aNumGeometries Number of geometries
         * @param aBulkPhases Phase table vector
         * @param aPhaseFunction Custom phase function
         * @return Phase table
         */
        static Phase_Table create_phase_table(
                uint                    aNumGeometries,
                const Matrix< DDUMat >& aBulkPhases,
                PHASE_FUNCTION          aPhaseFunction = nullptr,
                uint                    aNumPhases     = 1 );

    };
}
