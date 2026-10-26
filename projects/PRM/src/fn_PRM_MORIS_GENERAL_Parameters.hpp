/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_PRM_MORIS_GENERAL_Parameters.hpp
 *
 */

#pragma once

#include "cl_Parameter_List.hpp"

namespace moris::prm
{

    //------------------------------------------------------------------------------

    // creates a parameter list with default inputs
    inline Parameter_List
    create_moris_general_parameter_list()
    {
        Parameter_List tParameterList( "General" );

        return tParameterList;
    }

    //------------------------------------------------------------------------------

    inline void
    create_refinement_parameterlist( Parameter_List& aParameterlist )
    {
        aParameterlist.insert( "field_names", "" );
        aParameterlist.insert( "levels_of_refinement", "" );
        aParameterlist.insert( "refinement_pattern", "" );

        aParameterlist.insert( "refinement_function_name", "" );

        aParameterlist.insert( "remeshing_copy_old_pattern_to_pattern", "" );
    }

    //------------------------------------------------------------------------------

    inline void
    create_remeshing_parameterlist( Parameter_List& aParameterlist )
    {
        // Remeshing mode. Options are "ab_initio", "former"
        aParameterlist.insert( "mode", "" );
        aParameterlist.insert( "remeshing_refinement_pattern", "" );
        aParameterlist.insert( "refinement_function_name", "" );

        aParameterlist.insert( "remeshing_frequency", MORIS_SINT_MAX );

        // mode "ab_initio"
        aParameterlist.insert( "remeshing_field_names", "" );
        aParameterlist.insert( "remeshing_levels_of_refinement", "" );
        aParameterlist.insert( "remeshing_refinement_pattern", "" );

        // mode "based_on_previous"
        aParameterlist.insert( "remeshing_maximum_refinement_level", "" );
        aParameterlist.insert( "remeshing_minimum_refinement_level", "" );

        aParameterlist.insert( "remeshing_copy_old_pattern_to_pattern", "" );

        // minimum refinement level per pattern and level
        // input: pattern, level, iter, level, iter,level, iter, ... ; pattern, level, iter, level, iter,level, iter, ...
        aParameterlist.insert( "minimum_refinement_level", "" );

        aParameterlist.insert( "output_meshes", false );
    }

    //------------------------------------------------------------------------------

    inline void
    create_mapping_parameterlist( Parameter_List& aParameterlist )
    {
        aParameterlist.insert( "adv_field", "" );
        aParameterlist.insert( "dof_type", "" );
        aParameterlist.insert( "reinitialization_frequency", 1 );
        aParameterlist.insert( "output_mesh_file", "" );
        aParameterlist.insert( "time_offset", 0.0 );

        // SDF reinitialization mode: "field_remap" (default) or "interface"
        // "field_remap": Original mode - map solution field to target mesh
        // "interface": Compute exact SDF from XTK interface facets
        aParameterlist.insert( "sdf_reinit_mode", "field_remap" );

        // Material phase index for sign determination in interface mode
        // Nodes in this phase get negative SDF (inside), others get positive (outside)
        aParameterlist.insert( "sdf_reinit_material_phase", 0 );
    }

    //------------------------------------------------------------------------------

    /**
     * Parameters of the topological-derivative design operator (wrk::TD_Design_Operator), MORISGENERAL sublist
     * TD_OPERATOR. The operator maps per-cell topological-derivative values g_e (ELEMENTAL FEM fields) to the design
     * B-spline space of a level-set field, b_A = sum_e g_e int_e B_A, and solves the Helmholtz filter
     * ( gamma K + M ) g~ = b, with M the B-spline mass matrix over the domain and K its Laplacian (natural BC).
     * Usage in a deck:
     *   aParameterLists( MORISGENERAL_Submodule::TD_OPERATOR ).add_parameter_list( prm::create_td_operator_parameter_list() );
     *   aParameterLists.set( "adv_field", "<GEN design name>" );
     */
    inline Parameter_List
    create_td_operator_parameter_list()
    {
        Parameter_List tParameterList( "TD_Operator" );

        // name of the GEN design (geometry or property) whose B-spline coefficients are the level-set ADVs
        tParameterList.insert( "adv_field", "" );

        // comma-separated ELEMENTAL FEM fields (field_elemental_normalization = ip_cell), summed per background cell
        // to g_e, e.g. one topological-derivative field per phase
        tParameterList.insert( "td_field_names", "" );

        // optional ELEMENTAL ip_cell field holding the cut ratio rho_e (a VOLUME IQI on the material sets); when set,
        // ( 1 - rho_e ) * void_td_value is added to g_e (the void value of a model that does not mesh the void)
        tParameterList.insert( "cut_ratio_field", "" );
        tParameterList.insert( "void_td_value", 0.0 );

        // Helmholtz filter coefficient gamma (length scale sqrt( gamma )); 0 gives the L2 projection M^-1 b
        tParameterList.insert( "filter_gamma", 1.0e-4 );

        // optional hdf5 base name: every evaluation writes <base>_<n>.hdf5 (cell values, load, direction, ADVs)
        tParameterList.insert( "output_file", "" );

        return tParameterList;
    }

    //------------------------------------------------------------------------------

}    // namespace moris::prm
