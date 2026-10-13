# Anatomy of an input deck

A MORIS input deck is one C++ file that fills a parameter list for each module and, where
needed, defines user callbacks: spatially varying properties, user-defined geometry fields,
output criteria, and the objective and constraint functions of an optimization problem. The file
is compiled to a shared object and loaded by the `moris` executable (see
{doc}`../getting-started/running`).

## The modules a deck configures

| Module | Parameter list configures |
|---|---|
| **OPT** | whether this is an optimization problem; the algorithm (GCMMA, MMA, ROL, LBFGS, SQP, or SWEEP) and its settings; restart file |
| **HMR** | the hierarchical B-spline background mesh: domain, element counts, refinement patterns, Lagrange and B-spline meshes |
| **XTK** | the immersed cut: decomposition, basis enrichment, ghost stabilization, integration-mesh output |
| **GEN** | geometry (level-set fields, field arrays, surface meshes), the phase table, design variables, and the list of criteria (`IQI_types`) |
| **FEM** | the physics, built bottom-up: properties → constitutive models → stabilization parameters → weak-form integrands (IWGs) → quantities of interest (IQIs) |
| **SOL** | the solver chain: linear algorithms and solvers, nonlinear (Newton, NLBGS, arc-length), time solvers, preconditioners |
| **MSI** | which degree-of-freedom types exist and which B-spline mesh discretizes each |
| **VIS** | the Exodus output: mesh type, sets, fields, file name |
| **MIG** | periodic side-set pairs (optional) |
| **MORISGENERAL** | general settings (optional) |

## The legacy style

A legacy deck defines one `extern "C"` function per module, which MORIS looks up by name in the
shared object:

```cpp
extern "C"
{
    void OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "is_optimization_problem", true );
        aParameterLists( OPT::ALGORITHMS ).add_parameter_list( opt::Optimization_Algorithm_Type::GCMMA );
        aParameterLists.set( "max_its", 100 );
    }

    void HMRParameterList( Module_Parameter_Lists& aParameterLists ) { /* … */ }
    void XTKParameterList( Module_Parameter_Lists& aParameterLists ) { /* … */ }
    // GEN, FEM, SOL, MSI, VIS, MIG, MORISGENERAL …

    // plus the OPT callbacks: compute_objectives, compute_constraints,
    // compute_dobjective_dadv, compute_dobjective_dcriteria, … and any
    // property / field functions referenced by name in the parameter lists
}
```

Within a module, `add_parameter_list()` opens a new entry in a sub-list (for example a new IWG
in `FEM::IWG`), and the following `set()` calls fill that entry. Calls that set an unknown key
abort at load time, so a typo in a key name is caught; a typo in a *value* that names
something in another module is not (see {doc}`interlinks`).

## The Deck API v2 style

The fork adds a single entry point. One function receives a `Deck` builder; modules are
activated by touching them, callbacks are registered in-process, and the objective and
constraints are written as expressions:

```cpp
MORIS_DECK( aDeck )
{
    aDeck.register_function( "Func_Neumann_U", &Func_Neumann_U );

    Module_Parameter_Lists& tHmr = aDeck.hmr();
    tHmr.set( "number_of_elements_per_dimension", 45, 15 );

    // …

    aDeck.objective( deck::criterion( "IQIStrainEnergy" ) / tInitialStrainEnergy );
    aDeck.constraint( deck::criterion( "IQIVolume" ) / tMaxMass - 1.0 <= 0.0 );
}
```

Both styles load through the same machinery and produce the same `Parameter_Receipt.xml`;
mixing them in one shared object is an error. See {doc}`deck-api-v2` for the full API.

## Reading an unfamiliar deck

Read the modules in dependency order rather than file order:

1. **GEN** — what geometry, how many phases (`phase_table`), and which criteria (`IQI_types`,
   whose order is the index order of the objective callbacks).
2. **FEM** — which physics on which phases and sets, and what each IQI computes.
3. **XTK** and **HMR** — how the mesh is built and cut.
4. **OPT** — algorithm and the objective/constraint definitions.
5. **SOL**, **MSI**, **VIS** — solver settings, DOF bookkeeping, output.
