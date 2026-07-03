/*
 * Minimal pybind11 module to run MORIS workflows from Python.
 * Usage (multi-rank): mpirun -np 4 python -c "import moris_py; moris_py.run_so('./my_case.so')"
 */

#include <algorithm>
#include <cctype>
#include <memory>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "cl_Communication_Manager.hpp"
#include "cl_Logger.hpp"
#include "banner.hpp"

#include "cl_Library_Factory.hpp"
#include "cl_Library_IO.hpp"
#include "cl_WRK_Workflow_Factory.hpp"
#include "cl_WRK_Performer_Manager.hpp"
#include "cl_WRK_Workflow.hpp"
#include "cl_OPT_Manager.hpp"
#include "cl_Module_Parameter_Lists.hpp"
#include "cl_Submodule_Parameter_Lists.hpp"
#include "cl_Parameter_List.hpp"
#include "cl_Parameter.hpp"
#include "cl_Variant.hpp"
#include "cl_Design_Variable.hpp"

#ifdef MORIS_HAVE_ARBORX
#include <Kokkos_Core.hpp>
#endif

// Define MORIS globals for this shared library
moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

namespace py = pybind11;

struct RunResult {
    int                         rank   = 0;
    int                         size   = 1;
    int                         ret    = 0;
    std::string                 parameter_receipt;
    std::vector< moris::real >  iqi; // filled only for non-optimization runs
};

namespace
{
    std::string normalize_token( std::string token )
    {
        std::string result;
        result.reserve( token.size() );
        for ( unsigned char ch : token )
        {
            if ( std::isalnum( ch ) )
            {
                result.push_back( static_cast< char >( std::toupper( ch ) ) );
            }
        }
        return result;
    }

    moris::Module_Type module_from_string( const std::string& name )
    {
        const std::string normalized = normalize_token( name );
        for ( uint i = 0; i < static_cast< uint >( moris::Module_Type::END_ENUM ); ++i )
        {
            moris::Module_Type module = static_cast< moris::Module_Type >( i );
            if ( module == moris::Module_Type::END_ENUM )
            {
                continue;
            }
            if ( normalize_token( moris::convert_parameter_list_enum_to_string( module ) ) == normalized )
            {
                return module;
            }
        }

        std::ostringstream oss;
        oss << "Unknown MORIS module '" << name << "'.";
        throw pybind11::value_error( oss.str() );
    }

    uint submodule_index_from_string( moris::Module_Type module, const std::string& name )
    {
        const std::string normalized = normalize_token( name );
        moris::Vector< std::string > submodules = moris::get_submodule_names( module );

        for ( uint idx = 0; idx < submodules.size(); ++idx )
        {
            if ( normalize_token( submodules( idx ) ) == normalized )
            {
                return idx;
            }
        }

        std::ostringstream oss;
        oss << "Unknown submodule '" << name << "' for module '"
            << moris::convert_parameter_list_enum_to_string( module ) << "'.";
        throw pybind11::value_error( oss.str() );
    }

    template< typename T >
    T convert_py_value( const pybind11::handle& value, const std::string& context )
    {
        try
        {
            if constexpr ( std::is_same_v< T, bool > )
            {
                return value.cast< bool >();
            }
            else if constexpr ( std::is_same_v< T, moris::uint > )
            {
                return static_cast< moris::uint >( value.cast< unsigned long long >() );
            }
            else if constexpr ( std::is_same_v< T, moris::sint > )
            {
                return static_cast< moris::sint >( value.cast< long long >() );
            }
            else if constexpr ( std::is_same_v< T, moris::real > )
            {
                return value.cast< moris::real >();
            }
            else if constexpr ( std::is_same_v< T, std::string > )
            {
                return value.cast< std::string >();
            }
            else if constexpr ( std::is_same_v< T, std::pair< std::string, std::string > > )
            {
                pybind11::sequence seq = value.cast< pybind11::sequence >();
                if ( seq.size() != 2 )
                {
                    throw pybind11::value_error( context + ": expected a pair of two strings" );
                }
                return std::make_pair( seq[ 0 ].cast< std::string >(), seq[ 1 ].cast< std::string >() );
            }
            else if constexpr ( std::is_same_v< T, moris::Vector< moris::uint > > )
            {
                std::vector< moris::uint > data = value.cast< std::vector< moris::uint > >();
                moris::Vector< moris::uint > result;
                result = data;
                return result;
            }
            else if constexpr ( std::is_same_v< T, moris::Vector< moris::sint > > )
            {
                std::vector< moris::sint > data = value.cast< std::vector< moris::sint > >();
                moris::Vector< moris::sint > result;
                result = data;
                return result;
            }
            else if constexpr ( std::is_same_v< T, moris::Vector< moris::real > > )
            {
                if ( pybind11::isinstance< pybind11::float_ >( value ) )
                {
                    moris::Vector< moris::real > single( 1, value.cast< moris::real >() );
                    return single;
                }
                std::vector< moris::real > data = value.cast< std::vector< moris::real > >();
                moris::Vector< moris::real > result;
                result = data;
                return result;
            }
            else if constexpr ( std::is_same_v< T, moris::Vector< std::string > > )
            {
                std::vector< std::string > data = value.cast< std::vector< std::string > >();
                moris::Vector< std::string > result;
                result = data;
                return result;
            }
            else if constexpr ( std::is_same_v< T, moris::Design_Variable > )
            {
                if ( value.is_none() )
                {
                    throw pybind11::value_error( context + ": Design_Variable override cannot be None" );
                }

                if ( pybind11::isinstance< pybind11::dict >( value ) )
                {
                    pybind11::dict dict = value.cast< pybind11::dict >();

                    if ( !( dict.contains( "value" ) && dict.contains( "lower" ) && dict.contains( "upper" ) ) )
                    {
                        throw pybind11::value_error( context +
                                ": Design_Variable dict must contain 'lower', 'value', and 'upper' keys" );
                    }

                    moris::real lower = dict[ "lower" ].cast< moris::real >();
                    moris::real current = dict[ "value" ].cast< moris::real >();
                    moris::real upper = dict[ "upper" ].cast< moris::real >();
                    std::string name = dict.contains( "name" ) ? dict[ "name" ].cast< std::string >() : "";

                    return moris::Design_Variable( lower, current, upper, name );
                }

                return moris::Design_Variable( value.cast< moris::real >() );
            }
            else
            {
                static_assert( moris::variant_index< T >() >= 0, "Unsupported variant type" );
                throw pybind11::value_error( context + ": override for this parameter type is not supported yet" );
            }
        }
        catch ( const pybind11::cast_error& )
        {
            std::ostringstream oss;
            oss << context << ": failed to convert override to required type";
            throw pybind11::type_error( oss.str() );
        }
    }

    void set_parameter_from_python(
            moris::Parameter_List&      list,
            const std::string&          param_name,
            const pybind11::handle&     value,
            const std::string&          context )
    {
        if ( !list.exists( param_name ) )
        {
            std::ostringstream oss;
            oss << context << ": parameter '" << param_name << "' does not exist";
            throw pybind11::value_error( oss.str() );
        }

        const moris::Variant& current = list.get_variant( param_name );

        std::visit(
                [&]( auto&& current_value ) {
                    using ValueType = std::decay_t< decltype( current_value ) >;
                    ValueType converted = convert_py_value< ValueType >( value, context + "." + param_name );
                    try
                    {
                        list.set( param_name, converted );
                    }
                    catch ( const std::runtime_error& err )
                    {
                        std::ostringstream oss;
                        oss << context << ": failed to set parameter '" << param_name << "': " << err.what();
                        throw pybind11::value_error( oss.str() );
                    }
                },
                current );
    }

    void apply_parameter_list_overrides(
            moris::Submodule_Parameter_Lists& submodule,
            moris::uint                       list_index,
            const pybind11::dict&             overrides,
            const std::string&                context )
    {
        if ( list_index >= submodule.size() )
        {
            std::ostringstream oss;
            oss << context << ": parameter list index " << list_index << " is out of range ("
                << submodule.size() << " available)";
            throw pybind11::value_error( oss.str() );
        }

        moris::Parameter_List& list = submodule( list_index );
        for ( auto item : overrides )
        {
            const std::string param_name = item.first.cast< std::string >();
            pybind11::handle value = item.second;
            if ( value.is_none() )
            {
                std::ostringstream oss;
                oss << context << "." << param_name << ": None overrides are not supported";
                throw pybind11::value_error( oss.str() );
            }

            set_parameter_from_python( list, param_name, value, context );
        }
    }

    void apply_submodule_overrides(
            moris::Module_Parameter_Lists& module_lists,
            moris::Module_Type             module,
            moris::uint                    submodule_index,
            const pybind11::handle&        overrides,
            const std::string&             module_label,
            const std::string&             submodule_label )
    {
        std::string context = module_label + "." + submodule_label;

        moris::Submodule_Parameter_Lists& submodule = module_lists( submodule_index );

        if ( pybind11::isinstance< pybind11::dict >( overrides ) )
        {
            pybind11::dict dict = overrides.cast< pybind11::dict >();
            bool treat_as_index_map = true;
            for ( auto item : dict )
            {
                if ( !pybind11::isinstance< pybind11::int_ >( item.first ) )
                {
                    treat_as_index_map = false;
                    break;
                }
            }

            if ( treat_as_index_map )
            {
                for ( auto item : dict )
                {
                    moris::uint list_index = item.first.cast< moris::uint >();
                    pybind11::handle payload = item.second;
                    if ( payload.is_none() )
                    {
                        continue;
                    }
                    pybind11::dict payload_dict = payload.cast< pybind11::dict >();
                    apply_parameter_list_overrides( submodule, list_index, payload_dict, context + "[" + std::to_string( list_index ) + "]" );
                }
            }
            else
            {
                apply_parameter_list_overrides( submodule, 0, dict, context + "[0]" );
            }
        }
        else if ( pybind11::isinstance< pybind11::sequence >( overrides ) )
        {
            pybind11::sequence seq = overrides.cast< pybind11::sequence >();
            const pybind11::ssize_t n = seq.size();
            for ( pybind11::ssize_t i = 0; i < n; ++i )
            {
                pybind11::handle entry = seq[ i ];
                if ( entry.is_none() )
                {
                    continue;
                }
                pybind11::dict entry_dict = entry.cast< pybind11::dict >();
                apply_parameter_list_overrides( submodule, static_cast< moris::uint >( i ), entry_dict, context + "[" + std::to_string( i ) + "]" );
            }
        }
        else
        {
            std::ostringstream oss;
            oss << context << ": Expected dict or sequence for overrides";
            throw pybind11::type_error( oss.str() );
        }
    }

    void apply_module_overrides(
            moris::Library_IO&       lib,
            moris::Module_Type       module,
            const pybind11::handle&  overrides,
            const std::string&       module_label )
    {
        moris::Vector< moris::Module_Parameter_Lists >& modules = lib.get_parameter_lists();
        moris::Module_Parameter_Lists& module_lists = modules( static_cast< moris::uint >( module ) );

        if ( module_lists.size() == 0 )
        {
            std::ostringstream oss;
            oss << module_label << ": module has no parameter lists in this Library";
            throw pybind11::value_error( oss.str() );
        }

        if ( pybind11::isinstance< pybind11::dict >( overrides ) )
        {
            pybind11::dict dict = overrides.cast< pybind11::dict >();
            for ( auto item : dict )
            {
                std::string submodule_name;
                moris::uint submodule_index = 0;
                if ( pybind11::isinstance< pybind11::int_ >( item.first ) )
                {
                    submodule_index = item.first.cast< moris::uint >();
                    submodule_name = "#" + std::to_string( submodule_index );
                }
                else
                {
                    submodule_name = item.first.cast< std::string >();
                    submodule_index = submodule_index_from_string( module, submodule_name );
                }

                apply_submodule_overrides( module_lists, module, submodule_index, item.second, module_label, submodule_name );
            }
        }
        else
        {
            std::ostringstream oss;
            oss << module_label << ": Expected dict for module overrides";
            throw pybind11::type_error( oss.str() );
        }
    }

    void apply_overrides( moris::Library_IO& lib, const pybind11::dict& overrides )
    {
        if ( overrides.empty() )
        {
            return;
        }

        for ( auto item : overrides )
        {
            std::string module_name = item.first.cast< std::string >();
            moris::Module_Type module = module_from_string( module_name );
            apply_module_overrides( lib, module, item.second, module_name );
        }
    }

    bool dry_run_requested()
    {
        if ( const char* flag = std::getenv( "MORIS_PY_DRY_RUN" ) )
        {
            return flag[ 0 ] != '\0' && flag[ 0 ] != '0';
        }
        return false;
    }
}

// Helper to build argc/argv for Logger initialization
static void build_argv(const std::vector<std::string>& args, int& argc, char**& argv, std::vector<std::string>& storage, std::vector<char*>& ptrs)
{
    storage.clear();
    ptrs.clear();
    storage.reserve(args.size() + 1);
    ptrs.reserve(args.size() + 1);
    storage.emplace_back("moris_py");
    for (const auto& s : args) storage.emplace_back(s);
    for (auto& s : storage) ptrs.push_back(const_cast<char*>(s.c_str()));
    argc = static_cast<int>(ptrs.size());
    argv = ptrs.data();
}

static RunResult run_impl(const std::string& so_path,
                          const std::string& xml_path,
                          bool                meshgen,
                          const std::vector<std::string>& logger_args,
                          const std::string& parameter_receipt,
                          const py::dict&    overrides)
{
    RunResult res;

    // Initialize MPI/Comm and Logger
    int argc = 0; char** argv = nullptr; std::vector<std::string> argv_storage; std::vector<char*> argv_ptrs;
    build_argv(logger_args, argc, argv, argv_storage, argv_ptrs);

    gMorisComm = moris::Comm_Manager(&argc, &argv);
    gLogger.initialize(argc, argv);

    // Optionally print banner (keep quiet by default)
    // moris::print_banner(argc, argv);

    moris::Library_Factory lf;
    std::shared_ptr< moris::Library_IO > lib = meshgen ? lf.create_Library(moris::Library_Type::MESHGEN)
                                                       : lf.create_Library(moris::Library_Type::STANDARD);

    if (!so_path.empty()) {
        lib->load_parameter_list(so_path, moris::File_Type::SO_FILE);
    }
    if (!xml_path.empty()) {
        lib->load_parameter_list(xml_path, moris::File_Type::XML_FILE);
    }

    if (!overrides.empty()) {
        apply_overrides(*lib, overrides);
    }

    lib->finalize(parameter_receipt);

    if ( dry_run_requested() )
    {
        res.rank = static_cast<int>(moris::par_rank());
        res.size = static_cast<int>(moris::par_size());
        res.ret  = 0;
        res.parameter_receipt = parameter_receipt;
        gMorisComm.finalize();
        return res;
    }

#ifdef MORIS_HAVE_ARBORX
    std::unique_ptr<Kokkos::ScopeGuard> guard = (!Kokkos::is_initialized() && !Kokkos::is_finalized())
        ? std::make_unique<Kokkos::ScopeGuard>() : nullptr;
#endif

    // Get OPT params and workflow
    moris::Module_Parameter_Lists optPL = lib->get_parameters_for_module(moris::Module_Type::OPT);
    MORIS_ERROR(optPL.size() > 0, "moris_py: OPT parameter list required but not set.");

    moris::wrk::Performer_Manager perf(lib);
    std::string wrk = optPL(0)(0).get< std::string >("workflow");
    moris::Vector< std::shared_ptr< moris::opt::Criteria_Interface > > flows( 1 );
    flows( 0 ) = moris::wrk::create_workflow( wrk, &perf );

    if (optPL(0)(0).get<bool>("is_optimization_problem")) {
        moris::opt::Manager mgr( optPL, flows );
        mgr.perform();
    } else {
        moris::Vector<moris::real> advs(0,0);
        moris::Vector<moris::real> bounds;
        moris::Matrix<moris::IdMat> dummy(1,1,0.0);
        flows( 0 )->initialize( advs, bounds, bounds, dummy );
        moris::Vector<moris::real> iqi = flows( 0 )->get_criteria( advs );
        const auto& iq_vals = iqi.data();
        res.iqi.assign( iq_vals.begin(), iq_vals.end() );
    }

    res.rank = static_cast<int>(moris::par_rank());
    res.size = static_cast<int>(moris::par_size());
    res.ret  = 0;
    res.parameter_receipt = parameter_receipt;

    gMorisComm.finalize();
    return res;
}

PYBIND11_MODULE(moris_py, m) {
    m.doc() = "Minimal MORIS Python bindings to run workflows using .so or .xml inputs";

    py::class_<RunResult>(m, "RunResult")
        .def_readonly("rank", &RunResult::rank)
        .def_readonly("size", &RunResult::size)
        .def_readonly("ret", &RunResult::ret)
        .def_readonly("parameter_receipt", &RunResult::parameter_receipt)
        .def_readonly("iqi", &RunResult::iqi);

    m.def("run_so",
        [](const std::string& so_path,
           const std::vector<std::string>& logger_args,
           bool meshgen,
           const std::string& parameter_receipt,
           py::object overrides) {
            py::dict override_dict = overrides.is_none() ? py::dict{} : overrides.cast< py::dict >();
            return run_impl(so_path, std::string{}, meshgen, logger_args, parameter_receipt, override_dict);
        },
        py::arg("so_path"),
        py::arg("logger_args") = std::vector<std::string>{},
        py::arg("meshgen") = false,
        py::arg("parameter_receipt") = std::string("./Parameter_Receipt.xml"),
        py::arg("overrides") = py::none(),
        R"pbdoc(Run a MORIS workflow from a shared object (.so) parameter file.

Notes:
- For multi-rank execution, launch Python under mpirun, e.g.:
  mpirun -np 4 python -c "import moris_py; moris_py.run_so('./case.so')"
- Pass Logger CLI args in logger_args, e.g. ["--outputlog", "run.log", "--severity", "1"]
- Provide optional overrides as nested dicts keyed by module/submodule names.
)pbdoc");

    m.def("run_xml",
        [](const std::string& xml_path,
           const std::vector<std::string>& logger_args,
           bool meshgen,
           const std::string& parameter_receipt,
           py::object overrides) {
            py::dict override_dict = overrides.is_none() ? py::dict{} : overrides.cast< py::dict >();
            return run_impl(std::string{}, xml_path, meshgen, logger_args, parameter_receipt, override_dict);
        },
        py::arg("xml_path"),
        py::arg("logger_args") = std::vector<std::string>{},
        py::arg("meshgen") = false,
        py::arg("parameter_receipt") = std::string("./Parameter_Receipt.xml"),
        py::arg("overrides") = py::none(),
        R"pbdoc(Run a MORIS workflow from an XML parameter file.

Notes:
- For multi-rank execution, launch Python under mpirun.
- Pass Logger CLI args in logger_args, e.g. ["--outputlog", "run.log", "--severity", "1"]
- Provide optional overrides as nested dicts keyed by module/submodule names.
)pbdoc");
}
