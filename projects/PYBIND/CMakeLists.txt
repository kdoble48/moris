#
# Minimal pybind11 bindings for MORIS workflows
#

cmake_minimum_required(VERSION 3.17)

if (NOT BUILD_PYBINDINGS)
    return()
endif()

# Prefer the active Spack Python over system fallbacks.
set(PYBIND11_FINDPYTHON ON)
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 QUIET)

if (NOT pybind11_FOUND)
    message(FATAL_ERROR
        "BUILD_PYBINDINGS is ON, but pybind11 was not found. "
        "Install pybind11 or rerun with -DBUILD_PYBINDINGS=OFF.")
endif()

# Include dirs similar to main executable
set(PYBIND_INCLUDES
    ${COM}
    ${IOS}
    ${WRK}
    ${LINALG}
)

foreach(PYBIND_INCLUDE ${PYBIND_INCLUDES})
    include_directories(../${PYBIND_INCLUDE}/src)
endforeach()

set(PYBIND_DEPENDENCIES
    ${HMR}-lib
    ${COM}-lib
    ${IOS}-lib
    ${GEN}-lib
    ${WRK}-lib
    ${MORIS_PETSC_LIBS}
    ${MORIS_BOOST_LIBS}
    ${MORIS_ACML_LAPACK_MKL_OPENBLAS_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_BASE_LIBS}
)

pybind11_add_module(moris_py moris_py.cpp)
target_link_libraries(moris_py PRIVATE ${PYBIND_DEPENDENCIES})
target_link_libraries(moris_py PRIVATE ${MORIS}::trilinos)

set_target_properties(moris_py PROPERTIES OUTPUT_NAME "moris_py")
set_target_properties(moris_py PROPERTIES
    CXX_VISIBILITY_PRESET "default"
    VISIBILITY_INLINES_HIDDEN OFF)
