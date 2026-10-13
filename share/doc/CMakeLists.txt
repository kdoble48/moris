#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# moris Doxygen -----------------------------------------------------------
# -------------------------------------------------------------------------

# Find Doxygen
find_package(Doxygen)

if(DOXYGEN_FOUND)
    # Input proper paths #> need to replace CMAKE_SOURCE_DIR in there
    configure_file(Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    
    # Create a doxygen make target
    add_custom_target(doxy
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating moris documentation with Doxygen."
        VERBATIM )
endif()

# moris documentation site (Sphinx) ---------------------------------------
# -------------------------------------------------------------------------
# Not part of ALL. `cmake --build <tree> --target docs` runs doc/Makefile and
# writes the site to <tree>/share/doc/site/html.

find_program(SPHINX_BUILD_EXECUTABLE NAMES sphinx-build)
find_program(MORIS_DOCS_MAKE_EXECUTABLE NAMES make gmake)

if(DOXYGEN_FOUND AND SPHINX_BUILD_EXECUTABLE AND MORIS_DOCS_MAKE_EXECUTABLE)
    add_custom_target(docs
        COMMAND ${MORIS_DOCS_MAKE_EXECUTABLE} -C ${PROJECT_SOURCE_DIR}/doc html
                DOXYGEN=${DOXYGEN_EXECUTABLE}
                SPHINXBUILD=${SPHINX_BUILD_EXECUTABLE}
                BUILDDIR=${CMAKE_CURRENT_BINARY_DIR}/site
        COMMENT "Building the moris documentation site with Doxygen and Sphinx."
        VERBATIM )
endif()

