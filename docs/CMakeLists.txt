find_package(Doxygen REQUIRED)
find_package(Sphinx REQUIRED)

# Find all the public headers
get_target_property(N_BODY_PUBLIC_HEADER_DIR N_Body INTERFACE_INCLUDE_DIRECTORIES)
file(GLOB_RECURSE N_BODY_PUBLIC_HEADERS ${N_BODY_PUBLIC_HEADER_DIR}/*.h)

set(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR})
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/docs/doxygen)
set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/html/index.html)
set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

message(STATUS "print_all_variables------------------------------------------{")
message(STATUS "DOXYGEN_INPUT_DIR=${DOXYGEN_INPUT_DIR}")
message(STATUS "DOXYGEN_OUTPUT_DIR=${DOXYGEN_OUTPUT_DIR}")
message(STATUS "DOXYGEN_INDEX_FILE=${DOXYGEN_INDEX_FILE}")
message(STATUS "DOXYFILE_IN=${DOXYFILE_IN}")
message(STATUS "DOXYFILE_OUT=${DOXYFILE_OUT}")
message(STATUS "print_all_variables------------------------------------------}")

#Replace variables inside @@ with the current values
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR}) #Doxygen won't create this for us

add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
                   DEPENDS ${N_BODY_PUBLIC_HEADERS}
                   COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
                   MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
                   COMMENT "Generating docs")

add_custom_target(Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})

set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/docs/sphinx)

add_custom_target(Sphinx ALL
                  COMMAND
                  ${SPHINX_EXECUTABLE} -b html -c ./
		  -Dbreathe_projects.N_Body=${DOXYGEN_OUTPUT_DIR}/xml
                  ${SPHINX_SOURCE} ${SPHINX_BUILD}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                  COMMENT "Generating documentation with Sphinx")
