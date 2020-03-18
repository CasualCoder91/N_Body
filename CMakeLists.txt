cmake_minimum_required(VERSION 3.10.2)
project (N_Body)

# Add the cmake folder so the FindSphinx module is found
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

#checks if OpenMP is available
find_package(OpenMP REQUIRED)

set(SOURCES src/main.cpp src/Analysis.cpp src/InitialConditions.cpp src/InOut.cpp src/Integrator.cpp src/Node.cpp src/Parameters.cpp src/Star.cpp src/Vec3D.cpp src/Configuration.cpp src/Database.cpp src/Simulation.cpp)


include_directories(include)
#include_directories(lib)
#add_subdirectory(lib)

#add_executable(main main.cpp)
add_executable(N_Body ${SOURCES})

SET( LIBS OpenMP::OpenMP_CXX sqlite3 )
target_link_libraries(N_Body ${LIBS}  )

add_subdirectory (docs)

configure_file("${PROJECT_SOURCE_DIR}/docs/stars0.jpg" "${PROJECT_SOURCE_DIR}/build/docs/stars0.jpg" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/docs/conf.py" "${PROJECT_SOURCE_DIR}/build/docs/conf.py" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/sqlite3.dll" "${PROJECT_SOURCE_DIR}/build/sqlite3.dll" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/sqlite3.lib" "${PROJECT_SOURCE_DIR}/build/sqlite3.lib" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/sqlite3.dll" "${PROJECT_SOURCE_DIR}/build/Release/sqlite3.dll" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/sqlite3.lib" "${PROJECT_SOURCE_DIR}/build/Release/sqlite3.lib" COPYONLY)
file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/build/Release/Output/Database)
