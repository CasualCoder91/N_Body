cmake_minimum_required(VERSION 3.10.2)
project (N_Body)

# Add the cmake folder so the FindSphinx module is found
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

#checks if OpenMP is available
find_package(OpenMP REQUIRED)

# set everything up for c++ 17 features
set(CMAKE_CXX_STANDARD 17)
# Don't add this line if you will try_compile with boost.
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SOURCES src/Analysis.cpp src/Configuration.cpp src/Constants.cpp src/Database.cpp src/InitialConditions.cpp src/InOut.cpp
src/Integrator.cpp src/LookupTable.cpp src/main.cpp src/Matrix.cpp src/MWPotential.cpp src/Node.cpp src/Parameters.cpp
src/Plot.cpp src/Point src/ProgressBar.cpp src/Projection.cpp src/Simulation.cpp src/Star.cpp src/Vec3D.cpp
src/VDBSCAN.cpp src/Test.cpp)


include_directories(include)
#include_directories(lib)
#add_subdirectory(lib)

#add_executable(main main.cpp)
add_executable(N_Body ${SOURCES})

SET( LIBS OpenMP::OpenMP_CXX sqlite )
target_link_libraries(N_Body ${LIBS}  )

add_subdirectory (docs)

configure_file("${PROJECT_SOURCE_DIR}/docs/stars0.jpg" "${PROJECT_SOURCE_DIR}/build/docs/stars0.jpg" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/docs/conf.py" "${PROJECT_SOURCE_DIR}/build/docs/conf.py" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/sqlite3.dll" "${PROJECT_SOURCE_DIR}/build/sqlite3.dll" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/sqlite3.lib" "${PROJECT_SOURCE_DIR}/build/sqlite3.lib" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/sqlite3.dll" "${PROJECT_SOURCE_DIR}/build/Release/sqlite3.dll" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/sqlite3.lib" "${PROJECT_SOURCE_DIR}/build/Release/sqlite3.lib" COPYONLY)

configure_file("${PROJECT_SOURCE_DIR}/src/Plots/potentialCircularVelocity.py" "${PROJECT_SOURCE_DIR}/docs/pyplots/potentialCircularVelocity.py" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/src/Plots/initialConditionsMassBulge.py" "${PROJECT_SOURCE_DIR}/docs/pyplots/initialConditionsMassBulge.py" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/src/Plots/potentialPositions.py" "${PROJECT_SOURCE_DIR}/docs/pyplots/potentialPositions.py" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/src/Plots/massDistribution.py" "${PROJECT_SOURCE_DIR}/docs/pyplots/massDistribution.py" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/src/Plots/potentialSurfaceDensity.py" "${PROJECT_SOURCE_DIR}/docs/pyplots/potentialSurfaceDensity.py" COPYONLY)

configure_file("${PROJECT_SOURCE_DIR}/src/Plots/initialConditionsMassSalpeter.py" "${PROJECT_SOURCE_DIR}/docs/pyplots/initialConditionsMassSalpeter.py" COPYONLY)

configure_file("${PROJECT_SOURCE_DIR}/src/Plots/initialConditionsMassDisk.py" "${PROJECT_SOURCE_DIR}/docs/pyplots/initialConditionsMassDisk.py" COPYONLY)


file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/build/Release/Output/Database)
