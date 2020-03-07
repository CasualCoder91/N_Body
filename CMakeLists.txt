cmake_minimum_required(VERSION 3.10.2)
project (N_Body)

#checks if OpenMP is available
find_package(OpenMP REQUIRED)

include_directories(include)

set(SOURCES src/main.cpp src/Analysis.cpp src/InitialConditions.cpp src/InOut.cpp src/Integrator.cpp src/Node.cpp src/Parameters.cpp src/Star.cpp src/Vec3D.cpp)
#add_executable(main main.cpp)

add_executable(N_Body ${SOURCES}) 
target_link_libraries(N_Body PRIVATE OpenMP::OpenMP_CXX)
