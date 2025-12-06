cmake_minimum_required(VERSION 3.8)

# Create a library from simulation-related source files
add_library(simulation boundary_manager.cpp simulation.cpp)

# Specify link libraries for this library
target_link_libraries(${PROJECT_NAME} simulation)
