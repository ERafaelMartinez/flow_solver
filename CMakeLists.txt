cmake_minimum_required(VERSION 3.8)

# Define the project name.
project(numsim)

# generate files used by clangd 
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_subdirectory(src)
