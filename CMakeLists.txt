cmake_minimum_required(VERSION 3.8)
#project(numsim)
add_subdirectory(src)

# generate files used by clangd 
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
