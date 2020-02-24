# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.13)

project (PCFTScheme CXX)
set(CMAKE_CXX_STANDARD 14)

# Include sub-projects.
add_subdirectory ("src")
add_subdirectory("test")
