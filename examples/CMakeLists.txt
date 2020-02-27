#
# CMakeList.txt : CMake project for example executables.
#

MESSAGE("${PCFT_LIB} out")

# Add source to this project's executable.
add_executable (Example1 PCFT.cpp)
target_link_libraries(Example1 PRIVATE ${PCFT_LIB})