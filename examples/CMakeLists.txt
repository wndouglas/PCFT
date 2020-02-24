# CMakeList.txt : CMake project for PCFT, include source and define
# project specific logic here.
#

MESSAGE("${PCFT_LIB} out")

# Add source to this project's executable.
add_executable (Example1 PCFT.cpp)
target_link_libraries(Example1 PRIVATE ${PCFT_LIB})
target_include_directories(Example1 PRIVATE ${PROJECT_SOURCE_DIR}/include)