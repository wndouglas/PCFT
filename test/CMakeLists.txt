#
# CMakeList.txt : CMake project for test executables.
#
find_package(catch2)

if(catch2_FOUND)
    set(PCFT_TEST PCFTTest)
    add_executable (${PCFT_TEST} src/main.cpp src/numerics/numericsTest.cpp src/tools/toolsTest.cpp)

    target_link_libraries(${PCFT_TEST} PRIVATE ${PCFT_LIB})
    target_link_libraries(${PCFT_TEST} PRIVATE Catch2::Catch2)
    target_include_directories(${PCFT_TEST} PRIVATE ${PCFT_INCLUDE_DIR})
endif()