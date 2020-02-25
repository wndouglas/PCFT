find_package(catch2)

if(catch2_FOUND)
    set(PCFT_TEST "PCFTTest")
    add_executable (${PCFT_TEST} src/main.cpp src/numerics/numericsTests.cpp)

    target_link_libraries(${PCFT_TEST} ${PCFT_LIB})
    target_link_libraries(${PCFT_TEST} Catch2::Catch2)
    target_include_directories(${PCFT_TEST} PUBLIC ${PCFT_INCLUDE_DIR})
endif()