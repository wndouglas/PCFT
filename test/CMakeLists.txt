find_package(catch2)

if(catch2_FOUND)
    enable_testing()

    set(PCFT_TEST "PCFTTest")
    add_executable (${PCFT_TEST} main.cpp)

    target_link_libraries(${PCFT_TEST} ${PCFT_LIB})
    target_link_libraries(${PCFT_TEST} Catch2::Catch2)
    target_include_directories(${PCFT_TEST} PUBLIC ${PROJECT_SOURCE_DIR}/include)

endif()