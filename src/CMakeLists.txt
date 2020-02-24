# CMakeList.txt : CMake project for PCFT, include source and define
# project specific logic here.
#
find_package(FFTW3)
if(FFTW3_FOUND)
    add_compile_definitions(USING_FFTW3)
endif()

# Add source to this project's library.
set(PCFT_LIB "PCFT")
set(PCFT_LIB ${PCFT_LIB} PARENT_SCOPE)
add_library (${PCFT_LIB})

# TODO: Add tests and install targets if needed.
target_link_libraries(${PCFT_LIB} PRIVATE FFTW3::fftw3)
target_include_directories(${PCFT_LIB} PUBLIC ${PROJECT_SOURCE_DIR}/include)

# Add subdirectories containing code
add_subdirectory("numerics")
add_subdirectory("tools")
