# CMakeList.txt : CMake project for PCFT, include source and define
# project specific logic here.
#
find_package(FFTW3)
if(FFTW3_FOUND)
    add_compile_definitions(USING_FFTW3)
endif()

# Add source to this project's executable.
set(PCFT_LIB "PCFT")
add_library (${PCFT_LIB})

set(PCFT_EXE "PCFTexe")
add_executable(${PCFT_EXE} PCFT.cpp)

# TODO: Add tests and install targets if needed.
target_link_libraries(${PCFT_LIB} PRIVATE FFTW3::fftw3)
target_include_directories(${PCFT_LIB} PUBLIC ${PROJECT_SOURCE_DIR}/include)

target_link_libraries(${PCFT_EXE} PRIVATE ${PCFT_LIB})

# Add subdirectories containing code
add_subdirectory("numerics")
add_subdirectory("tools")