# CMakeList.txt : CMake project for PCFT, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

find_package(FFTW3 CONFIG REQUIRED)
find_package(OpenMP REQUIRED)

# Add source to this project's executable.
add_executable (PCFT 
"PCFT.cpp"
"PCFT.hpp"
"numerics/FFTWTransformer.cpp"
"numerics/FFTWTransformer.hpp"
"numerics/NaiveTransformer.hpp"
"numerics/NaiveTransformer.cpp"
"numerics/IFourierTransformer.hpp"
"tools/Timer.hpp")

# TODO: Add tests and install targets if needed.
target_link_libraries(PCFT PRIVATE FFTW3::fftw3)
target_link_libraries(PCFT PRIVATE OpenMP::OpenMP_CXX)