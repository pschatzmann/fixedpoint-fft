cmake_minimum_required(VERSION 3.14)
project(fixedpoint-fft CXX)

# Header-only library: just exposes src/ as an include path so other
# CMake projects can do
#   add_subdirectory(path/to/fixedpoint-fft)
#   target_link_libraries(your_target PRIVATE fixedpoint-fft)
add_library(fixedpoint-fft INTERFACE)
target_include_directories(fixedpoint-fft INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_compile_features(fixedpoint-fft INTERFACE cxx_std_17)

# Only build the host-side test suite when this is the top-level
# project (i.e. not when pulled in via add_subdirectory() by someone
# else), unless explicitly requested.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    option(FIXEDPOINT_FFT_BUILD_TESTS "Build the fixedpoint-fft host-side test suite" ON)
else()
    option(FIXEDPOINT_FFT_BUILD_TESTS "Build the fixedpoint-fft host-side test suite" OFF)
endif()

if(FIXEDPOINT_FFT_BUILD_TESTS)
    enable_testing()
    add_executable(fixedpoint_fft_test test/host_test.cpp)
    target_link_libraries(fixedpoint_fft_test PRIVATE fixedpoint-fft)
    add_test(NAME fixedpoint_fft_test COMMAND fixedpoint_fft_test)
endif()
