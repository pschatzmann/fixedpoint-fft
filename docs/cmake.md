# Testing & CMake usage

This covers building/running the host-side test suite and consuming
this library from another CMake project - none of this is needed for
normal Arduino use (see [README.md](README.md) for that).

## Running the tests

`test/host_test.cpp` is a plain host-side (non-Arduino) program that
checks the twiddle table against `<cmath>`'s `cos()`/`sin()`, verifies
FFT peak-bin detection, and compares FFT/IFFT output against a
double-precision reference DFT. `<cmath>` usage there is fine - it's a
verification tool, not part of the library.

Via CMake/CTest:

```sh
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

Or directly with g++:

```sh
g++ -std=c++17 -O2 -Isrc test/host_test.cpp -o /tmp/fft_test && /tmp/fft_test
```

## Using the library from another CMake project

The library itself is header-only; `CMakeLists.txt` exposes it as an
`INTERFACE` target named `fixedpoint-fft` for other CMake projects to
consume (e.g. desktop tooling, simulation, or a non-Arduino embedded
build):

```cmake
add_subdirectory(path/to/fixedpoint-fft)
target_link_libraries(your_target PRIVATE fixedpoint-fft)
```

This only adds `src/` to your include path - it does not build the test
suite when pulled in this way (`FIXEDPOINT_FFT_BUILD_TESTS` defaults to
`OFF` for non-top-level use, `ON` when building this repo directly).
