cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O1")

project(ReSpaSol LANGUAGES C CXX)

set(MKL_ARCH intel64)
set(MKL_THREADING sequential)
set(MKL_INTERFACE lp64)
set(MKL_DIR "/opt/intel/oneapi/mkl")
find_package(MKL CONFIG REQUIRED PATHS $ENV{MKLROOT})
message(STATUS "Imported oneMKL targets: ${MKL_IMPORTED_TARGETS}")

option(BENCHMARK_ENABLE_GTEST_TESTS "BENCHMARK_ENABLE_GTEST_TESTS" OFF)
option(BENCHMARK_ENABLE_TESTING "BENCHMARK_ENABLE_TESTING" OFF)
add_subdirectory(benchmark)

add_subdirectory(ReadMatrixMarket)

add_executable(bench_spmv bench_spmv.cc)
target_compile_options(bench_spmv PUBLIC
  $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>
)
target_include_directories(bench_spmv PUBLIC
  $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>
  "${PROJECT_BINARY_DIR}"
  "${PROJECT_SOURCE_DIR}/ReadMatrixMarket"
)
target_link_libraries(bench_spmv PUBLIC
  $<LINK_ONLY:MKL::MKL>
  benchmark::benchmark
  libloadmatrix
)

# add_executable(test_spmv test_spmv.c)
# target_include_directories(test_spmv PUBLIC
#   $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>
#   "${PROJECT_BINARY_DIR}"
#   "${PROJECT_SOURCE_DIR}/ReadMatrixMarket"
# )
# target_link_libraries(test_spmv PUBLIC
#   $<LINK_ONLY:MKL::MKL>
#   libloadmatrix
# )
