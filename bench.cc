#include <benchmark/benchmark.h>
#include <iostream>

#include "bench_consts.h"

static void DoSetup(const benchmark::State &state) {
  const std::string MAT_PATH =
      BENCH_CONSTS::MAT_DIR + BENCH_CONSTS::MATS[state.range(0)];
  std::cout << "setting up: " << MAT_PATH << std::endl;
}

static void DoTeardown(const benchmark::State &state) {
  const std::string MAT_PATH =
      BENCH_CONSTS::MAT_DIR + BENCH_CONSTS::MATS[state.range(0)];
  std::cout << "tearing down: " << MAT_PATH << std::endl;
}

static void BM_StringCreation(benchmark::State &state) {
  for (auto _ : state) {
    std::string empty_string;
  }
}
// Register the function as a benchmark
BENCHMARK(BM_StringCreation)
    ->DenseRange(0, 12, 1)
    ->Threads(1)
    ->Setup(DoSetup)
    ->Teardown(DoTeardown);

// template <class... Args>
// void BM_takes_args(benchmark::State &state, Args &&...args) {
//   auto args_tuple = std::make_tuple(std::move(args)...);
//   for (auto _ : state) {
//     std::string empty_string;
//   }
// }
// BENCHMARK_CAPTURE(BM_takes_args, 2cubes_sphere, std::string("2cubes_sphere"))
//     ->Setup(DoSetup)
//     ->Teardown(DoTeardown);

BENCHMARK_MAIN();
