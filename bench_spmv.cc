#include "mkl_spblas.h"
#include <mkl.h>

#include "bench_consts.h"
#include "loadMatrixMarket.h"
#include <benchmark/benchmark.h>
#include <iostream>
#include <string>
#include <vector>

CSR A;
sparse_matrix_t mklA_d, mklA_s;
sparse_status_t status = SPARSE_STATUS_EXECUTION_FAILED;
matrix_descr descA;
std::vector<double> x_d, y_d;
std::vector<float> x_s, y_s;

static void DoSetupDouble(const benchmark::State &state) {
  const std::string MAT_PATH =
      BENCH_CONSTS::MAT_DIR + BENCH_CONSTS::MATS[state.range(0)] + ".mtx";
  loadMatrixMarket(MAT_PATH.c_str(), &A, 0, 0);

  x_d.resize(A.m);
  y_d.resize(A.m);
  int seed[] = {0, 0, 0, 1};
  LAPACKE_dlarnv(1, seed, A.m, x_d.data());

  status = mkl_sparse_d_create_csr(&mklA_d, SPARSE_INDEX_BASE_ZERO, A.m, A.n,
                                   A.rowptr, A.rowptr + 1, A.colidx, A.values);
  descA.type = SPARSE_MATRIX_TYPE_GENERAL;
  descA.diag = SPARSE_DIAG_NON_UNIT;
}

static void DoTeardownDouble(const benchmark::State &state) {
  x_d.clear();
  y_d.clear();
  mkl_sparse_destroy(mklA_d);
}

static void BMMatrixDouble(benchmark::State &state) {
  if (status != SPARSE_STATUS_SUCCESS) {
    state.SkipWithError("failed in creation of mkl csr!");
    return;
  }
  for (auto _ : state) {
    mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, mklA_d, descA,
                    x_d.data(), 0, y_d.data());
  }
}

BENCHMARK(BMMatrixDouble)
    ->DenseRange(0, 36, 1)
    ->Unit(benchmark::kMillisecond)
    ->Setup(DoSetupDouble)
    ->Teardown(DoTeardownDouble);

static void DoSetupSingle(const benchmark::State &state) {
  const std::string MAT_PATH =
      BENCH_CONSTS::MAT_DIR + BENCH_CONSTS::MATS[state.range(0)] + ".mtx";
  loadMatrixMarket(MAT_PATH.c_str(), &A, 0, 0);

  x_d.resize(A.m);
  y_d.resize(A.m);
  int seed[] = {0, 0, 0, 1};
  LAPACKE_dlarnv(1, seed, A.m, x_d.data());

  const double RMAX = static_cast<double>(LAPACKE_slamch('O'));
  x_s.resize(A.m);
  y_s.resize(A.m);
  std::vector<float> val_s(A.nnz);

  for (int i = 0; i < A.nnz; i++) {
    if (A.values[i] <= RMAX) {
      val_s[i] = static_cast<float>(A.values[i]);
    }
  }

  for (int i = 0; i < A.m; i++) {
    if (A.values[i] <= RMAX) {
      x_s[i] = static_cast<float>(x_d[i]);
    }
  }

  status =
      mkl_sparse_s_create_csr(&mklA_s, SPARSE_INDEX_BASE_ZERO, A.m, A.n,
                              A.rowptr, A.rowptr + 1, A.colidx, val_s.data());
  if (status != SPARSE_STATUS_SUCCESS) {
    std::cout << "matrix creation failed" << std::endl;
  }
  descA.type = SPARSE_MATRIX_TYPE_GENERAL;
  descA.diag = SPARSE_DIAG_NON_UNIT;
}

static void DoTeardownSingle(const benchmark::State &state) {
  x_d.clear();
  y_d.clear();
  x_s.clear();
  y_s.clear();
  mkl_sparse_destroy(mklA_s);
}

static void BMMatrixSingle(benchmark::State &state) {
  if (status != SPARSE_STATUS_SUCCESS) {
    state.SkipWithError("failed in creation of mkl csr!");
    return;
  }
  for (auto _ : state) {
    mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, mklA_s, descA,
                    x_s.data(), 0, y_s.data());
  }
}

BENCHMARK(BMMatrixSingle)
    ->DenseRange(0, 36, 1)
    ->Unit(benchmark::kMillisecond)
    ->Setup(DoSetupSingle)
    ->Teardown(DoTeardownSingle);

BENCHMARK_MAIN();
