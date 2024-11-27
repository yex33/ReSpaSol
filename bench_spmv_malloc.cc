#include "mkl_spblas.h"
#include <cstdlib>
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
double *x_d, *y_d;
float *x_s, *y_s;

static void DoSetupDouble(const benchmark::State &state) {
  const std::string MAT_PATH =
      BENCH_CONSTS::MAT_DIR + BENCH_CONSTS::MATS[state.range(0)] + ".mtx";
  loadMatrixMarket(MAT_PATH.c_str(), &A, 1, 0);

  std::cout << A.m << "x" << A.n << ", nnz = " << A.nnz << std::endl;
  x_d = (double *)malloc(sizeof(double) * A.m);
  y_d = (double *)malloc(sizeof(double) * A.m);
  int seed[] = {0, 0, 0, 1};
  LAPACKE_dlarnv(1, seed, A.m, x_d);

  status = mkl_sparse_d_create_csr(&mklA_d, SPARSE_INDEX_BASE_ZERO, A.m, A.n,
                                   A.rowptr, A.rowptr + 1, A.colidx, A.values);

  descA.type = SPARSE_MATRIX_TYPE_GENERAL;
  descA.diag = SPARSE_DIAG_NON_UNIT;
}

static void DoTeardownDouble(const benchmark::State &state) {
  free(x_d);
  free(y_d);
  mkl_sparse_destroy(mklA_d);
}

static void DoSetupSingle(const benchmark::State &state) {
  const std::string MAT_PATH =
      BENCH_CONSTS::MAT_DIR + BENCH_CONSTS::MATS[state.range(0)] + ".mtx";
  loadMatrixMarket(MAT_PATH.c_str(), &A, 1, 0);

  x_d = (double *)malloc(sizeof(double) * A.m);
  y_d = (double *)malloc(sizeof(double) * A.m);
  int seed[] = {0, 0, 0, 1};
  LAPACKE_dlarnv(1, seed, A.m, x_d);

  const double RMAX = static_cast<double>(LAPACKE_slamch('O'));
  x_s = (float *)malloc(sizeof(float) * A.m);
  y_s = (float *)malloc(sizeof(float) * A.m);
  float *val_s = (float *)malloc(sizeof(float) * A.nnz);

#pragma omp parallel for
  for (int i = 0; i < A.nnz; i++) {
    if (A.values[i] <= RMAX) {
      val_s[i] = static_cast<float>(A.values[i]);
    }
  }

#pragma omp parallel for
  for (int i = 0; i < A.m; i++) {
    if (A.values[i] <= RMAX) {
      x_s[i] = static_cast<float>(x_d[i]);
    }
  }

  status = mkl_sparse_s_create_csr(&mklA_s, SPARSE_INDEX_BASE_ZERO, A.m, A.n,
                                   A.rowptr, A.rowptr + 1, A.colidx, val_s);

  descA.type = SPARSE_MATRIX_TYPE_GENERAL;
  descA.diag = SPARSE_DIAG_NON_UNIT;
}

static void DoTeardownSingle(const benchmark::State &state) {
  free(x_d);
  free(y_d);
  free(x_s);
  free(y_s);
  mkl_sparse_destroy(mklA_s);
}

static void BMMatrixDouble(benchmark::State &state) {
  if (status != SPARSE_STATUS_SUCCESS) {
    state.SkipWithError("failed in creation of mkl csr!");
    return;
  }
  std::cout << "start" << std::endl;
  for (auto _ : state) {
    mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, mklA_d, descA, x_d, 0,
                    y_d);
  }
  std::cout << "end" << std::endl;
}

static void BMMatrixSingle(benchmark::State &state) {
  if (status != SPARSE_STATUS_SUCCESS) {
    state.SkipWithError("failed in creation of mkl csr!");
    return;
  }
  for (auto _ : state) {
    mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, mklA_s, descA, x_s, 0,
                    y_s);
  }
}

BENCHMARK(BMMatrixDouble)
    ->DenseRange(0, 36, 1)
    ->Unit(benchmark::kMillisecond)
    ->Setup(DoSetupDouble)
    ->Teardown(DoTeardownDouble);
BENCHMARK(BMMatrixSingle)
    ->DenseRange(0, 36, 1)
    ->Unit(benchmark::kMillisecond)
    ->Setup(DoSetupSingle)
    ->Teardown(DoTeardownSingle);

BENCHMARK_MAIN();
