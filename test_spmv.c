#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fenv.h>
#include <algorithm>

#include<mkl.h>
#include "mkl_types.h"
#include "mkl_spblas.h"
#include "omp.h"

#include "loadMatrixMarket.h"



static const size_t LLC_CAPACITY = 32*1024*1024;
static const double *bufToFlushLlc = NULL;


static void printPerf(double *times, int REPEAT)
{
    //sort(times, times + REPEAT);
    double sum = 0;
    for (int i  = 0; i < REPEAT; i++){
        sum += times[i];
    }
    double t = sum/(REPEAT);

    printf("\t%.6f (us) \n",  t*1e6);
}

void flushLlc()
{
  double sum = 0;
#pragma omp parallel for reduction(+:sum)
  for (size_t i = 0; i < LLC_CAPACITY/sizeof(bufToFlushLlc[0]); ++i) {
    sum += bufToFlushLlc[i];
  }
  FILE *fp = fopen("/dev/null", "w");
  fprintf(fp, "%f\n", sum);
  fclose(fp);
}

int main(int argc, char **argv)
{
  if (argc < 2) {
    fprintf(
            stderr,
            "-- Usage examples --\n"
            "  %s inline_1.mtx type: run with inline_1 matrix in matrix market format\n",
            argv[0]);
    return -1;
  }
  
  
  CSR A;
  int outputbase = 1;
  loadMatrixMarket(argv[1], &A, outputbase, 0 /*transpose =false*/);
  int nnz = A.nnz;
  double bytes = sizeof(int)*(nnz  + A.m) + //rowprt + colidx
    sizeof(double)*nnz +   // values
    sizeof(double)*2*A.m; // vectors
  printf("m = %d nnz = %d %f bytes = %f\n", A.m, nnz, (double)nnz/A.m, bytes);

  
  double *x = (double*)malloc(sizeof(double)*A.m);
  double *y = (double*)malloc(sizeof(double)*A.m);

  // Random x
  int seed[] = {0, 0, 0, 1};
  LAPACKE_dlarnv(1, seed, (size_t)A.m, x);
  
  // allocate a large buffer to flush out cache
  int nbFlush = 2048;
  bufToFlushLlc = (double *)_mm_malloc(LLC_CAPACITY, 64);

  const int REPEAT = 5;
  double times[REPEAT];
  sparse_status_t stat;

/******************************** Preparing for Double precision*****************/
  
  // Inspector-Executor interface in MKL 11.3+
  sparse_matrix_t mklA;
  stat = mkl_sparse_d_create_csr(
    &mklA,
    SPARSE_INDEX_BASE_ZERO, A.m, A.n,
    A.rowptr, A.rowptr + 1,
    A.colidx, A.values);

  if (SPARSE_STATUS_SUCCESS != stat) {
    fprintf(stderr, "Failed to create mkl csr\n");
    return -1;
  }

  matrix_descr descA;
  descA.type = SPARSE_MATRIX_TYPE_GENERAL;
  descA.diag = SPARSE_DIAG_NON_UNIT;

/******************************END Preparation for Double precision********************************/
  

/******************************Preparing for single  SINGLE************************************/
  // Single precision variant
  double RMAX = (double) LAPACKE_slamch('O');
  int overflow = 0;

  // convert fp64 matrix to fp32
  
  float *x32 = (float*)malloc(sizeof(float)*A.m);
  float *y32 = (float*)malloc(sizeof(float)*A.m);
  double *y32to64 = (double*)malloc(sizeof(double)*A.m);
  float *values32 = (float*)malloc(sizeof(float)*nnz);

  #pragma omp parallel for
  for (int k = 0; k < nnz; k++){
      if (A.values[k] > RMAX) {
          overflow = 1;
      }else{
          values32[k] = (float) (A.values[k]);
      }
  }

  if (overflow) {
      printf("Conversion of A overflow\n");
      overflow = 0;
  }
  

  #pragma omp parallel for
  for (int k = 0; k < A.m; k++) {
      if (A.values[k] > RMAX) {
          overflow = 1;
      } else {
          x32[k] = (float)x[k];
      }
  }
  
  if (overflow) {
      printf("Conversion of X overflow\n");
      overflow = 0;
  }

  sparse_matrix_t mklA32;
  stat = mkl_sparse_s_create_csr(
      &mklA32,
      SPARSE_INDEX_BASE_ZERO, A.m, A.n,
      A.rowptr, A.rowptr + 1,
      A.colidx, values32);

  if (SPARSE_STATUS_SUCCESS != stat) {
      fprintf(stderr, "Failed to create mkl csr\n");
      return -1;
  }

  matrix_descr descA32;
  descA32.type = SPARSE_MATRIX_TYPE_GENERAL;
  descA32.diag = SPARSE_DIAG_NON_UNIT;
  
/*****************************END SINGLE*****************************************/
  
/******************* Run double and single***************************************/
  for (int i = 0; i < REPEAT; ++i) {
    for (int j = 0; j < nbFlush; ++j) flushLlc();
    double t = omp_get_wtime();
    mkl_sparse_d_mv(
      SPARSE_OPERATION_NON_TRANSPOSE, 1, mklA, descA, x, 0, y);
    times[i] = omp_get_wtime() - t;
  }

  printf("Double precision SpMV time:");
  printPerf(times, REPEAT);

  for (int i = 0; i < REPEAT; ++i) {
      for (int j = 0; j < nbFlush; ++j) flushLlc();
      double t = omp_get_wtime();
      mkl_sparse_s_mv(
          SPARSE_OPERATION_NON_TRANSPOSE, 1, mklA32, descA32, x32, 0, y32);
      times[i] = omp_get_wtime() - t;
  }

  printf("Single precision SpMV time:");
  printPerf(times, REPEAT);
/*****************************END RUNS******************************************/

  // Copy solution back to fp64

  #pragma omp parallel for 
  for (int k = 0; k < A.m; k++) {
      y32to64[k] = (double) y32[k];
  }


  #pragma omp parallel for 
  for (int k = 0; k < nnz; k++) {
      A.values[k] = (double) values32[k];
  }

  // Check results
  double error = 0.0;
  #pragma omp parallel for
  for(int i = 0; i < A.m; i++) {
      error += abs(y[i] - y32to64[i]);
  }

  printf ("Error= %e\n", error/A.m); 
  free(x);
  free(y);
  free(x32);
  free(y32);

  free(y32to64);
  free(values32);
}
