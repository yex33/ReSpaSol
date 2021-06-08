#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <cusparse.h>
#include <cuda.h>
#include <mkl.h>
#include <mkl_spblas.h>
#include "loadMatrixMarket.h"


//#define FLOAT


// Define some error checking macros.
#define cudaErrCheck(stat) { cudaErrCheck_((stat), __FILE__, __LINE__); }

void cudaErrCheck_(cudaError_t stat, const char *file, int line) {
  if (stat != cudaSuccess) {
    fprintf(stderr, "CUDA Error: %s %s %d\n", cudaGetErrorString(stat), file, line);
  }
}


#define cusparseErrCheck(stat) { cusparseErrCheck_((stat), __FILE__, __LINE__); }
void cusparseErrCheck_(cusparseStatus_t stat, const char *file, int line) {
  if (stat != CUSPARSE_STATUS_SUCCESS) {
    fprintf(stderr, "CUSPARSE Error: %d %s %d\n", stat, file, line);
  }
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


    CSR matrixA;
    int outputbase = 0;
    loadMatrixMarket(argv[1], &matrixA, outputbase, 0 /*transpose =false*/);
    
    
    const int A_num_rows = matrixA.n;
    const int A_num_cols = matrixA.n;
    const int A_num_nnz  = matrixA.nnz;
    
    int  *hA_csrOffsets = matrixA.rowptr;
    int  *hA_columns  = matrixA.colidx;  
    
    // index pointer on device
    int   *dA_csrOffsets, *dA_columns;
    
#ifdef FLOAT
    float alpha = (float)1.0;
    float beta = (float)0.0;
    float *hA_values = (float*)malloc(A_num_nnz*sizeof(float));
    #pragma omp parallel for
    for(int i =0; i < A_num_nnz; i++){
      hA_values[i] = (float) matrixA.values[i];   
    }
    float *Y        = (float*)malloc(A_num_rows*sizeof(float)); // result
    float *X        = (float*)malloc(A_num_cols*sizeof(float));
    #pragma omp parallel for
    for (int i = 0; i < A_num_cols; i++){ X[i] = (float) 1.0;}
    // device 
    float *dA_values;
    float *dX;
    float *dY;
#else
    double alpha = (double)1.0;
    double beta =  (double)0.0;
    double *hA_values = matrixA.values;
    double *Y        = (double*)malloc(A_num_rows*sizeof(double));
    double *X        = (double*)malloc(A_num_cols*sizeof(double));
     #pragma omp parallel for
    for (int i = 0; i < A_num_cols; i++){ X[i] = (double)1.0;}
    //device 
    double *dA_values;
    double *dX;
    double *dY;
#endif
    
    // Allocate device memory to store the sparse CSR representation of A
    cudaErrCheck(cudaMalloc((void**) &dA_csrOffsets, (A_num_rows + 1) * sizeof(int)));
    cudaErrCheck(cudaMalloc((void **)&dA_columns, A_num_nnz*sizeof(int)));
    
#ifdef FLOAT
    cudaErrCheck(cudaMalloc((void **)&dA_values, A_num_nnz*sizeof(float)));
#else
    cudaErrCheck(cudaMalloc((void **)&dA_values, A_num_nnz*sizeof(double)));
#endif 

    // Allocate device memory to store the X and Y
#ifdef FLOAT
    cudaErrCheck(cudaMalloc((void **)&dX, A_num_cols*sizeof(float)));
    cudaErrCheck(cudaMalloc((void **)&dY, A_num_rows*sizeof(float)));
#else
    cudaErrCheck(cudaMalloc((void **)&dX, A_num_cols*sizeof(double)));
    cudaErrCheck(cudaMalloc((void **)&dY, A_num_rows*sizeof(double)));
#endif 

    // transfer data to device 
    // Transfer the input vectors and dense matrix A to the device
    cudaErrCheck(cudaMemcpy(dA_csrOffsets, hA_csrOffsets, (A_num_rows+1)*sizeof(int), cudaMemcpyHostToDevice));
    cudaErrCheck(cudaMemcpy(dA_columns, hA_columns, A_num_nnz*sizeof(int), cudaMemcpyHostToDevice));
#ifdef FLOAT
    cudaErrCheck(cudaMemcpy(dA_values, hA_values, A_num_nnz*sizeof(float),  cudaMemcpyHostToDevice));
    cudaErrCheck(cudaMemcpy(dX,      X         , A_num_cols*sizeof(float), cudaMemcpyHostToDevice));
#else
    cudaErrCheck(cudaMemcpy(dA_values, hA_values, A_num_nnz*sizeof(double),  cudaMemcpyHostToDevice));
    cudaErrCheck(cudaMemcpy(dX,      X,           A_num_cols*sizeof(double), cudaMemcpyHostToDevice));
#endif


    // CUSPARSE APIs
    cusparseHandle_t     handle = 0;
    cusparseSpMatDescr_t matA;
    cusparseDnVecDescr_t vecX, vecY;
    void*  dBuffer    = NULL;
    size_t bufferSize = 0;
    cusparseErrCheck(cusparseCreate(&handle));

    // Create sparse matrix A in CSR format
#ifdef FLOAT
    cusparseErrCheck( cusparseCreateCsr(&matA, A_num_rows, A_num_cols, A_num_nnz,
					dA_csrOffsets, dA_columns, dA_values,
					CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
					CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F));
    // Create dense vector X
    cusparseErrCheck( cusparseCreateDnVec(&vecX, A_num_cols, dX, CUDA_R_32F) );
    
    // Create dense vector y
    cusparseErrCheck( cusparseCreateDnVec(&vecY, A_num_rows, dY, CUDA_R_32F) );

    // allocate an external buffer if needed
    cusparseErrCheck( cusparseSpMV_bufferSize(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
					      &alpha, matA, vecX, &beta, vecY, CUDA_R_32F,
					      CUSPARSE_MV_ALG_DEFAULT, &bufferSize) );

#else
    cusparseErrCheck( cusparseCreateCsr(&matA, A_num_rows, A_num_cols, A_num_nnz,
					dA_csrOffsets, dA_columns, dA_values,
					CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
					CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F));
    // Create dense vector X
    cusparseErrCheck( cusparseCreateDnVec(&vecX, A_num_cols, dX, CUDA_R_64F) );
    
    // Create dense vector y
    cusparseErrCheck( cusparseCreateDnVec(&vecY, A_num_rows, dY, CUDA_R_64F) );
    
    // allocate an external buffer if needed
    cusparseErrCheck( cusparseSpMV_bufferSize(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
					      &alpha, matA, vecX, &beta, vecY, CUDA_R_64F,
					      CUSPARSE_MV_ALG_DEFAULT, &bufferSize) );
#endif

    cudaErrCheck( cudaMalloc(&dBuffer, bufferSize));
    
    
    // execute SpMV
    // Timing the solve 
    cudaEvent_t start;
    cudaEvent_t stop;
    cudaErrCheck(cudaEventCreate(&start));
    cudaErrCheck(cudaEventCreate(&stop));
    
    int REPEAT = 50;
    float times [REPEAT];
    for (int i = 0; i < REPEAT; i++) {
      cudaEventRecord(start);
#ifdef FLOAT
      cusparseErrCheck( cusparseSpMV(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
				     &alpha, matA, vecX, &beta, vecY, CUDA_R_32F,
				     CUSPARSE_MV_ALG_DEFAULT, dBuffer));
      
#else
      cusparseErrCheck( cusparseSpMV(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
				     &alpha, matA, vecX, &beta, vecY, CUDA_R_64F,
				     CUSPARSE_MV_ALG_DEFAULT, dBuffer) );
#endif	 
      
      cudaEventRecord(stop);  
      cudaEventSynchronize(stop);
      float milliseconds = 0;
      cudaEventElapsedTime(&milliseconds, start, stop);
      
      times[i] = milliseconds;
    }
    
    float sum_times = 0.0f;
    for (int i = 0; i < REPEAT; i++) {
      sum_times += times[i];
    }
    
#ifdef FLOAT
    printf ("SINGLE PRECISION SPMV ");
#else
    printf ("DOUBLE PRECISION SPMV ");
#endif
    printf ("solve time (microseconds) = %f\n", (sum_times/REPEAT)*1000);


    // device result check
#ifdef FLOAT
    cudaErrCheck( cudaMemcpy(Y, dY, A_num_rows * sizeof(float),
    cudaMemcpyDeviceToHost));
    
#else
    cudaErrCheck( cudaMemcpy(Y, dY, A_num_rows * sizeof(double),
    cudaMemcpyDeviceToHost) )
#endif


    // use mkl to check result
    sparse_matrix_t mklA;
    sparse_status_t stat;

#ifdef FLOAT
    stat = mkl_sparse_s_create_csr(&mklA, SPARSE_INDEX_BASE_ZERO, A_num_rows, A_num_cols,
				   hA_csrOffsets, hA_csrOffsets + 1,
				   hA_columns, hA_values);
#else
    stat = mkl_sparse_d_create_csr(&mklA, SPARSE_INDEX_BASE_ZERO, A_num_rows, A_num_cols,
				   hA_csrOffsets, hA_csrOffsets + 1,
				   hA_columns, hA_values);
#endif
    if (SPARSE_STATUS_SUCCESS != stat) {
      fprintf(stderr, "Failed to create mkl csr\n");
      return -1;
    }
#ifdef FLOAT    
    float *result = (float*)malloc(sizeof(float)*A_num_rows);
    float error = 0;
#else
    double *result = (double*)malloc(sizeof(double)*A_num_rows);
    double error = 0;
#endif
    
    matrix_descr descA;
    descA.type = SPARSE_MATRIX_TYPE_GENERAL;
    descA.diag = SPARSE_DIAG_NON_UNIT;

#ifdef FLOAT
    mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, mklA, descA, X, 0, result);
#else
    mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, mklA, descA, X, 0, result);
#endif

  #pragma omp parallel for
  for(int i = 0; i < A_num_rows; i++) {
    error += abs(result[i] - Y[i]);
  }
  printf ("Error= %e\n", error/A_num_cols);

    // destroy matrix/vector descriptors    
    cudaErrCheck(cudaEventDestroy(start));             
    cudaErrCheck(cudaEventDestroy(stop));
    free(hA_values);
    free(hA_csrOffsets);
    free(hA_columns);
    free(X);
    free(Y);
    free(result);
    cudaErrCheck(cudaFree(dY));
    cudaErrCheck(cudaFree(dX));
    cudaErrCheck(cudaFree(dA_values));
    cudaErrCheck(cudaFree(dA_csrOffsets));
    cudaErrCheck(cudaFree(dA_columns));
    cudaErrCheck(cudaFree(dBuffer));
    

    cusparseErrCheck( cusparseDestroySpMat(matA));
    cusparseErrCheck( cusparseDestroyDnVec(vecX));
    cusparseErrCheck( cusparseDestroyDnVec(vecY));
    cusparseErrCheck( cusparseDestroy(handle) );
    return 0;
}
