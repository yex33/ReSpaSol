#include <stdio.h>
#include "omp.h"
#include "mkl_types.h"
#include "mkl_spblas.h"
#include <string.h>
#include "mpi.h"
#include <math.h>
#include <float.h>

//#define FLOAT

#ifdef FLOAT
#include "smumps_c.h"
#else
#include "dmumps_c.h"
#endif

#define JOB_INIT -1
#define JOB_END -2


void set_ftz(void)
{
   asm("stmxcsr -0x4(%rsp)\n\t"          // store CSR register on stack
       "orl     $0x8040,-0x4(%rsp)\n\t"  // set bits 15 (flush-to-zero) and 7 (denormals-are-zero)
       "ldmxcsr -0x4(%rsp)");            // load CSR register from stack
}


#include "loadMatrixMarket.h"

int main(int argc, char ** argv)
{

    #ifdef FLOAT
    //   set_ftz();
    #endif
    if (argc < 2) {
        fprintf(
            stderr,
            "-- Usage examples --\n"
            "  %s inline_1.mtx type: run with inline_1 matrix in matrix market format\n",
            argv[0]);
        return -1;
    }
    
    #ifdef FLOAT
    SMUMPS_STRUC_C id;
    #else
    DMUMPS_STRUC_C id;
    #endif
    
    // Load sparse matrix
    COO matrixA;
    int outputbase = 0;
    // Load the matrix transpose to COO format
    loadCooMatrix(argv[1], &matrixA, outputbase, 0 /*transpose =false*/);

    // Set matrix parameters
    MUMPS_INT n   = matrixA.n;
    MUMPS_INT8 nnz = matrixA.nnz;

   // Convert value of A to single
    #ifdef FLOAT
    float *A = (float*)malloc(matrixA.nnz*sizeof(float));
    #pragma omp parallel for
    for(int i =0; i < matrixA.nnz; i++){
        A[i] = (float) matrixA.values[i];
        // if ((fpclassify(A[i]) == FP_SUBNORMAL ))
    }
    #else
    double *A = matrixA.values;
    #endif

    /* RHS and solution vectors. */
    #ifdef FLOAT
    float 
        *b = (float*)malloc(sizeof(float)*n),
        *bs = (float*)malloc(sizeof(float)*n),
        res;
    #else
    double 
        *b = (double*)malloc(sizeof(double)*n),
        *bs = (double*)malloc(sizeof(double)*n),
        res;
    #endif
    for(int i=0; i<n; i++) {
        b[i] = 1.0;
        bs[i] = 0.0;
    }

    MUMPS_INT myid;
    int error = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);


    /* Initialize a MUMPS instance. Use MPI_COMM_WORLD */
    // id.comm_fortran=USE_COMM_WORLD;
    id.comm_fortran=MPI_COMM_WORLD;
    id.par=1; id.sym=0;
    id.job=JOB_INIT;
    
    #ifdef FLOAT
    smumps_c(&id);
    #else
    dmumps_c(&id);
    #endif

    /* Define the problem on the host */
    if (myid == 0) {
        id.n = n; id.nnz =nnz; id.irn= matrixA.Rowidx; id.jcn=matrixA.Colidx;
        id.a = A; id.rhs = b;
    }
    #define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
    /* No outputs */
    id.ICNTL(1)=6; id.ICNTL(2)=1; id.ICNTL(3)=6; id.ICNTL(4)=2; id.ICNTL(11)=2;

    // Call the MUMPS package (analyse and factorization).
    id.job=4;
    double time_fact = omp_get_wtime();
    #ifdef FLOAT
    smumps_c(&id);
    #else
    dmumps_c(&id);
    #endif
    time_fact = omp_get_wtime() - time_fact;
    
    if (id.infog[0]<0) {
        printf(" (PROC %d) ERROR RETURN: \tINFOG(1)= %d\n\t\t\t\tINFOG(2)= %d\n",
               myid, id.infog[0], id.infog[1]);
        error = 1;
    }

    // Call the MUMPS package for solve.
    id.job=3;
    double time_solve = omp_get_wtime();
    #ifdef FLOAT
    smumps_c(&id);
    #else
    dmumps_c(&id);
    #endif
    time_solve = omp_get_wtime() - time_solve;
    if (id.infog[0]<0) {
        printf(" (PROC %d) ERROR RETURN: \tINFOG(1)= %d\n\t\t\t\tINFOG(2)= %d\n",
               myid, id.infog[0], id.infog[1]);
        error = 1;
    }

    if (myid == 0) {
        printf("Factor time = %lf\n Solve time = %lf\n", time_fact, time_solve);
    }
    // Terminate instance.
  id.job=JOB_END;
    #ifdef FLOAT
    smumps_c(&id);
    #else
    dmumps_c(&id);
    #endif
    MPI_Finalize();
    return 0;
}
