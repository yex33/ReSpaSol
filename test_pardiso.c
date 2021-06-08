#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fenv.h>
#include <algorithm>

#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl_spblas.h"
#include "omp.h"

#include "loadMatrixMarket.h"

using namespace std;
#define FLOAT

void set_ftz(void)
{
   asm("stmxcsr -0x4(%rsp)\n\t"          // store CSR register on stack
       "orl     $0x8040,-0x4(%rsp)\n\t"  // set bits 15 (flush-to-zero) and 7 (denormals-are-zero)
       "ldmxcsr -0x4(%rsp)");            // load CSR register from stack
}



static const size_t LLC_CAPACITY = 32*1024*1024;
static const double *bufToFlushLlc = NULL;

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


MKL_INT main(int argc, char *argv[])
{
    #ifdef FLOAT
    // set_ftz();
    #endif

    if (argc < 2) {
        fprintf(
            stderr,
            "-- Usage examples --\n"
            "  %s inline_1.mtx type: run with inline_1 matrix in matrix market format\n",
            argv[0]);
        return -1;
    }

    // allocate a large buffer to flush out cache
    int nbFlush = 2048;
    bufToFlushLlc = (double *)_mm_malloc(LLC_CAPACITY, 64);
    
   CSR matrixA;
   int outputbase = 1;
   loadMatrixMarket(argv[1], &matrixA, outputbase, 0 /*transpose =false*/);


    MKL_INT    n  = matrixA.n;
    MKL_INT   *ia =  matrixA.rowptr;
    MKL_INT   *ja =  matrixA.colidx; 

    #ifdef FLOAT
    float *a = (float*)malloc(matrixA.nnz*sizeof(float));
    for(int i =0; i < matrixA.nnz; i++){
        a[i] = (float) matrixA.values[i];   
    }
    #else
    double    *a  = matrixA.values;
    #endif
    

    MKL_INT mtype = 11; /* Real unsymmetric matrix */
    // Descriptor of main sparse matrix properties
    struct matrix_descr descrA;
    // Structure with sparse matrix stored in CSR format
    sparse_matrix_t csrA;
    /* RHS and solution vectors. */
    #ifdef FLOAT
    float 
       *b = (float*)malloc(sizeof(float)*n),
       *x = (float*)malloc(sizeof(float)*n),
       *bs = (float*)malloc(sizeof(float)*n),
       res, res0;    
    #else
    double 
       *b = (double*)malloc(sizeof(double)*n),
       *x = (double*)malloc(sizeof(double)*n),
       *bs = (double*)malloc(sizeof(double)*n),
       res, res0;
    #endif
    for(int i=0; i<n; i++) {
        b[i] = 1;
        bs[i] = 0;
        x[i] = 0;
    }


    MKL_INT *perm = (MKL_INT*)malloc(n*sizeof(MKL_INT));
    for(int i=0; i<n; i++) perm[i] = 0;

    MKL_INT nrhs = 1; /* Number of right hand sides. */
    /* Internal solver memory pointer pt, */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
    /* or void *pt[64] should be OK on both architectures */
    void *pt[64];
    /* Pardiso control parameters. */
    MKL_INT iparm[64];
    MKL_INT maxfct, mnum, phase, error, msglvl;
    /* Auxiliary variables. */
    MKL_INT i, j;
    
    double t1, t2, t3, t4;
    
    /* -------------------------------------------------------------------- */
    /* .. Setup Pardiso control parameters. */
    /* -------------------------------------------------------------------- */
    for (i = 0; i < 64; i++)
    {
        iparm[i] = 0;
    }
    iparm[0] = 1;   /* No solver default */
    iparm[1] = 3;   /* 2 Fill-in reordering from METIS , use 3 for parallel analyze */ 
    iparm[3] = 0;   /* No iterative-direct algorithm */
    iparm[4] = 0;   /* No user fill-in reducing permutation */
    iparm[5] = 0;   /* Write solution into x */
    iparm[6] = 0;   /* Not in use */
    iparm[7] = 0;   /* Max numbers of iterative refinement steps */
    iparm[8] = 0;   /* Not in use */
    #ifdef FLOAT
    iparm[9] = 4;  /* Perturb the pivot elements with 1E-13 */
    #else
    iparm[9] = 13;  /* Perturb the pivot elements with 1E-13 */
    #endif
    iparm[10] = 1;  /* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 0;  /* Conjugate/transpose solve */
    iparm[12] = 1;  /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
    iparm[13] = 0;  /* Output: Number of perturbed pivots */
    iparm[14] = 0;  /* Not in use */
    iparm[15] = 0;  /* Not in use */
    iparm[16] = 0;  /* Not in use */
    iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1; /* Output: Mflops for LU factorization */
    iparm[19] = 0;  /* Output: Numbers of CG Iterations */
    iparm[23] = 1;  /* Parallel factorization control  */
    #ifdef FLOAT
    iparm[27] = 1;  // Set to single precision.
    #else
    iparm[27] = 0;  // Set to double precision.
    #endif

    maxfct = 1;     /* Maximum number of numerical factorizations. */
    mnum = 1;       /* Which factorization to use. */
    msglvl =0;     /* Print statistical information  */
    error = 0;      /* Initialize error flag */
                    /* -------------------------------------------------------------------- */
                    /* .. Initialize the internal solver memory pointer. This is only */
                    /* necessary for the FIRST call of the PARDISO solver. */
                    /* -------------------------------------------------------------------- */
    for (i = 0; i < 64; i++)
    {
        pt[i] = 0;
    }

    t1 = omp_get_wtime();

    /* -------------------------------------------------------------------- */
    /* .. Reordering and Symbolic Factorization. This step also allocates */
    /* all memory that is necessary for the factorization. */
    /* -------------------------------------------------------------------- */
    phase = 11;    
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
            &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, b, x, &error);

    t2 = omp_get_wtime();

     printf("TIMING: symbolic factorization = %gs\n", t2 - t1);

    if (error != 0)
    {
        printf("\nERROR during symbolic factorization: %d", error);
        exit(1);
    }
    //printf("\nReordering completed ... ");
    //printf("\nNumber of nonzeros in factors = %d", iparm[17]);
    //printf("\nNumber of factorization MFLOPS = %d\n", iparm[18]);
    /* -------------------------------------------------------------------- */
    /* .. Numerical factorization. */
    /* -------------------------------------------------------------------- */
 
    phase = 22;
    for (int j = 0; j < nbFlush; ++j) flushLlc();
    t3 = omp_get_wtime();
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
            &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, b, x, &error);
    double time_fact =  omp_get_wtime() -t3;

    printf("TIMING: numerical factorization = %gs\n", time_fact);
    if (error != 0)
    {
        printf("\nERROR during numerical factorization: %d", error);
        exit(2);
    }
    //printf("\nFactorization completed ... ");

    /* -------------------------------------------------------------------- */
    /* .. Solution phase. */
    /* -------------------------------------------------------------------- */
    phase = 33;

    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
     descrA.mode = SPARSE_FILL_MODE_UPPER;
    descrA.diag = SPARSE_DIAG_NON_UNIT;
    #ifdef FLOAT
    mkl_sparse_s_create_csr(&csrA, SPARSE_INDEX_BASE_ONE, n, n, ia, ia + 1, ja, a);
    #else
    mkl_sparse_d_create_csr(&csrA, SPARSE_INDEX_BASE_ONE, n, n, ia, ia + 1, ja, a);
    #endif


    /* Set right hand side to one. */
    for (i = 0; i < n; i++)
    {
        b[i] = 1;
    }

    
    t3 = omp_get_wtime();
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
	    &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, b, x, &error);
    t4 = omp_get_wtime();
    
    printf("Solve : error = %d\n", error);
    // printf("TIMING: solving the system = %gs\n", t4 - t3);
    printf("TIMING: solving the system = %gs\n", t4 -t3);

    if (error != 0)
    {
        printf("\nERROR during solution: %d", error);
        exit(3);
    }

    printf("TIMING: total FACTORIZATION time = %gs\n", t2 - t1 + time_fact);
    
    printf("TIMING: total time = %gs\n", t2 - t1 + time_fact + t4 - t3);


    // Compute residual
    // the CSC format for A is the CSR format for A transposed
    #ifdef FLOAT
    mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, csrA, descrA, x, 0.0, bs);
    #else
    mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, csrA, descrA, x, 0.0, bs);
    #endif
    res = 0.0;
    res0 = 0.0;
    for (j = 1; j <= n; j++)
    {
        res += (bs[j - 1] - b[j - 1]) * (bs[j - 1] - b[j - 1]);
        res0 += b[j - 1] * b[j - 1];
    }
    res = sqrt(res) / sqrt(res0);
    printf("\n Relative residual = %e\n", res);
    // Check residual
    if (res > 1e-10)
    {
        printf("Error: residual is too high!\n");
        exit(10 + i);
    }
    mkl_sparse_destroy(csrA);

    /* -------------------------------------------------------------------- */
    /* .. Termination and release of memory. */
    /* -------------------------------------------------------------------- */
    phase = -1; /* Release internal memory. */
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
            &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, b, x, &error);

    return 0;
}
