#include<unistd.h>
#include "omp.h"
#include <stdlib.h>

#define FLOAT
#ifdef FLOAT // set to 0 for single precision 
extern "C" {
#include "slu_mt_ddefs.h"
}
#else
extern "C" {
    #include "slu_mt_sdefs.h"
}
#endif

#include "loadMatrixMarket.h"


void set_ftz(void)
{
   asm("stmxcsr -0x4(%rsp)\n\t"          // store CSR register on stack
       "orl     $0x8040,-0x4(%rsp)\n\t"  // set bits 15 (flush-to-zero) and 7 (denormals-are-zero)
       "ldmxcsr -0x4(%rsp)");            // load CSR register from stack
}


#ifdef FLOAT
void sCreate_CompRow_Matrix(SuperMatrix *A, int_t m, int_t n, int_t nnz, double *nzval,
		      int_t *colind, int_t *rowptr,
		      Stype_t stype, Dtype_t dtype, Mtype_t mtype)
#else
void dCreate_CompRow_Matrix(SuperMatrix *A, int_t m, int_t n, int_t nnz, double *nzval,
		      int_t *colind, int_t *rowptr,
		      Stype_t stype, Dtype_t dtype, Mtype_t mtype)
#endif
{
    NRformat *Astore;

    A->Stype = stype;
    A->Dtype = dtype;
    A->Mtype = mtype;
    A->nrow = m;
    A->ncol = n;
    A->Store = (void *) SUPERLU_MALLOC( sizeof(NRformat) );
    Astore = (NRformat *) A->Store;
    Astore->nnz = nnz;
    Astore->nzval = nzval;
    Astore->colind = colind;
    Astore->rowptr = rowptr;
}


int main(int argc, char *argv[])
{

    #ifdef FLOAT
    // set_ftz();
    #endif

    int      *perm_r; /* row permutations from partial pivoting */
    int      *perm_c; /* column permutation vector */
    SuperMatrix   L;       /* factor L */
    SuperMatrix   U;       /* factor U */
    SuperMatrix   B;
    int      permc_spec;


    int    n, m, nnz;
    int    nprocs            = 1;

    if (argc < 2) {
        fprintf(
            stderr,
            "-- Usage examples --\n"
            "  %s inline_1.mtx type: run with inline_1 matrix in matrix market format\n",
            argv[0]);
        return -1;
    }

    // Load sparse matrix
    CSR matrixA;
    int outputbase = 0;
    // Load the matrix transpose to get CSC format even though using CSR
    loadMatrixMarket(argv[1], &matrixA, outputbase, 1 /*transpose =true*/);

    // Set matrix parameters
    n   = matrixA.n;
    m   = matrixA.m;
    nnz = matrixA.nnz;

    SuperMatrix A;
    // Convert value of A to single
    #ifdef FLOAT
        float *Afp32 = (float*)malloc(matrixA.nnz*sizeof(float));
        #pragma omp parallel for
        for(int i =0; i < matrixA.nnz; i++){
            Afp32[i] = (float) matrixA.values[i];
        }
    #endif

    #ifdef FLOAT
        sCreate_CompCol_Matrix(&A, n, m, nnz,
        Afp32, matrixA.colidx, matrixA.rowptr,
        SLU_NC, SLU_S, SLU_GE);
    #else
         dCreate_CompCol_Matrix(&A, n, m, nnz,
        matrixA.values, matrixA.colidx, matrixA.rowptr,
        SLU_NC, SLU_D, SLU_GE);
    #endif

    NCformat *Astore = (NCformat *)A.Store;
    printf("Dimension " IFMT "x" IFMT "; # nonzeros " IFMT "\n", A.nrow, A.ncol, Astore->nnz);

    // Generate rhs
    int nrhs   = 1;

    #ifdef FLOAT
        float *rhs = (float*) malloc(m*sizeof(float));
        sCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_S, SLU_GE);
        float * xact = (float*) malloc((n * nrhs)*sizeof(float));
        int ldx = n;
        sGenXtrue(n, nrhs, xact, ldx);
        sFillRHS(NOTRANS, nrhs, xact, ldx, &A, &B);
    #else
        double *rhs = (double*) malloc(m*sizeof(double));
        dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);
        double * xact = doubleMalloc(n * nrhs);
        int ldx = n;
        dGenXtrue(n, nrhs, xact, ldx);
        dFillRHS(NOTRANS, nrhs, xact, ldx, &A, &B);
    #endif
    // Make a copy of b
     DNformat    *Bstore;
    #ifdef FLOAT
        float *b = (float*) malloc (sizeof(float)*n);
        Bstore = (DNformat *) B.Store;
        float  *bvalues = (float*) Bstore->nzval;
    #else
        double *b = (double*) malloc (sizeof(double)*n);
        Bstore = (DNformat *) B.Store;
        double  *bvalues = (double*) Bstore->nzval;
    #endif
    
    for(int i = 0; i< n; i++) {
        b[i] = bvalues[i];
    }

    //if (!(perm_r = intMalloc(m))) SUPERLU_ABORT("Malloc fails for perm_r[].");
    perm_r = (int*) malloc(m*sizeof(int));
    if (!(perm_c = intMalloc(n))) SUPERLU_ABORT("Malloc fails for perm_c[].");

    /*
     * Get column permutation vector perm_c[], according to permc_spec:
     *   permc_spec = 0: natural ordering 
     *   permc_spec = 1: minimum degree ordering on structure of A'*A
     *   permc_spec = 2: minimum degree ordering on structure of A'+A
     *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
     */
    
    double t1 = omp_get_wtime();
    permc_spec = 3;
    get_perm_c(permc_spec, &A, perm_c);

    double time_analysis = omp_get_wtime() - t1;
    int info;
    
    #ifdef FLOAT
        psgssv(nprocs, &A, perm_c, perm_r, &L, &U, &B, &info);
     #else
        pdgssv(nprocs, &A, perm_c, perm_r, &L, &U, &B, &info);
     #endif
     double t2 = omp_get_wtime();

    

    /* -------------------------------------------------------------------- */
    /* .. Compute residual */
    /* -------------------------------------------------------------------- */
    // get the computed solution
    DNformat    *Xstore;
    Xstore = (DNformat *) B.Store;

    #ifdef FLOAT
        float * Xvalues = (float*) Xstore->nzval;
        float *x = (float*) malloc (sizeof(float)*n);
    #else
        double * Xvalues = (double*) Xstore->nzval;
        double *x = (double*) malloc (sizeof(double)*n);
    #endif
    for (int i=0; i<n; i++){
        x[i] = Xvalues[i];
    }

    // Compute b = A*X
    #ifdef FLOAT
        sFillRHS(NOTRANS, nrhs, x, ldx, &A, &B);
    #else
        dFillRHS(NOTRANS, nrhs, x, ldx, &A, &B);
    #endif
    DNformat    *bstore;
    bstore = (DNformat *) B.Store;
    #ifdef FLOAT
        float *bs = (float*) malloc (sizeof(float)*n);
        bvalues = (float*) bstore->nzval;
    #else
        double *bs = (double*) malloc (sizeof(double)*n);
        bvalues = (double*) bstore->nzval;
    #endif

    for (int i=0; i<n; i++){
        bs[i] = bvalues[i];
    }

    #ifdef FLOAT
        float res = 0.0;
        float res0 = 0.0;
    #else
        double res = 0.0;
        double res0 = 0.0;
    #endif
    for (int j = 0; j < n; j++)
    {
        res += (bs[j] - b[j]) * (bs[j] - b[j]);
        res0 += b[j] * b[j];
    }
    res = sqrt(res) / sqrt(res0);
    printf("\nRelative residual = %e\n", res);
    // Check residual
    if (res > 1e-10)
    {
        printf("Error: residual is too high!\n");
        // exit(-1);
    }

    
    /* ***********************************************************************/

    free(b);
    free(bs);
    SUPERLU_FREE (rhs);
    SUPERLU_FREE (xact);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_SCP(&L);
    Destroy_CompCol_NCP(&U);

    printf("\n Reordering TIMING:  total  = %gs\n", time_analysis);
    printf("\n  Total TIMING   = %gs\n", t2 - t1);

}

