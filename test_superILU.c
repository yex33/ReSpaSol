//#define FLOAT

#ifdef FLOAT
#include "slu_sdefs.h"
#else
#include "slu_ddefs.h"
#endif

#include "loadMatrixMarket.h"


void set_ftz(void)
{
   asm("stmxcsr -0x4(%rsp)\n\t"          // store CSR register on stack
       "orl     $0x8040,-0x4(%rsp)\n\t"  // set bits 15 (flush-to-zero) and 7 (denormals-are-zero)
       "ldmxcsr -0x4(%rsp)");            // load CSR register from stack
}

int main(int argc, char *argv[])
{

    //set_ftz();
    char           equed[1];
    SuperMatrix    A, L, U;
    SuperMatrix    B, X;
    NCformat       *Astore;
    GlobalLU_t	   Glu; /* facilitate multiple factorizations with 
                           SamePattern_SameRowPerm                  */
    int            *perm_r; /* row permutations from partial pivoting */
    int            *perm_c; /* column permutation vector */
    int            *etree;
    void           *work;
    int            info, lwork, nrhs;
    mem_usage_t    mem_usage;
    superlu_options_t options;
    SuperLUStat_t stat;

    if (argc < 2) {
        fprintf(
            stderr,
            "-- Usage examples --\n"
            "  %s inline_1.mtx type: run with inline_1 matrix in matrix market format\n",
            argv[0]);
        return -1;
    }
    /* Defaults */
    lwork = 0;

    // Load sparse matrix
    CSR matrixA;
    int outputbase = 0;

    ilu_set_default_options(&options);
    loadMatrixMarket(argv[1], &matrixA, outputbase, 0 /* transpose =false*/);

    // Set matrix parameters
    int n   = matrixA.n;
    int m   = matrixA.m;
    int nnz = matrixA.nnz;

    // Convert value of A to single
    #ifdef FLOAT
        float *Afp32 = (float*)malloc(matrixA.nnz*sizeof(float));
        #pragma omp parallel for
        for(int i =0; i < matrixA.nnz; i++){
            Afp32[i] = (float) matrixA.values[i];
        }
    #endif
    
    #ifdef FLOAT
        sCreate_CompRow_Matrix(&A, m, n, nnz,
                           Afp32, matrixA.colidx, matrixA.rowptr,
                           SLU_NC, SLU_S, SLU_GE);
    #else
        dCreate_CompRow_Matrix(&A, m, n, nnz,
                           matrixA.values, matrixA.colidx, matrixA.rowptr,
                           SLU_NC, SLU_D, SLU_GE);
    #endif
    Astore = (NCformat *)A.Store;
    printf("Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, Astore->nnz);

    // Number of column in B
    nrhs   = 1;
    int ldx = n;
    #ifdef FLOAT
        float          rpg, rcond;
        float *R = (float*) malloc(m*sizeof(float));
        float *C = (float*) malloc(m*sizeof(float));
        float *rhs = (float*) malloc(m*sizeof(float));
        sCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_S, SLU_GE);
        float * xact = (float*) malloc((n * nrhs)*sizeof(float));
        sCreate_Dense_Matrix(&X, m, nrhs, xact, m, SLU_DN, SLU_S, SLU_GE);

        sGenXtrue(n, nrhs, xact, ldx);
        sFillRHS(options.Trans, nrhs, xact, ldx, &A, &B);
    #else
        double          rpg, rcond;
        double *R = (double*) malloc(m*sizeof(double));
        double *C = (double*) malloc(m*sizeof(double));
        double *rhs = (double*) malloc(m*sizeof(double));
        dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);
        double * xact = doubleMalloc(n * nrhs);
        dCreate_Dense_Matrix(&X, m, nrhs, xact, m, SLU_DN, SLU_D, SLU_GE);
        dGenXtrue(n, nrhs, xact, ldx);
        dFillRHS(options.Trans, nrhs, xact, ldx, &A, &B);
    #endif

    if ( !(etree = intMalloc(n)) ) ABORT("Malloc fails for etree[].");
    if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
    if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");

    /* Initialize the statistics variables. */
    StatInit(&stat);

    printf("starting factorization\n");

    #ifdef FLOAT
        sgsisx(&options, &A, perm_c, perm_r, etree, equed, R, C,
           &L, &U, work, lwork, &B, &X, &rpg, &rcond,
           &Glu, &mem_usage, &stat, &info);
    #else
        dgsisx(&options, &A, perm_c, perm_r, etree, equed, R, C,
           &L, &U, work, lwork, &B, &X, &rpg, &rcond,
           &Glu, &mem_usage, &stat, &info);
    #endif
    if ( info == 0 ) {

        /* Compute the infinity norm of the error. */
        #ifdef FLOAT
            sinf_norm_error(nrhs, &B, xact);
        #else
            dinf_norm_error(nrhs, &B, xact);
        #endif 

        #ifdef FLOAT
            sQuerySpace(&L, &U, &mem_usage);
        #else
            dQuerySpace(&L, &U, &mem_usage);
        #endif
        /* printf("L\\U MB %.3f\ttotal MB needed %.3f\n", */
        /*        mem_usage.for_lu/1e6, mem_usage.total_n eeded/1e6); */
	
    } else {
        printf("dgsisx() error returns INFO= %d\n", info);
        if ( info <= n ) { /* factorization completes */
            #ifdef FLOAT
                sQuerySpace(&L, &U, &mem_usage);
            #else
                dQuerySpace(&L, &U, &mem_usage);
            #endif
            printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
            mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
        }
    }

    if ( options.PrintStat ) StatPrint(&stat);
    StatFree(&stat);

    SUPERLU_FREE (rhs);
    SUPERLU_FREE (xact);
    SUPERLU_FREE (etree);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    SUPERLU_FREE (R);
    SUPERLU_FREE (C);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperMatrix_Store(&X);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC("Exit main()");
#endif
}

