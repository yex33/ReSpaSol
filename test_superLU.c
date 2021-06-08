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

    // set_ftz();
    int      *perm_c; /* column permutation vector */
    int      *perm_r; /* row permutations from partial pivoting */
    SuperMatrix L;      /* factor L */
    SCformat *Lstore;
    SuperMatrix U;      /* factor U */
    NCformat *Ustore;
    SuperMatrix B;
    mem_usage_t   mem_usage;
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

    // Load sparse matrix
    CSR matrixA;
    int outputbase = 0;
    set_default_options(&options);
    loadMatrixMarket(argv[1], &matrixA, outputbase, 0 /* transpose =false*/);

    // Set matrix parameters
    int n   = matrixA.n;
    int m   = matrixA.m;
    int nnz = matrixA.nnz;

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
        sCreate_CompRow_Matrix(&A, m, n, nnz,
                           Afp32, matrixA.colidx, matrixA.rowptr,
                           SLU_NC, SLU_S, SLU_GE);
    #else
        dCreate_CompRow_Matrix(&A, m, n, nnz,
                           matrixA.values, matrixA.colidx, matrixA.rowptr,
                           SLU_NC, SLU_D, SLU_GE);
    #endif
    NCformat *Astore = (NCformat *)A.Store;
    printf("Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, Astore->nnz);

    // Number of column in B
    int nrhs   = 1;

    #ifdef FLOAT
        float *rhs = (float*) malloc(m*sizeof(float));
        sCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_S, SLU_GE);
        float * xact = (float*) malloc((n * nrhs)*sizeof(float));
        int ldx = n;
        sGenXtrue(n, nrhs, xact, ldx);
        sFillRHS(options.Trans, nrhs, xact, ldx, &A, &B);
    #else
        double *rhs = (double*) malloc(m*sizeof(double));
        dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);
        double * xact = doubleMalloc(n * nrhs);
        int ldx = n;
        dGenXtrue(n, nrhs, xact, ldx);
        dFillRHS(options.Trans, nrhs, xact, ldx, &A, &B);
    #endif
    
    if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
    if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");

    /* Initialize the statistics variables. */
    StatInit(&stat);

    /*
     * Get column permutation vector perm_c[], according to permc_spec:
     *   permc_spec = 0: natural ordering 
     *   permc_spec = 1: minimum degree ordering on structure of A'*A
     *   permc_spec = 2: minimum degree ordering on structure of A'+A
     *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
     */    	
    int permc_spec = 3;
    get_perm_c(permc_spec, &A, perm_c);

    printf("starting factorization\n");
    int info;

    #ifdef FLOAT
        sgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
    #else
        dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
    #endif
    if ( info == 0 ) {

        /* This is how you could access the solution matrix. */
        #ifdef FLOAT
            float *sol = (float*) ((DNformat*) B.Store)->nzval; 
        #else
            double *sol = (double*) ((DNformat*) B.Store)->nzval; 
        #endif
            
            /* Compute the infinity norm of the error. */
        #ifdef FLOAT
            sinf_norm_error(nrhs, &B, xact);
        #else
            dinf_norm_error(nrhs, &B, xact);
        #endif 
            
        Lstore = (SCformat *) L.Store;
        /* Ustore = (NCformat *) U.Store; */
    	/* printf("#NZ  in factor L = %d\n", Lstore->nnz); */
    	/* printf("#NZ  in factor U = %d\n", Ustore->nnz); */
    	/* printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - n); */
    	/* printf("FILL ratio = %.1f\n", (float)(Lstore->nnz + Ustore->nnz - n)/nnz); */

        #ifdef FLOAT
            sQuerySpace(&L, &U, &mem_usage);
        #else
            dQuerySpace(&L, &U, &mem_usage);
        #endif
        /* printf("L\\U MB %.3f\ttotal MB needed %.3f\n", */
        /*        mem_usage.for_lu/1e6, mem_usage.total_n eeded/1e6); */
	
    } else {
        printf("dgssv() error returns INFO= %d\n", info);
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
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);

    /* #ifdef FLOAT */
    /* free(Afp32); */
    /* #endif */
    
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC("Exit main()");
#endif
}

