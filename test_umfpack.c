#include <omp.h>
#include<mkl.h>
#include <assert.h>
#include "amd.h"
#include "loadMatrixMarket.h"




#include "umfpack.h"

/* -------------------------------------------------------------------------- */
/* error: print a message and exit */
/* -------------------------------------------------------------------------- */

static void error (const char *message ) {
    printf ("\n\n====== error: %s =====\n\n", message) ;
    return;
}



int main(int argc, char *argv[])
{

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

    // transpose A in CSR format 
    loadMatrixMarket(argv[1], &matrixA, outputbase, 1 /*transpose =1*/);

    int *Ap = matrixA.rowptr;
    int *Ai = matrixA.colidx;
    double *Ax = matrixA.values;

 
    // Set matrix parameters
    int n   = matrixA.n;
    int m   = matrixA.m;
    int nnz = matrixA.nnz;
    double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO];
    void *Symbolic, *Numeric;
    int status;
    int sys = UMFPACK_A; // solve Ax = b

    // declare and set rhs
    int seed[] = {0, 0, 0, 1};
    lapack_int retval;
    double *B = (double*) malloc(m*sizeof(double));
    retval = LAPACKE_dlarnv(1, seed, m, B);
    double * X = (double*) malloc(m*sizeof(double));
    assert(retval == 0);
    
    
    printf("Dimension %dx%d; # nonzeros %d\n", m, n, nnz);
    printf ("\nUMFPACK V%d.%d (%s) demo: _di_ version\n",
            UMFPACK_MAIN_VERSION, UMFPACK_SUB_VERSION, UMFPACK_DATE) ;

    /* change the default print level for this demo */
    /* (otherwise, nothing will print) */

    /* get the default control parameters */
    umfpack_di_defaults (Control);
    Control [UMFPACK_PRL] = 2;
    Control [UMFPACK_PIVOT_TOLERANCE] = 0.1;
    Control [UMFPACK_ORDERING] = UMFPACK_ORDERING_CHOLMOD;
    //Control [UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
    /* print the control parameters */
    umfpack_di_report_control (Control);


    double t1 = omp_get_wtime();

    status = umfpack_di_symbolic (m, n, Ap, Ai, Ax, &Symbolic, Control, Info);
    if (status < 0)
    {
        error ("umfpack_di_symbolic failed");
        umfpack_di_report_info (Control, Info) ;
        umfpack_di_report_status (Control, status);
        return -1;
    }

    status = umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, Control, Info);
    if (status < 0)
    {
        error ("umfpack_di_numeric failed");
        umfpack_di_report_info (Control, Info) ;
        umfpack_di_report_status (Control, status);
        return -1;
    } 
    
    status = umfpack_di_solve (sys, Ap, Ai, Ax, X, B, Numeric, Control, Info);
    if (status < 0)
    {
        error ("umfpack_di_solve failed");
        umfpack_di_report_info (Control, Info) ;
        umfpack_di_report_status (Control, status);
        return -1;
    }
    double t2 = omp_get_wtime();
    umfpack_di_report_info (Control, Info);
    
    double res0 = 0.0;
    #pragma omp for
    for (int j = 0; j < n; j++)
    {
        res0 += B[j] * B[j];
    }

    // B = B - A*x
    int p, rowidx;
    double temp;
    #pragma omp for private(p, rowidx, temp)
    for(int j=0; j < m; j++) {
        for( p=Ap[j]; p <Ap[j+1]; p++) {
            rowidx= Ai[p];
            temp = X[j] * Ax[p];
            #pragma omp atomic
            B[rowidx] -= temp;
        }
    }

    printf("end B = B - AX\n");
    double res = 0.0;
    #pragma omp for
    for (int j = 0; j < n; j++)
    {
        res += B[j] * B[j];
    }

    res = sqrt(res) / sqrt(res0);
    printf("\nRelative residual = %e\n", res);
    // Check residual
    if (res > 1e-10)
    {
        printf("Error: residual is too high!\n");
        // exit(-1);
    }
    printf("TIMING: total  = %gs\n", t2 - t1);
    
    umfpack_di_free_symbolic (&Symbolic);
    umfpack_di_free_numeric (&Numeric);
    free(Ai);
    free(Ap);
    free(Ax);
    free(X);
    free(B);

}

