#include "../loadMatrixMarket.h"

int main(int argc, char **argv)
{

    CSR matrix;

    if (argc < 2) {
        fprintf(
            stderr,
            "-- Usage examples --\n"
            "  %s inline_1.mtx: run with inline_1 matrix in matrix market format\n\n",
            argv[0]);
        return -1;
    }

    int outputbase = 0;
    int transpose = 0;
    loadMatrixMarket(argv[1], &matrix, outputbase, transpose);

    // Let's print some values
    int m = matrix.m;
    int n = matrix.n;
    int nnz = matrix.nnz;

    printf("m=%d\t n=%d\t nnz=%d\n", m, n, nnz);

    printf("Printing values\n");

    for (int i = 0; i < nnz; i++ )
        printf("%.9lf\n", matrix.values[i]);
    return 0;
}
