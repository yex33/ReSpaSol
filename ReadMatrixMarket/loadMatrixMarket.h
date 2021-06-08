#ifndef LOAD_MATRIX_MARKET_H
#define LOAD_MATRIX_MARKET_H

#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <malloc.h>
#include <x86intrin.h>
#include "mm_io.h"


#define MALLOC(type, len) (type *)_mm_malloc(sizeof(type)*(len), 64)
#define FREE(x) { if (x) _mm_free(x); x = NULL; }

typedef struct {
    int isSymmetric; 
	int m;
	int n;
    int nnz;
    int *rowptr;
    int *colidx;
    double *values;
} CSR;


typedef struct {
    int isSymmetric; 
	int m;
	int n;
    int nnz;
    int *Colidx;
    int *Rowidx;
    double *values;
} COO;




int loadMatrixMarket(const char *file, CSR *matrix, int outputBase, int transpose);
int loadCooMatrix(const char *file, COO *matrix, int outputBase, int transpose);


void qsort(int *idx, double *w, int left, int right);

#endif /* LOAD_MATRIX_MARKET_H */
