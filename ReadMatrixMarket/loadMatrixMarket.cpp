#include "loadMatrixMarket.h"

using namespace std; 

void qsort(int *idx, double *w, int left, int right)
{
    if (left >= right) return;
    
    std::swap(idx[left], idx[left + (right - left)/2]);
    std::swap(w[left], w[left + (right - left)/2]);
    
    int last = left;
    for (int i = left+1; i <= right; i++) {
        if (idx[i] < idx[left]) {
            ++last;
            std::swap(idx[last], idx[i]);
            std::swap(w[last], w[i]);
        }
    }
    
    std::swap(idx[left], idx[last]);
    std::swap(w[left], w[last]);
    
    qsort(idx, w, left, last-1);
    qsort(idx, w, last+1, right);
}

/********************************************************************//**
 * loadMatrixMarket load matrix in a matrix market format
 * 
 * @param[in] file
 *            The input file to read the matrix from
 * 
 * @param[out] matrix
 *           The data structure that contains the sparse  matrix
 *           in either CSR or CSC format 
 *
 * @param[in] outputbase
 *          - 0: the indexing starts from 0 (for fortran)
 *          - 1: the indexing starts from 1 (for C/C++)
 *
 * @param[in] transpose
 *           - 0: the matrix is in CSR format
 *           - 1: the matrix is in CSC format
 *
 ***********************************************************************/
int loadMatrixMarket(const char *file, CSR *matrix, int outputbase, int transpose)
{
  FILE *fp=fopen(file, "r");
  if (NULL == fp) {
    fprintf(stderr, "Failed to open file %s\n", file);
    exit(-1);
  }

  // read banner
  MM_typecode matcode;
  if (mm_read_banner (fp, &matcode) != 0) {
    fprintf(stderr, "Error: could not process Matrix Market banner.\n");
    fclose(fp);
    return false;
  }

  if (!mm_is_valid (matcode) || mm_is_array (matcode) || mm_is_dense (matcode) ) {
    fprintf(stderr, "Error: only support sparse and real matrices.\n");
    fclose(fp);
    return false;
  }
  bool pattern = mm_is_pattern (matcode);

  // read sizes
  int m, n;
  int nnz; // # of non-zeros specified in the file
  if (mm_read_mtx_crd_size(fp, &m, &n, &nnz) !=0) {
    fprintf(stderr, "Error: could not read matrix size.\n");
    fclose(fp);
    return false;
  }

  if (transpose) {
      swap(m, n);
  }

  size_t count;
  if (mm_is_symmetric (matcode) == 1) {
    matrix->isSymmetric = 1;
    count = 2L*nnz;
  }
  else {
      matrix->isSymmetric = 0;
      count = nnz;
  }

  
  // allocate memory
  double *cooValues = MALLOC(double, count);
  int *cooColidx = MALLOC(int, count);
  int *cooRowidx = MALLOC(int, count);
  if (!cooValues || !cooColidx || !cooRowidx) {
    fprintf(stderr, "Failed to allocate memory\n");
    fclose(fp);
    return false;
  }

  int *colidx_temp=NULL, *rowcnt = NULL;

  if (matrix->isSymmetric) {
      colidx_temp = MALLOC(int, count);
      rowcnt = MALLOC(int, m + 1);
      if (!colidx_temp || !rowcnt) {
          fprintf(stderr, "Failed to allocate memory\n");
          fclose(fp);
          return false;
      }
      memset(rowcnt, 0, sizeof(int)*(m + 1));
  }

  // read values
  count = 0;
  int lines = 0;
  int x, y;
  double real, imag;
  int base = 1;
  while (mm_read_mtx_crd_entry (fp, &x, &y, &real, &imag, matcode) == 0) {
      if (transpose) swap(x, y);

    if (x > m || y > n) {
      fprintf(stderr, "Error: (%d %d) coordinate is out of range.\n", x, y);
      fclose(fp);
      return false;
    }

    cooRowidx[count] = x;
    cooColidx[count] = y;
    cooValues[count] = pattern ? 1 : real;
    if (0 == x || 0 == y) base = 0;

    ++count;
    ++lines;
    if (matrix->isSymmetric) rowcnt[x]++;
      // this is not a bug. we're intentionally indexing rowcnt[x] instead of rowcnt[x-1]
  }
  fclose(fp);
  
  if (0 == base) {
    for (size_t i = 0; i < count; ++i) {
      cooRowidx[i]++;
      cooColidx[i]++;
    }
    if (matrix->isSymmetric) {
        for (int i = m; i > 0; --i) {
            rowcnt[i] = rowcnt[i - 1];
        }
    } 
  }

  if (lines != nnz) {
    fprintf(stderr, "Error: nnz (%d) specified in the header doesn't match with # of lines (%d) in file %s\n",
            nnz, lines, file);
    return false;
  }
  
  if (matrix->isSymmetric) {
    // add transposed elements only if it doesn't exist
    size_t real_count = count;
    // preix-sum
    for (int i = 0; i < m; ++i) {
      rowcnt[i + 1] += rowcnt[i];
    }
    for (size_t i = 0; i < count; ++i) {
      int j = rowcnt[cooRowidx[i] - 1];
      colidx_temp[j] = cooColidx[i];
      rowcnt[cooRowidx[i] - 1]++;
    }
    for (int i = m; i > 0; --i) {
      rowcnt[i] = rowcnt[i - 1];
    }
    rowcnt[0] = 0;

#pragma omp parallel for
    for (int i = 0; i < m; ++i) {
      sort(colidx_temp + rowcnt[i], colidx_temp + rowcnt[i + 1]);
    }

    for (size_t i = 0; i < count; ++i) {
      int x = cooRowidx[i], y = cooColidx[i];
      if (x != y) {
        if (!binary_search(
                colidx_temp + rowcnt[y - 1], colidx_temp + rowcnt[y], x)) {
            cooRowidx[real_count] = y;
            cooColidx[real_count] = x;
            cooValues[real_count] = cooValues[i];
            ++real_count;
        }
      }
    }
    count = real_count;

    free(rowcnt);
    free(colidx_temp);
  }

/*********************** Convert COO to CSR*****************/
  // Allocate memory for CSR
  matrix->values = MALLOC(double, count);
  matrix->colidx = MALLOC(int, count);
  matrix->rowptr = MALLOC(int, count);
  if (!matrix->values || !matrix->colidx || !matrix->rowptr) {
      fprintf(stderr, "Failed to allocate memory\n");
      return false;
  }
  
  int i, l;
  int outBase = outputbase;
  int sort = 1;
  
  #pragma omp parallel for
  for (i = 0; i <= m; i++) matrix->rowptr[i] = 0;
  
  // determine row lengths
  for (i = 0; i < nnz; i++) matrix->rowptr[cooRowidx[i]]++;

  for (i = 0; i < m; i++) matrix->rowptr[i+1] += matrix->rowptr[i];

  // go through the structure  once more. Fill in output matrix.
  for (l = 0; l < nnz; l++) {
    i = matrix->rowptr[cooRowidx[l] - 1];
    matrix->values[i] = cooValues[l];
    matrix->colidx[i] = cooColidx[l] - 1 + outBase;
    matrix->rowptr[cooRowidx[l] - 1]++;
  }

  // shift back rowptr
  for (i = m; i > 0; i--) matrix->rowptr[i] = matrix->rowptr[i-1] + outBase;

  matrix->rowptr[0] = outBase;
  
  if (sort) {
      #pragma omp parallel for
      for (i=0; i < m; i++){
          qsort(matrix->colidx, matrix->values, matrix->rowptr[i] - outBase, matrix->rowptr[i+1] - 1 - outBase);
      }
  }

  matrix->m = m;
  matrix->n = n;
  matrix->nnz = count;
  
// Free memory
  FREE(cooRowidx);
  FREE(cooColidx);
  FREE(cooValues);
  return true;
}




/********************************************************************//**
 * loadCooMatrix load matrix in a matrix market format
 * 
 * @param[in] file
 *            The input file to read the matrix from
 * 
 * @param[out] matrix
 *           The data structure that contains the sparse  matrix
 *           in COO format 
 *
 * @param[in] outputbase
 *          - 0: the indexing starts from 0 (for fortran)
 *          - 1: the indexing starts from 1 (for C/C++)
 *
 * @param[in] transpose
 *           - 0: no transpose
 *           - 1:transpose
 *
 ***********************************************************************/
int loadCooMatrix(const char *file, COO *matrix, int outputbase, int transpose)
{
  FILE *fp=fopen(file, "r");
  if (NULL == fp) {
    fprintf(stderr, "Failed to open file %s\n", file);
    exit(-1);
  }

  // read banner
  MM_typecode matcode;
  if (mm_read_banner (fp, &matcode) != 0) {
    fprintf(stderr, "Error: could not process Matrix Market banner.\n");
    fclose(fp);
    return false;
  }

  if (!mm_is_valid (matcode) || mm_is_array (matcode) || mm_is_dense (matcode) ) {
    fprintf(stderr, "Error: only support sparse and real matrices.\n");
    fclose(fp);
    return false;
  }
  bool pattern = mm_is_pattern (matcode);

  // read sizes
  int m, n;
  int nnz; // # of non-zeros specified in the file
  if (mm_read_mtx_crd_size(fp, &m, &n, &nnz) !=0) {
    fprintf(stderr, "Error: could not read matrix size.\n");
    fclose(fp);
    return false;
  }

  if (transpose) {
      swap(m, n);
  }

  size_t count;
  if (mm_is_symmetric (matcode) == 1) {
    matrix->isSymmetric = 1;
    count = 2L*nnz;
  }
  else {
      matrix->isSymmetric = 0;
      count = nnz;
  }

  
  // allocate memory
  matrix->values = MALLOC(double, count);
  matrix->Colidx = MALLOC(int, count);
  matrix->Rowidx = MALLOC(int, count);
  if (!matrix->values || !matrix->Colidx || !matrix->Rowidx) {
    fprintf(stderr, "Failed to allocate memory\n");
    fclose(fp);
    return false;
  }

  int *colidx_temp=NULL, *rowcnt = NULL;

  if (matrix->isSymmetric) {
      colidx_temp = MALLOC(int, count);
      rowcnt = MALLOC(int, m + 1);
      if (!colidx_temp || !rowcnt) {
          fprintf(stderr, "Failed to allocate memory\n");
          fclose(fp);
          return false;
      }
      memset(rowcnt, 0, sizeof(int)*(m + 1));
  }

  // read values
  count = 0;
  int lines = 0;
  int x, y;
  double real, imag;
  int base = 1;
  while (mm_read_mtx_crd_entry (fp, &x, &y, &real, &imag, matcode) == 0) {
      if (transpose) swap(x, y);

    if (x > m || y > n) {
      fprintf(stderr, "Error: (%d %d) coordinate is out of range.\n", x, y);
      fclose(fp);
      return false;
    }

    matrix->Rowidx[count] = x;
    matrix->Colidx[count] = y;
    matrix->values[count] = pattern ? 1 : real;
    if (0 == x || 0 == y) base = 0;

    ++count;
    ++lines;
    if (matrix->isSymmetric) rowcnt[x]++;
      // this is not a bug. we're intentionally indexing rowcnt[x] instead of rowcnt[x-1]
  }
  fclose(fp);
  
  if (0 == base) {
    for (size_t i = 0; i < count; ++i) {
      matrix->Rowidx[i]++;
      matrix->Colidx[i]++;
    }
    if (matrix->isSymmetric) {
        for (int i = m; i > 0; --i) {
            rowcnt[i] = rowcnt[i - 1];
        }
    } 
  }

  if (lines != nnz) {
    fprintf(stderr, "Error: nnz (%d) specified in the header doesn't match with # of lines (%d) in file %s\n",
            nnz, lines, file);
    return false;
  }
  
  if (matrix->isSymmetric) {
    // add transposed elements only if it doesn't exist
    size_t real_count = count;
    // preix-sum
    for (int i = 0; i < m; ++i) {
      rowcnt[i + 1] += rowcnt[i];
    }
    for (size_t i = 0; i < count; ++i) {
      int j = rowcnt[matrix->Rowidx[i] - 1];
      colidx_temp[j] = matrix->Colidx[i];
      rowcnt[matrix->Rowidx[i] - 1]++;
    }
    for (int i = m; i > 0; --i) {
      rowcnt[i] = rowcnt[i - 1];
    }
    rowcnt[0] = 0;

#pragma omp parallel for
    for (int i = 0; i < m; ++i) {
      sort(colidx_temp + rowcnt[i], colidx_temp + rowcnt[i + 1]);
    }

    for (size_t i = 0; i < count; ++i) {
      int x = matrix->Rowidx[i], y = matrix->Colidx[i];
      if (x != y) {
        if (!binary_search(
                colidx_temp + rowcnt[y - 1], colidx_temp + rowcnt[y], x)) {
            matrix->Rowidx[real_count] = y;
            matrix->Colidx[real_count] = x;
            matrix->values[real_count] = matrix->values[i];
            ++real_count;
        }
      }
    }
    count = real_count;

    free(rowcnt);
    free(colidx_temp);
  }
  
  matrix->m = m;
  matrix->n = n;
  matrix->nnz = count;
  return true;
}




