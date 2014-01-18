/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
 * https://github.com/nesro/sparse-matrices
 */

#include "utils.h"
#include "matrix.h"
#include "csr_matrix.h"

#ifndef COO_MATRIX_H_
#define COO_MATRIX_H_

typedef struct coo_matrix {
	matrix_t info;
	datatype_t *values;
	int *col;
	int *row;
} coo_matrix_t;

int coo_matrix_init(coo_matrix_t *coo_matrix, int nnz, int width, int height);
void coo_matrix_free(coo_matrix_t *coo_matrix);

int coo_matrix_load_mm(coo_matrix_t *coo_matrix, const char *filename);

void coo_to_csr(coo_matrix_t *coo_matrix, csr_matrix_t *csr_matrix);
void coo_to_vector(coo_matrix_t *coo_matrix, vector_t *vector);
void coo_to_dense(coo_matrix_t *coo_matrix, dense_matrix_t *dense_matrix);



#endif /* COO_MATRIX_H_ */
