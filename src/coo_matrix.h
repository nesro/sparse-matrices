/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
 * https://github.com/nesro/sparse-matrices
 */

#include "utils.h"
#include "matrix.h"
#include "virtual_matrix.h"
#include "csr_matrix.h"
#include "den_matrix.h"

#ifndef COO_MATRIX_H_
#define COO_MATRIX_H_

typedef struct coo_matrix {
	vm_t _;
	datatype_t *v;
	int *c;
	int *r;
} coo_matrix_t;

void coo_vm_init(coo_matrix_t **, va_list);
void coo_init(coo_matrix_t **, int, int, int);
void coo_free(coo_matrix_t *);
vm_t *coo_convert(coo_matrix_t *, vm_type_t);
double coo_mul(const coo_matrix_t *, const coo_matrix_t *, den_matrix_t **,
		char);
double coo_from_mm(coo_matrix_t **, const char *, va_list);

//int coo_matrix_init(coo_matrix_t *coo_matrix, int nnz, int width, int height);
//void coo_matrix_free(coo_matrix_t *coo_matrix);
//int coo_matrix_load_mm(coo_matrix_t *coo_matrix, const char *filename);
//void coo_to_csr(coo_matrix_t *coo_matrix, csr_t *csr_matrix);
//void coo_to_vector(coo_matrix_t *coo_matrix, vector_t *vector);
//void coo_to_dense(coo_matrix_t *coo_matrix, den_matrix_t *dense_matrix);

#endif /* COO_MATRIX_H_ */
