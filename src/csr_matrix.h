/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
 * https://github.com/nesro/sparse-matrices
 */

#ifndef CSR_MATRIX_H_
#define CSR_MATRIX_H_

#include "matrix.h"
#include "dense_matrix.h"
#include "vector.h"

typedef struct csr_matrix {
	matrix_t info;
	datatype_t *v;
	int *ci;
	int *rp;
} csr_matrix_t;

void csr_matrix_init(csr_matrix_t *csr_matrix, int nnz, int width, int height);

double csr_matrix_matrix_mul_unrolled_parallel(csr_matrix_t *a, csr_matrix_t *b,
	dense_matrix_t *c);
double csr_matrix_matrix_mul_unrolled(csr_matrix_t *a, csr_matrix_t *b,
	dense_matrix_t *c);
double csr_matrix_matrix_mul_parallel(csr_matrix_t *a, csr_matrix_t *b,
	dense_matrix_t *c);
double csr_matrix_matrix_mul(csr_matrix_t *a, csr_matrix_t *b,
	dense_matrix_t *c);

double csr_matrix_vector_mul_unrolled_parallel(csr_matrix_t *a, vector_t *b,
	vector_t *c);
double csr_matrix_vector_mul_unrolled(csr_matrix_t *a, vector_t *b, vector_t *c);
double csr_matrix_vector_mul_parallel(csr_matrix_t *a, vector_t *b, vector_t *c);
double csr_matrix_vector_mul(csr_matrix_t *a, vector_t *b, vector_t *c);

void csr_to_html(csr_matrix_t *a, char *filename);
void csr_matrix_print(csr_matrix_t *a);

int csr_matrix_save(csr_matrix_t *csr_matrix, const char *filename);
int csr_matrix_load(csr_matrix_t *csr_matrix, const char *filename);
int csr_matrix_copy(csr_matrix_t *a, csr_matrix_t *b);

void csr_matrix_free(csr_matrix_t *csr_matrix);

void csr_to_dense(csr_matrix_t *csr_matrix, dense_matrix_t *dense_matrix);

int csr_matrix_load_mm(csr_matrix_t *csr_matrix, const char *filename);

int csr_matrix_generate(csr_matrix_t *csr_matrix, int width, int height,
	int nnz);

time_record_t csr_matrix_mm_mul(const char * matrix_a, const char *matrix_b);

#endif /* CSR_MATRIX_H_ */
