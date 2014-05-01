/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
 * https://github.com/nesro/sparse-matrices
 */

#ifndef CSR_MATRIX_H_
#define CSR_MATRIX_H_

#include "virtual_matrix.h"

#define CSR_DEBUG 1

typedef struct csr_matrix {
	vm_t _;
	datatype_t *v;
	int *ci;
	int *rp;
} csr_t;

void csr_init(csr_t **csr, int width, int height, int nnz);
void csr_free(csr_t *csr);
void csr_vm_init(csr_t **csr, va_list va);
void csr_from_mm(csr_t **csr, const char *mm_filename, va_list va);
double csr_mul(const csr_t *a, const vm_t *b, vm_t **c, char flag /* unused */);

/*
 void csr_matrix_init(csr_t *csr_matrix, int nnz, int width, int height);

 double csr_matrix_matrix_mul_unrolled_parallel(csr_t *a, csr_t *b,
 den_matrix_t *c);
 double csr_matrix_matrix_mul_unrolled(csr_t *a, csr_t *b,
 den_matrix_t *c);
 double csr_matrix_matrix_mul_parallel(csr_t *a, csr_t *b,
 den_matrix_t *c);
 double csr_matrix_matrix_mul(csr_t *a, csr_t *b,
 den_matrix_t *c);

 double csr_matrix_vector_mul_unrolled_parallel(csr_t *a, vector_t *b,
 vector_t *c);
 double csr_matrix_vector_mul_unrolled(csr_t *a, vector_t *b, vector_t *c);
 double csr_matrix_vector_mul_parallel(csr_t *a, vector_t *b, vector_t *c);
 double csr_matrix_vector_mul(csr_t *a, vector_t *b, vector_t *c);

 void csr_to_html(csr_t *a, char *filename);
 void csr_matrix_print(csr_t *a);

 int csr_matrix_save(csr_t *csr_matrix, const char *filename);
 int csr_matrix_load(csr_t *csr_matrix, const char *filename);
 int csr_matrix_copy(csr_t *a, csr_t *b);

 void csr_matrix_free(csr_t *csr_matrix);

 void csr_to_dense(csr_t *csr_matrix, den_matrix_t *dense_matrix);

 int csr_matrix_load_mm(csr_t *csr_matrix, const char *filename);

 int csr_matrix_generate(csr_t *csr_matrix, int width, int height,
 int nnz);

 time_record_t csr_matrix_mm_mul(const char * matrix_a, const char *matrix_b);
 */
#endif /* CSR_MATRIX_H_ */
