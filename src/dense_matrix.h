/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
 * https://github.com/nesro/sparse-matrices
 */

#include "vector.h"
#include "matrix.h"

#ifndef DENSE_MATRIX_H_
#define DENSE_MATRIX_H_

typedef struct dense_matrix {
	matrix_t info;
	datatype_t **v;

	/* private */
	datatype_t *rows_block;
} dense_matrix_t;

int dense_matrix_init(dense_matrix_t *dense_matrix, int width, int height);
void dense_matrix_matrix_mul(dense_matrix_t *a, dense_matrix_t *b,
	dense_matrix_t *c);
void dense_matrix_vector_mul(dense_matrix_t *a, vector_t *b, vector_t *c);
void dense_to_html(dense_matrix_t *a, char *filename);
void dense_matrix_free(dense_matrix_t *dense_matrix);

int dense_matrix_load_mm(dense_matrix_t *dense_matrix, const char *filename);

void dense_matrix_hash_save(dense_matrix_t *dense_matrix, const char *filename);

int dense_matrix_compare(dense_matrix_t *a, dense_matrix_t *b);

#endif /* DENSE_MATRIX_H_ */
