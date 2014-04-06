/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
 * https://github.com/nesro/sparse-matrices
 */

#include <stdio.h>
#include <stdlib.h>

#include "vector.h"
#include "matrix.h"
#include "virtual_matrix.h"

#ifndef DENSE_MATRIX_H_
#define DENSE_MATRIX_H_

typedef struct den_matrix_init {
	int width;
	int height;
	int zero;
} den_matrix_init_t;

typedef struct den_matrix {
	vm_t _;
	datatype_t **v; /* values */

	/* private */
	datatype_t *rows_block;
	int strassen_block_treshold;
} den_matrix_t;

/***************************************************************************/

void den_vm_init(den_matrix_t **, va_list);
void den_from_mm(den_matrix_t **, const char *, va_list);
void den_matrix_init(den_matrix_t **, int, int, int);

int den_compare(den_matrix_t *, vm_t *);
double den_distance(den_matrix_t *, den_matrix_t *);

/***************************************************************************/

void den_matrix_print(den_matrix_t *);
void den_matrix_free(den_matrix_t *);

double den_mul_recursion(const den_matrix_t *a, const den_matrix_t *b,
		den_matrix_t *c);

double den_mul_unrolled_parallel(const den_matrix_t *, const den_matrix_t *,
		den_matrix_t *);
double den_mul_unrolled(const den_matrix_t *, const den_matrix_t *,
		den_matrix_t *);
double den_mul_parallel(const den_matrix_t *, const den_matrix_t *,
		den_matrix_t *);
double den_mul_naive(const den_matrix_t *, const den_matrix_t *, den_matrix_t *);
double den_mul_strassen_unrolled_parallel(const den_matrix_t *,
		const den_matrix_t *, den_matrix_t *);
double den_mul_strassen_unrolled(const den_matrix_t *, const den_matrix_t *,
		den_matrix_t *);
double den_mul_strassen_parallel(const den_matrix_t *, const den_matrix_t *,
		den_matrix_t *);
double den_mul_strassen(const den_matrix_t *, const den_matrix_t *,
		den_matrix_t *);

/*
 * This is called from virtual matrix.
 */
double mul(const den_matrix_t *, const den_matrix_t *, den_matrix_t **, char);

int den_count_nnz(const den_matrix_t *);

/***************************************************************************/

void dense_matrix_matrix_mul(den_matrix_t *a, den_matrix_t *b, den_matrix_t *c);
void dense_matrix_vector_mul(den_matrix_t *a, vector_t *b, vector_t *c);
void dense_to_html(den_matrix_t *a, char *filename);

int dense_matrix_load_mm(den_matrix_t *dense_matrix, const char *filename);

void dense_matrix_hash_save(den_matrix_t *dense_matrix, const char *filename);

int dense_matrix_compare(den_matrix_t *a, den_matrix_t *b);

/*****************************************************************************/

void den_offset_add(const den_matrix_t *a, const den_matrix_t *b,
		den_matrix_t *c, int ar, int ac, int br, int bc, int cr, int cc, int s);

void den_offset_addto(const den_matrix_t *a, den_matrix_t *c, int ar, int ac,
		int cr, int cc, int s);

void den_offset_sub(const den_matrix_t *a, const den_matrix_t *b,
		den_matrix_t *c, int ar, int ac, int br, int bc, int cr, int cc, int s);

#endif /* DENSE_MATRIX_H_ */
