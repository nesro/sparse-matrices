/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
 * https://github.com/nesro/sparse-matrices
 */

#ifndef QT_MATRIX_H_
#define QT_MATRIX_H_

#include "csr_matrix.h"
#include "dense_matrix.h"
#include "mmio.h"
#include "utils.h"

#define PRINT 0

typedef struct qt_submatrix {
	int x;
	int y;

	int iv; /* index in the values array */
	int irp; /* index in the row pointer array */

	int nnz;
	int rows;

	/* for loading purposes only */
	int lr; /* last row */
	int i; /* have i items from nnz */
} qt_submatrix_t;

typedef struct qt_node {
	qt_submatrix_t *sm; /* submatrix */
	struct qt_node *tl; /* top left node */
	struct qt_node *tr; /* top right node */
	struct qt_node *bl; /* bottom left node */
	struct qt_node *br; /* bottom right node */
} qt_node_t;

typedef struct qt_matrix {

	const char *filename;

	datatype_t *v;
	int *ci;
	int *rp;

	int blocks;

	int height; /* height of the tree */

	int sm_size;
	matrix_t info;
	qt_node_t *root;

	/* for loading purposes only */
	int last_index_v;
	int last_index_rp;
} qt_matrix_t;

void qt_matrix_init(qt_matrix_t *qt_matrix, int width, int height, int nnz,
	int sm_size);
void qt_matrix_free(qt_matrix_t *qt_matrix);
double qt_matrix_matrix_mul(qt_matrix_t *a, qt_matrix_t *b, dense_matrix_t *c);

double qt_matrix_load_mm(qt_matrix_t *qt_matrix, const char *filename,
	int sm_size);

void qt_matrix_to_dense(qt_matrix_t *qt_matrix, dense_matrix_t *dense_matrix);

qt_submatrix_t *qt_submatrix_get(qt_matrix_t *qt_matrix, int i, int j);

void qt_matrix_print(qt_matrix_t *qt_matrix);

time_record_t qt_matrix_mm_mul(const char *matrix_a, const char *matrix_b, int sm_size);

#endif /* QT_MATRIX_H_ */
