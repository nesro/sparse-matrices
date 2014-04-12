/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
 * https://github.com/nesro/sparse-matrices
 */

#ifndef QT_MATRIX_H_
#define QT_MATRIX_H_

#include "virtual_matrix.h"
#include "csr_matrix.h"
#include "den_matrix.h"
#include "mmio.h"
#include "utils.h"

#define PRINT 0
#define QDT_CSR 1

/*****************************************************************************/
typedef struct qdt_submatrix {
	int x;
	int y;

	/* FIXME: use pointers instead of indexes? */
#if QDT_CSR
	int iv; /* index in the values array */
	int irp; /* index in the row pointer array */
	int nnz;
	int rows;
	/* for loading purposes only */
	int lr; /* last row */
	int i; /* have i items from nnz */
#else /* QDT_CSR */
vm_t m;
#endif /* QDT_CSR */
} qdt_submatrix_t;

typedef struct qdt_node {
qdt_submatrix_t *sm; /* submatrix */
struct qdt_node *tl; /* top left node */
struct qdt_node *tr; /* top right node */
struct qdt_node *bl; /* bottom left node */
struct qdt_node *br; /* bottom right node */
} qdt_node_t;

typedef struct qdt_matrix {
vm_t _;

const char *filename;

#if QDT_CSR
datatype_t *v;
int *ci;
int *rp;
#endif /* QDT_CSR */

int blocks;
int height; /* height of the tree */

int sm_size;
//	matrix_t info;
qdt_node_t *root;

/* for loading purposes only */
int last_index_v;
int last_index_rp;
} qdt_matrix_t;

void qdt_vm_init(qdt_matrix_t **qdt, va_list va);
void qdt_init(qdt_matrix_t **qt_matrix, int width, int height, int nnz,
	int sm_size);
void qdt_free(qdt_matrix_t *qdt_matrix);
double qdt_mul(qdt_matrix_t *a, qdt_matrix_t *b, den_matrix_t **c, char);

void qdt_from_mm(qdt_matrix_t **, const char *, va_list);
double qdt_load_mm(qdt_matrix_t **qt_matrix, const char *filename, int sm_size);

void qdt_to_dense(qdt_matrix_t *qt_matrix, den_matrix_t *dense_matrix);

qdt_submatrix_t *qdt_submatrix_get(qdt_matrix_t *qt_matrix, int i, int j);

void qdt_print(qdt_matrix_t *qt_matrix);

vm_t *qdt_convert(qdt_matrix_t*, vm_type_t);

time_record_t qt_matrix_mm_mul(const char *matrix_a, const char *matrix_b,
	int sm_size);

#endif /* QT_MATRIX_H_ */
