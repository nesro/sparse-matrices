/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
 * https://github.com/nesro/sparse-matrices
 */

#ifndef SPARSE_MATRIX_H_
#define SPARSE_MATRIX_H_

typedef struct matrix {
	int w; /* width */
	int h; /* height */
	int nnz; /* number of nonzero element */
} matrix_t;

#endif /* SPARSE_MATRIX_H_ */
