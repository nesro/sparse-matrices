/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
 * https://github.com/nesro/sparse-matrices
 */

#include "virtual_matrix.h"

#ifndef BSR_MATRIX_H_
#define BSR_MATRIX_H_

#define BSR_DEBUG 1

typedef struct bsr_matrix {
	vm_t _;
	int bs; /* block size */
	int bc; /* block count */
	datatype_t *v;
	int *rp;
	int *ci;
} bsr_t;

void bsr_vm_init(bsr_t **bsr, va_list va);
void bsr_from_mm(bsr_t **bsr, const char *mm_filename, va_list va);
void bsr_init(bsr_t **bsr, int width, int height, int nnz, int b_size,
		int b_cnt);
void bsr_free(bsr_t *bsr);

#endif /* BSR_MATRIX_H_ */
