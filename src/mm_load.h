/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
 * https://github.com/nesro/sparse-matrices
 */

#include "virtual_matrix.h"

#ifndef MM_LOAD_H_
#define MM_LOAD_H_

typedef struct mm_item {
	int row;
	int col;
	datatype_t value;
} mm_item_t;

typedef struct mm_file {
	int width; /* cols */
	int height; /* rows */
	int nnz;
	mm_item_t *data;
	int data_size;
} mm_file_t;

mm_file_t *mm_load(const char *);
void mm_free(mm_file_t *);

#endif /* MM_LOAD_H_ */
