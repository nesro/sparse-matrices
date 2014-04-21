/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
 * https://github.com/nesro/sparse-matrices
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#include "utils.h"
#include "bsr_matrix.h"
#include "virtual_matrix.h"
#include "mm_load.h"

static vm_vmt_t bsr_vmt = { /**/
(reset_t) NULL, /**/
(free_t) bsr_free, /**/
(mm_load_t) NULL,/**/
(mm_save_t) NULL, /**/
(print_t) NULL, /**/
(compare_t) NULL, /**/
(distance_t) NULL, /**/
(convert_t) NULL, /**/
(mul_t) NULL, /**/
};

void bsr_vm_init(bsr_t **bsr, va_list va) {

	int va_flag = 0;
	int width;
	int height;
	int nnz;
	int b_size;
	int b_cnt;

	b_cnt = va_get_int(va, 1, &va_flag);
	b_size = va_get_int(va, 1, &va_flag);
	nnz = va_get_int(va, 1, &va_flag);
	height = va_get_int(va, -1, &va_flag);
	width = va_get_int(va, -1, &va_flag);

	bsr_init(bsr, width, height, nnz, b_size, b_cnt);
}

static double bsr_load_mm(bsr_t **bsr, const char *filename, int b_size) {

	double start_time;
	double end_time;
	mm_file_t *mm_file;
	int i;
	int blocks;
	int block_row;
	int block_col;
	int *tmp_last_col;
	int blocks_in_row;
	datatype_t **tmp_block_beginnigs;

	mm_file = mm_load(filename);

	assert(mm_file->width % b_size == 0);
	assert(mm_file->height % b_size == 0);

	start_time = omp_get_wtime();

	blocks_in_row = ((mm_file->width / b_size) + 1);

	/* compute how many blocks is in the mtx matrix */
	blocks = 0;
	tmp_last_col = malloc(blocks_in_row * sizeof(int));
	assert(tmp_last_col);
	memset(tmp_last_col, -1, blocks_in_row * sizeof(int));
	for (i = 0; i < mm_file->nnz; i++) {

		block_row = mm_file->data[i].row / b_size;
		block_col = mm_file->data[i].col / b_size;

		_s_debugf(BSR_DEBUG, "i=%d v=%lf bc=%d ff=%d\n", i,
				mm_file->data[i].value, block_col, blocks_in_row);

		if (tmp_last_col[block_col] != block_row) {
			blocks++;
			_s_debug(BSR_DEBUG, "found a block!\n");
			tmp_last_col[block_col] = block_row;
		}
	}
	free(tmp_last_col);

	/* construct a bsr matrix */
	bsr_init(bsr, mm_file->width, mm_file->height, mm_file->nnz, b_size,
			blocks);
	blocks = 0;
	tmp_block_beginnigs = calloc(blocks_in_row, sizeof(datatype_t *));
	assert(tmp_block_beginnigs != NULL);
	for (i = 0; i < mm_file->nnz; i++) {

		block_row = mm_file->data[i].row / b_size;
		block_col = mm_file->data[i].col / b_size;

		if (tmp_block_beginnigs[block_col] == NULL) {
			tmp_block_beginnigs[block_col] = &((*bsr)->v[b_size * b_size
					* blocks]);

			printf("%d\n", blocks);
			(*bsr)->ci[blocks] = block_col;
			blocks++;
		}

		/*
		 * Place the item in the right place in a block.
		 * FIXME: not right place
		 */
		printf("v=%p tbb=%p r=%d c=%d br=%d bc=%d\n", (void *) (*bsr)->v,
				(void *) tmp_block_beginnigs[block_col], (block_row % b_size),
				(block_col % b_size), block_row, block_col);

		*(tmp_block_beginnigs[block_col] + /**/
		b_size * (mm_file->data[i].row % b_size) + /**/
		(mm_file->data[i].col % b_size)) += mm_file->data[i].value;

		(*bsr)->rp[mm_file->data[i].row + 1] = blocks;
	}
	free(tmp_block_beginnigs);

	end_time = omp_get_wtime();
	mm_free(mm_file);
	return end_time - start_time;
}

void bsr_from_mm(bsr_t **bsr, const char *mm_filename, va_list va) {

	int va_flag = 0;
	int b_size;

	b_size = va_get_int(va, -1, &va_flag);

	assert(b_size > 0);

	bsr_load_mm(bsr, mm_filename, b_size);
}

void bsr_init(bsr_t **bsr, int width, int height, int nnz, int b_size,
		int b_cnt) {

	if (bsr != NULL && *bsr != NULL)
		bsr_free(*bsr);

	*bsr = calloc(1, sizeof(bsr_t));

	(*bsr)->_.type = BSR;
	(*bsr)->_.f = bsr_vmt;
	(*bsr)->_.w = width;
	(*bsr)->_.h = height;
	(*bsr)->_.nnz = nnz;

	(*bsr)->bs = b_size;
	(*bsr)->bc = b_cnt;

	(*bsr)->v = calloc(((*bsr)->bc * b_size * b_size), sizeof(datatype_t));
	assert((*bsr)->v != NULL);

	(*bsr)->rp = malloc(((height / b_size) + 1) * sizeof(int));
	assert((*bsr)->rp != NULL);

	(*bsr)->ci = calloc((*bsr)->bc, sizeof(int));
	assert((*bsr)->ci != NULL);
}

void bsr_free(bsr_t *bsr) {
	free(bsr->ci);
	free(bsr->rp);
	free(bsr->v);
	free(bsr);
}
