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
#include "den_matrix.h"
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
(convert_t) bsr_convert, /**/
(mul_t) bsr_mul, /**/
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
	int last_row;

	mm_file = mm_load(filename, 1);

	printf("bsr w=%d b_size=%d\n", mm_file->width, b_size);

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

		_s_debugf(BSR_DEBUG, "i=%d v="DPF" bc=%d ff=%d\n", i,
				mm_file->data[i].value, block_col, blocks_in_row);

		if (tmp_last_col[block_col] != block_row) {
			blocks++;
			_s_debug(BSR_DEBUG, "found a block!\n");
			tmp_last_col[block_col] = block_row;
		}
	}
	free(tmp_last_col);

	_s_debug(BSR_DEBUG, "Let's construct a BSR matrix:\n");

	/* construct a bsr matrix */
	bsr_init(bsr, mm_file->width, mm_file->height, mm_file->nnz, b_size,
			blocks);
	blocks = 0;
	tmp_block_beginnigs = calloc(blocks_in_row, sizeof(datatype_t *));
	assert(tmp_block_beginnigs != NULL);
	(*bsr)->rp[0] = 0;
	last_row = -1;
	for (i = 0; i < mm_file->nnz; i++) {

		block_row = mm_file->data[i].row / b_size;
		block_col = mm_file->data[i].col / b_size;

		if (block_row > last_row) {
			if (last_row == -1) {
				last_row = block_row;
			} else {
				last_row = block_row;
				_s_debugf(BSR_DEBUG, "deleting! block_row=%d\n", block_row);
				memset(tmp_block_beginnigs, 0,
						blocks_in_row * sizeof(datatype_t *));
			}
		}

		if (tmp_block_beginnigs[block_col] == NULL) {
			tmp_block_beginnigs[block_col] = &((*bsr)->v[b_size * b_size
					* blocks]);
			_s_debug(BSR_DEBUG, "found a block!\n");
			_s_debugf(BSR_DEBUG, "blocks=%d\n", blocks);
			(*bsr)->ci[blocks] = block_col;
			blocks++;
		}

		/*
		 * Place the item in the right place in a block.
		 */
		_s_debugf(BSR_DEBUG, "it="DPF" v=%p tbb=%p r=%d c=%d br=%d bc=%d\n",
				mm_file->data[i].value, (void * ) (*bsr)->v,
				(void * ) tmp_block_beginnigs[block_col], (block_row % b_size),
				(block_col % b_size), block_row, block_col);

		*(tmp_block_beginnigs[block_col] + /**/
		(mm_file->data[i].row % b_size) * b_size + /**/
		(mm_file->data[i].col % b_size)) += mm_file->data[i].value;

		(*bsr)->rp[block_row + 1] = blocks;
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
	(*bsr)->_.object_size += sizeof(bsr_t);

	(*bsr)->_.type = BSR;
	(*bsr)->_.f = bsr_vmt;
	(*bsr)->_.w = width;
	(*bsr)->_.h = height;
	(*bsr)->_.nnz = nnz;

	(*bsr)->bs = b_size;
	(*bsr)->bc = b_cnt;

	/**************************************************************************/


//	printf("test a\n");
//
//	double *d;
//	d = malloc(45988380672);
//	assert(d != NULL);
//
//	printf("test b\n");


	(*bsr)->v = calloc(((*bsr)->bc * b_size * b_size), sizeof(datatype_t));
	assert((*bsr)->v != NULL);

	_s_debugf(0, "calloc = %ld maxsizet=%zu iwant=%zu\n",
			((*bsr)->bc * b_size * b_size), ((size_t )-1),
			(size_t)(((*bsr)->bc * b_size * b_size) * sizeof(datatype_t)));


	(*bsr)->_.object_size += ((*bsr)->bc * b_size * b_size)
			* sizeof(datatype_t);
	/**************************************************************************/

	/* XXX: this should be malloc'd */
	(*bsr)->rp = calloc(((height / b_size) + 1), sizeof(int));
	(*bsr)->_.object_size += ((height / b_size) + 1) * sizeof(int);
	assert((*bsr)->rp != NULL);

	(*bsr)->ci = calloc((*bsr)->bc, sizeof(int));
	(*bsr)->_.object_size += (*bsr)->bc * sizeof(int);
	assert((*bsr)->ci != NULL);
}

void bsr_free(bsr_t *bsr) {
	free(bsr->ci);
	free(bsr->rp);

	free(bsr->v);

	free(bsr);
}

vm_t *bsr_convert(bsr_t *bsr, vm_type_t type) {

	vm_t *vm = NULL;
	int i;
	int j;
	int k;
	int l;

	switch (type) {
	/* Convert a BSR matrix to a DENSE matrix. */
	case DEN:
		vm_create(&vm, DEN, 1, bsr->_.w, bsr->_.h);

		for (i = 0; i < bsr->bc * bsr->bs * bsr->bs; i++) {
			_s_debugf(BSR_DEBUG, "i=%d v="DPF"\n", i, bsr->v[i]);
		}

		for (i = 0; i < ((bsr->_.h / bsr->bs) + 1); i++) {
			_s_debugf(BSR_DEBUG, "i=%d, rp=%d\n", i, bsr->rp[i]);
		}

		for (i = 0; i < bsr->bc; i++) {
			_s_debugf(BSR_DEBUG, "i=%d, ci=%d\n", i, bsr->ci[i]);
		}

		for (i = 0; i < (bsr->_.h / bsr->bs); i++) {
			_s_debugf(BSR_DEBUG, "cvt i=%d\n", i);
			for (j = bsr->rp[i]; j < bsr->rp[i + 1]; j++) {
				_s_debugf(BSR_DEBUG, "cvt j=%d col=%d\n", j, bsr->ci[j]);
				for (k = 0; k < bsr->bs; k++) {
					for (l = 0; l < bsr->bs; l++) {
						_s_debugf(BSR_DEBUG, "c[%d][%d]="DPF"\n",
								(i * bsr->bs) + k, (j * bsr->bs) + l,
								bsr->v[j * (bsr->bs * bsr->bs) + (k * bsr->bs)
										+ l]);
						((den_matrix_t *) vm)->v /**/
						[(i * bsr->bs) + k] /**/
						[(bsr->ci[j] * bsr->bs) + l] = /**/
						bsr->v[j * (bsr->bs * bsr->bs) + (k * bsr->bs) + l];
					}
				}
			}
		}
		break;
	default:
		fdie("unknown format to convert %d\n", type);
		break;
	}

	return vm;
}

/******************************************************************************/

static inline double mul_bsr_vec(const bsr_t *a, const vec_t *b, vec_t *c) {

	double start_time;
	int i;
	int j;
	int l;
	int m;

	_s_debugf(1, "w=%d h=%d\n", a->_.w, b->_.h);
	assert(a->_.w == b->_.h);

	start_time = omp_get_wtime();

	_s_debug(BSR_DEBUG, ".\n");

	for (i = 0; i < (a->_.h / a->bs); i++) {
		for (j = a->rp[i]; j < a->rp[i + 1]; j++) {
			for (l = 0; l < a->bs; l++) {
				for (m = 0; m < a->bs; m++) {

					c->v /**/
					[(i * a->bs) + l] /**/
					+= /**/
					a->v[j * (a->bs * a->bs) + (l * a->bs) + m] * /**/
					b->v[(a->ci[j] * a->bs) + m];
				}
			}
		}
	}

	return omp_get_wtime() - start_time;
}

static inline double mul_bsr_bsr(const bsr_t *a, const bsr_t *b,
		den_matrix_t *c) {

	double start_time;
	int i;
	int j;
	int k;
	int l;
	int m;
	int n;

	assert(a->_.w == b->_.w);
	assert(a->_.h == b->_.h);
	assert(a->bs == b->bs);

	start_time = omp_get_wtime();

	_s_debug(BSR_DEBUG, "--- multiplication ---\n");

	/*
	 * For each row of blocks
	 */
	for (i = 0; i < (a->_.h / a->bs); i++) {
		_s_debugf(BSR_DEBUG, "blockrow=%d, from=%d, to=%d\n", i, a->rp[i],
				a->rp[i + 1]);
		/*
		 * For each block in the row
		 */
		for (j = a->rp[i]; j < a->rp[i + 1]; j++) {
			_s_debugf(BSR_DEBUG,
					"j=%d a->ci[j]=%d, b->rp[a->ci[j]]=%d, b->rp[a->ci[j] + 1]=%d\n",
					j, a->ci[j], b->rp[a->ci[j]], b->rp[a->ci[j] + 1]);

			/*
			 * For each row in the matrix B corresponding to the column
			 * from the matrix A
			 */
			for (k = b->rp[a->ci[j]]; k < b->rp[a->ci[j] + 1]; k++) {
				_s_debugf(BSR_DEBUG, "k=%d\n", k);
				/*
				 * Multiply block like two dense matrices.
				 */
				for (l = 0; l < a->bs; l++) {
					for (m = 0; m < a->bs; m++) {
						for (n = 0; n < a->bs; n++) {
							c->v /**/
							[(i * a->bs) + l] /**/
							[(b->ci[k] * a->bs) + m] += /**/
							a->v[j * (a->bs * a->bs) + (l * a->bs) + n] * /**/
							b->v[k * (a->bs * a->bs) + (n * a->bs) + m];
						}
					}
				}
			}
		}
	}

	return omp_get_wtime() - start_time;
}

double bsr_mul(const bsr_t *a, const vm_t *b, vm_t **c, char flag /* unused */) {

	switch (b->type) {
	case VEC:
		vec_init((vec_t **) c, a->_.w);
		return mul_bsr_vec(a, (vec_t *) b, (vec_t *) *c);
	case BSR:
		den_matrix_init((den_matrix_t **) c, a->_.w, b->h, 1);
		return mul_bsr_bsr(a, (bsr_t *) b, (den_matrix_t *) *c);
	default:
		fprintf(stderr, "Unknown matrix type: %d\n", b->type);
		exit(1);
		return 0.;
	}
}
