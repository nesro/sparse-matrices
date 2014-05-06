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
#include "csr_matrix.h"
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
	int i;
	int j;
	int k;
	int blk_i; /* block row */
	int blk_j; /* block col */
	int blk_c; /* block count */
	csr_t *csr = NULL;
	int *blk_start;
	datatype_t **blk_vp;
	start_time = omp_get_wtime();

	/*
	 * TODO: we are creating BSR from a CSR matrix. We should load COO from
	 * MM and convert it to BSR.
	 */
	csr_from_mm(&csr, filename, 0);

	/* count blocks */
	blk_start = malloc((csr->_.w / b_size + 1) * sizeof(int));
	memset(blk_start, -1, (csr->_.w / b_size + 1) * sizeof(int));
	blk_c = 0;
	for (i = 0; i < csr->_.h; i++) {
		blk_i = i / b_size;
		for (j = csr->rp[i]; j < csr->rp[i + 1]; j++) {
			blk_j = csr->ci[j] / b_size;

			if (blk_start[blk_j] != blk_i) {
				blk_start[blk_j] = blk_i;
				blk_c++;
			}
		}
	}
	free(blk_start);

	bsr_init(bsr, csr->_.w, csr->_.h, csr->_.nnz, b_size, blk_c);
	blk_c = 0;
	(*bsr)->rp[0] = 0;
	blk_vp = calloc((csr->_.w / b_size + 1), sizeof(datatype_t *));
	for (i = 0; i < (csr->_.h / b_size); i++) {
		for (j = 0; j < b_size; j++) {
			blk_i = b_size * i + j;
			for (k = csr->rp[blk_i]; k < csr->rp[blk_i + 1]; k++) {
				blk_j = csr->ci[k] / b_size;

				if (blk_vp[blk_j] == NULL) {
					blk_vp[blk_j] = (*bsr)->v + blk_c * b_size * b_size;
					(*bsr)->ci[blk_c] = blk_j;
					blk_c++;
				}

				*(blk_vp[blk_j] + (b_size * j) + (csr->ci[k] % b_size)) +=
						csr->v[k];
			}
		}

		for (j = csr->rp[i * b_size]; j < csr->rp[(i + 1) * b_size]; j++) {
			blk_vp[csr->ci[j] / b_size] = NULL;
		}

		(*bsr)->rp[i + 1] = blk_c;
	}

	free(blk_vp);
	csr->_.f.free((vm_t*) csr);
	end_time = omp_get_wtime();
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

	assert(bsr != NULL);

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

	(*bsr)->v = calloc(((*bsr)->bc * b_size * b_size), sizeof(datatype_t));
	assert((*bsr)->v != NULL);

	_s_debugf(0, "calloc = %ld maxsizet=%zu iwant=%zu\n",
			((*bsr)->bc * b_size * b_size), ((size_t )-1),
			(size_t )(((*bsr)->bc * b_size * b_size) * sizeof(datatype_t)));

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
