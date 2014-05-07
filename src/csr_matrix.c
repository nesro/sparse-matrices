/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013,2014
 * https://github.com/nesro/sparse-matrices
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <ctype.h>

#include "omp.h"
#include "utils.h"
#include "mm_load.h"
#include "csr_matrix.h"
#include "coo_matrix.h"
#include "den_matrix.h"

/***************************************************************************/

static vm_vmt_t csr_vmt = { /**/
(reset_t) NULL, /**/
(free_t) csr_free, /**/
(mm_load_t) NULL,/**/
(mm_save_t) NULL, /**/
(print_t) NULL, /**/
(compare_t) NULL, /**/
(distance_t) NULL, /**/
(convert_t) NULL, /**/
(mul_t) csr_mul, /**/
};

/***************************************************************************/

void csr_init(csr_t **csr, int width, int height, int nnz) {

	*csr = calloc(1, sizeof(csr_t));
	(*csr)->_.object_size += sizeof(csr_t);
	assert(*csr != NULL);

	(*csr)->_.type = CSR;
	(*csr)->_.f = csr_vmt;
	(*csr)->_.w = width;
	(*csr)->_.h = height;
	(*csr)->_.nnz = nnz;

	(*csr)->v = malloc(nnz * sizeof(datatype_t));
	(*csr)->_.object_size += nnz * sizeof(datatype_t);
	assert((*csr)->v != NULL);

	(*csr)->ci = malloc(nnz * sizeof(int));
	(*csr)->_.object_size += nnz * sizeof(int);
	assert((*csr)->ci != NULL);

	/*
	 * We are adding to rp items to compute offset. Thus we need
	 * calloc it.
	 */
	(*csr)->rp = calloc(((*csr)->_.h + 1), sizeof(int));
	(*csr)->_.object_size += ((*csr)->_.h + 1) * sizeof(int);
	assert((*csr)->rp != NULL);
}

void csr_free(csr_t *csr) {

	free(csr->ci);
	free(csr->v);
	free(csr->rp);
	free(csr);
}

void csr_vm_init(csr_t **csr, va_list va) {

	int va_flag = 0;
	int width;
	int height;
	int nnz;

	nnz = va_get_int(va, 1, &va_flag);
	height = va_get_int(va, -1, &va_flag);
	width = va_get_int(va, -1, &va_flag);

	csr_init(csr, width, height, nnz);
}

void csr_from_mm(csr_t **csr, const char *mm_filename, va_list va) {

	mm_file_t *mm_file;
	int i;
	int tmp_sum;
	int tmp_row;

	int row;
	int dest;
	int last;

	/* XXX: FIXME: ssh! - for test purposes of course (8 */
	mm_file = mm_load(mm_filename, 1);

	csr_init(csr, mm_file->width, mm_file->height, mm_file->nnz);

	for (i = 0; i < mm_file->nnz; i++)
		(*csr)->rp[mm_file->data[i].row]++;

	tmp_sum = 0;
	for (i = 0; i < (*csr)->_.h; i++) {
		tmp_row = (*csr)->rp[i];
		(*csr)->rp[i] = tmp_sum;
		tmp_sum += tmp_row;
	}
	(*csr)->rp[(*csr)->_.h] = mm_file->nnz;

	for (i = 0; i < mm_file->nnz; i++) {
		row = mm_file->data[i].row;
		dest = (*csr)->rp[row];

		(*csr)->ci[dest] = mm_file->data[i].col;
		(*csr)->v[dest] = mm_file->data[i].value;

		(*csr)->rp[row]++;
	}

	last = 0;
	for (i = 0; i < (*csr)->_.h; i++) {
		tmp_row = (*csr)->rp[i];
		(*csr)->rp[i] = last;
		last = tmp_row;
	}

	mm_free(mm_file);
}

/******************************************************************************/

static double mul_csr_vec(const csr_t *a, const vec_t *b, vec_t *c) {

	double start_time;
	int r;
	int ac;

	start_time = omp_get_wtime();

	for (r = 0; r < a->_.w; r++)
		for (ac = a->rp[r]; ac < a->rp[r + 1]; ac++)
			c->v[r] += a->v[ac] * b->v[a->ci[ac]];

	return omp_get_wtime() - start_time;
}

static double mul_csr_csr(const csr_t *a, const csr_t *b, den_matrix_t *c) {

	double start_time;

	int r; /* row */
	int ac; /* col of A */
	int bc; /* col of B */

	start_time = omp_get_wtime();

	for (r = 0; r < a->_.h; r++) {
		for (ac = a->rp[r]; ac < a->rp[r + 1]; ac++) {
			for (bc = b->rp[a->ci[ac]]; bc < b->rp[a->ci[ac] + 1]; bc++) {
				c->v[r][b->ci[bc]] += a->v[ac] * b->v[bc];
			}
		}
	}

	return omp_get_wtime() - start_time;
}

double csr_mul(const csr_t *a, const vm_t *b, vm_t **c, char flag /* unused */) {

	switch (b->type) {
	case VEC:
		vec_init((vec_t **) c, a->_.w);
		return mul_csr_vec(a, (vec_t *) b, (vec_t *) *c);
	case CSR:
		den_matrix_init((den_matrix_t **) c, a->_.w, b->h, 1);
		return mul_csr_csr(a, (csr_t *) b, (den_matrix_t *) *c);
	default:
		fprintf(stderr, "Unknown matrix type: %d\n", b->type);
		exit(1);
		return 0.;
	}
}

///**
// * Algorithm has been taken from:
// * https://fedcsis.org/proceedings/2013/pliks/19.pdf
// */
//double csr_matrix_matrix_mul_unrolled_parallel(csr_matrix_t *a, csr_matrix_t *b,
//	den_matrix_t *c) {
//
//	double start_time;
//	double end_time;
//	int row;
//	int acol;
//	int acol2;
//	int bcol;
//	int bcol2;
//	int acol_i;
//	int acol2_i;
//	int bcol_i;
//	int bcol2_i;
//	int x;
//	int x2;
//
//	if (!den_matrix_init(c, a->info.h, b->info.w)) {
//		return -1;
//	}
//	start_time = omp_get_wtime();
//
//#pragma omp parallel for default(none) shared(a,b,c) \?
//	private(row,acol,acol2,bcol,bcol2,acol_i,acol2_i,bcol_i,bcol2_i,x,x2) \?
//	__OMP_NUM_THREADS__
//	for (row = 0; row < a->info.h - 1; row += 2) {
//		acol = a->rp[row];
//		acol2 = a->rp[row + 1];
//		acol_i = a->rp[row + 1] - a->rp[row];
//		acol2_i = a->rp[row + 2] - a->rp[row + 1];
//
//		if (acol_i <= acol2_i) {
//			for (; acol < a->rp[row + 1]; acol++, acol2++) {
//				x = a->ci[acol];
//				x2 = a->ci[acol2];
//				bcol = b->rp[x];
//				bcol2 = b->rp[x2];
//				bcol_i = b->rp[x + 1] - b->rp[x];
//				bcol2_i = b->rp[x2 + 1] - b->rp[x2];
//
//				if (bcol_i <= bcol2_i) {
//					for (; bcol < b->rp[x + 1]; bcol++, bcol2++) {
//						c->v[row][b->ci[bcol]] += a->v[acol] * b->v[bcol];
//						c->v[row + 1][b->ci[bcol2]] += a->v[acol2]
//							* b->v[bcol2];
//					}
//					for (; bcol2 < b->rp[x2 + 1]; bcol2++) {
//						c->v[row + 1][b->ci[bcol2]] += a->v[acol2]
//							* b->v[bcol2];
//					}
//				} else {
//					for (; bcol2 < b->rp[x2 + 1]; bcol++, bcol2++) {
//						c->v[row][b->ci[bcol]] += a->v[acol] * b->v[bcol];
//						c->v[row + 1][b->ci[bcol2]] += a->v[acol2]
//							* b->v[bcol2];
//					}
//					for (; bcol < b->rp[x + 1]; bcol++) {
//						c->v[row][b->ci[bcol]] += a->v[acol] * b->v[bcol];
//					}
//				}
//			}
//			for (; acol2 < a->rp[row + 2]; acol2++) {
//				x = a->ci[acol];
//				x2 = a->ci[acol2];
//				bcol = b->rp[x];
//				bcol2 = b->rp[x2];
//				bcol_i = b->rp[x + 1] - b->rp[x];
//				bcol2_i = b->rp[x2 + 1] - b->rp[x2];
//
//				for (; bcol2 < b->rp[x2 + 1]; bcol2++) {
//					c->v[row + 1][b->ci[bcol2]] += a->v[acol2] * b->v[bcol2];
//				}
//			}
//		} else {
//			for (; acol2 < a->rp[row + 2]; acol++, acol2++) {
//				x = a->ci[acol];
//				x2 = a->ci[acol2];
//				bcol = b->rp[x];
//				bcol2 = b->rp[x2];
//				bcol_i = b->rp[x + 1] - b->rp[x];
//				bcol2_i = b->rp[x2 + 1] - b->rp[x2];
//
//				if (bcol_i < bcol2_i) {
//					for (; bcol < b->rp[x + 1]; bcol++, bcol2++) {
//						c->v[row][b->ci[bcol]] += a->v[acol] * b->v[bcol];
//						c->v[row + 1][b->ci[bcol2]] += a->v[acol2]
//							* b->v[bcol2];
//					}
//					for (; bcol2 < b->rp[x2 + 1]; bcol2++) {
//						c->v[row + 1][b->ci[bcol2]] += a->v[acol2]
//							* b->v[bcol2];
//					}
//				} else { /* bcol_i >= bcol2_i */
//					for (; bcol2 < b->rp[x2 + 1]; bcol++, bcol2++) {
//						c->v[row][b->ci[bcol]] += a->v[acol] * b->v[bcol];
//						c->v[row + 1][b->ci[bcol2]] += a->v[acol2]
//							* b->v[bcol2];
//					}
//					for (; bcol < b->rp[x + 1]; bcol++) {
//						c->v[row][b->ci[bcol]] += a->v[acol] * b->v[bcol];
//					}
//				}
//			}
//			for (; acol < a->rp[row + 1]; acol++) {
//				x = a->ci[acol];
//				bcol = b->rp[x];
//				bcol_i = b->rp[x + 1] - b->rp[x];
//
//				for (; bcol < b->rp[x + 1]; bcol++) {
//					c->v[row][b->ci[bcol]] += a->v[acol] * b->v[bcol];
//				}
//			}
//		}
//	}
//
//	if (a->info.h % 2 == 1) {
//		for (acol = a->rp[a->info.h - 1]; acol < a->rp[a->info.h]; acol++) {
//			x = a->ci[acol];
//
//			for (bcol = b->rp[x]; bcol < b->rp[x + 1]; bcol++) {
//				x2 = b->ci[bcol];
//				c->v[a->info.h - 1][x2] += a->v[acol] * b->v[bcol];
//			}
//		}
//	}
//
//	end_time = omp_get_wtime();
//	return end_time - start_time;
//}
//
//double csr_matrix_matrix_mul_unrolled(csr_matrix_t *a, csr_matrix_t *b,
//	den_matrix_t *c) {
//
//	double start_time;
//	double end_time;
//	int row;
//	int acol;
//	int acol2;
//	int bcol;
//	int bcol2;
//	int acol_i;
//	int acol2_i;
//	int bcol_i;
//	int bcol2_i;
//	int x;
//	int x2;
//
//	if (!den_matrix_init(c, a->info.h, b->info.w)) {
//		return -1;
//	}
//	start_time = omp_get_wtime();
//
//	for (row = 0; row < a->info.h - 1; row += 2) {
//		acol = a->rp[row];
//		acol2 = a->rp[row + 1];
//		acol_i = a->rp[row + 1] - a->rp[row];
//		acol2_i = a->rp[row + 2] - a->rp[row + 1];
//
//		if (acol_i <= acol2_i) {
//			for (; acol < a->rp[row + 1]; acol++, acol2++) {
//				x = a->ci[acol];
//				x2 = a->ci[acol2];
//				bcol = b->rp[x];
//				bcol2 = b->rp[x2];
//				bcol_i = b->rp[x + 1] - b->rp[x];
//				bcol2_i = b->rp[x2 + 1] - b->rp[x2];
//
//				if (bcol_i <= bcol2_i) {
//					for (; bcol < b->rp[x + 1]; bcol++, bcol2++) {
//						c->v[row][b->ci[bcol]] += a->v[acol] * b->v[bcol];
//						c->v[row + 1][b->ci[bcol2]] += a->v[acol2]
//							* b->v[bcol2];
//					}
//					for (; bcol2 < b->rp[x2 + 1]; bcol2++) {
//						c->v[row + 1][b->ci[bcol2]] += a->v[acol2]
//							* b->v[bcol2];
//					}
//				} else {
//					for (; bcol2 < b->rp[x2 + 1]; bcol++, bcol2++) {
//						c->v[row][b->ci[bcol]] += a->v[acol] * b->v[bcol];
//						c->v[row + 1][b->ci[bcol2]] += a->v[acol2]
//							* b->v[bcol2];
//					}
//					for (; bcol < b->rp[x + 1]; bcol++) {
//						c->v[row][b->ci[bcol]] += a->v[acol] * b->v[bcol];
//					}
//				}
//			}
//			for (; acol2 < a->rp[row + 2]; acol2++) {
//				x = a->ci[acol];
//				x2 = a->ci[acol2];
//				bcol = b->rp[x];
//				bcol2 = b->rp[x2];
//				bcol_i = b->rp[x + 1] - b->rp[x];
//				bcol2_i = b->rp[x2 + 1] - b->rp[x2];
//
//				for (; bcol2 < b->rp[x2 + 1]; bcol2++) {
//					c->v[row + 1][b->ci[bcol2]] += a->v[acol2] * b->v[bcol2];
//				}
//			}
//		} else {
//			for (; acol2 < a->rp[row + 2]; acol++, acol2++) {
//				x = a->ci[acol];
//				x2 = a->ci[acol2];
//				bcol = b->rp[x];
//				bcol2 = b->rp[x2];
//				bcol_i = b->rp[x + 1] - b->rp[x];
//				bcol2_i = b->rp[x2 + 1] - b->rp[x2];
//
//				if (bcol_i < bcol2_i) {
//					for (; bcol < b->rp[x + 1]; bcol++, bcol2++) {
//						c->v[row][b->ci[bcol]] += a->v[acol] * b->v[bcol];
//						c->v[row + 1][b->ci[bcol2]] += a->v[acol2]
//							* b->v[bcol2];
//					}
//					for (; bcol2 < b->rp[x2 + 1]; bcol2++) {
//						c->v[row + 1][b->ci[bcol2]] += a->v[acol2]
//							* b->v[bcol2];
//					}
//				} else {
//					for (; bcol2 < b->rp[x2 + 1]; bcol++, bcol2++) {
//						c->v[row][b->ci[bcol]] += a->v[acol] * b->v[bcol];
//						c->v[row + 1][b->ci[bcol2]] += a->v[acol2]
//							* b->v[bcol2];
//					}
//					for (; bcol < b->rp[x + 1]; bcol++) {
//						c->v[row][b->ci[bcol]] += a->v[acol] * b->v[bcol];
//					}
//				}
//			}
//			for (; acol < a->rp[row + 1]; acol++) {
//				x = a->ci[acol];
//				bcol = b->rp[x];
//				bcol_i = b->rp[x + 1] - b->rp[x];
//
//				for (; bcol < b->rp[x + 1]; bcol++) {
//					c->v[row][b->ci[bcol]] += a->v[acol] * b->v[bcol];
//				}
//			}
//		}
//	}
//	if (row < a->info.h) {
//		for (acol = a->rp[row]; acol < a->rp[row + 1]; acol++) {
//			x = a->ci[acol];
//
//			for (bcol = b->rp[x]; bcol < b->rp[x + 1]; bcol++) {
//				x2 = b->ci[bcol];
//				c->v[row][x2] += a->v[acol] * b->v[bcol];
//			}
//		}
//	}
//
//	end_time = omp_get_wtime();
//	return end_time - start_time;
//}
//
//double csr_matrix_matrix_mul_parallel(csr_matrix_t *a, csr_matrix_t *b,
//	den_matrix_t *c) {
//
//	double start_time;
//	double end_time;
//
//	int r;
//	int ac;
//	int bc;
//
//	if (!den_matrix_init(c, a->info.h, b->info.w)) {
//		return -1;
//	}
//	start_time = omp_get_wtime();
//
//#pragma omp parallel for \?
//	default(none) \?
//	shared(a,b,c) \?
//	private(ac, bc) \?
//	__OMP_NUM_THREADS__
//	for (r = 0; r < a->info.h; r++)
//		for (ac = a->rp[r]; ac < a->rp[r + 1]; ac++)
//			for (bc = b->rp[a->ci[ac]]; bc < b->rp[a->ci[ac] + 1]; bc++)
//				c->v[r][b->ci[bc]] += a->v[ac] * b->v[bc];
//
//	end_time = omp_get_wtime();
//	return end_time - start_time;
//}
//
//double csr_matrix_matrix_mul(csr_matrix_t *a, csr_matrix_t *b,
//	den_matrix_t *c) {
//
//	double start_time;
//	double end_time;
//
//	int r;
//	int ac;
//	int bc;
//
//	if (!den_matrix_init(c, a->info.h, b->info.w)) {
//		return -1;
//	}
//
//	start_time = omp_get_wtime();
//
//	for (r = 0; r < a->info.h; r++)
//		for (ac = a->rp[r]; ac < a->rp[r + 1]; ac++)
//			for (bc = b->rp[a->ci[ac]]; bc < b->rp[a->ci[ac] + 1]; bc++)
//				c->v[r][b->ci[bc]] += a->v[ac] * b->v[bc];
//
//	end_time = omp_get_wtime();
//	return end_time - start_time;
//}
//
//double csr_matrix_vector_mul_unrolled_parallel(csr_matrix_t *a, vector_t *b,
//	vector_t *c) {
//
//	double start_time;
//	double end_time;
//	int row_ptr;
//	int frs; /* first row start */
//	int fre; /* first row end */
//	int sre; /* second row end */
//	int fri; /* first row items */
//	int sri; /* second row items*/
//	int frci; /* first row column items */
//	int srci; /* second row column items */
//	datatype_t sfr; /* sum - first row */
//	datatype_t ssr; /* sum - second row */
//
//	assert(a->info.w == b->size);
//	vector_init(c, b->size);
//	start_time = omp_get_wtime();
//
//#pragma omp parallel default(none) shared(a,b,c) \?
//	private(row_ptr,frs,fre,sre,fri,sri,frci,srci,sfr,ssr) 	__OMP_NUM_THREADS__
//	{
//
//#pragma omp for
//		for (row_ptr = 0; row_ptr < a->info.w - 1; row_ptr += 2) {
//
//			frs = a->rp[row_ptr];
//			fre = a->rp[row_ptr + 1];
//			sre = a->rp[row_ptr + 2];
//			fri = fre - frs;
//			sri = sre - fre;
//			sfr = 0;
//			ssr = 0;
//			frci = frs;
//			srci = fre;
//
//			if (fri < sri) {
//				while (frci < fre) {
//					sfr += a->v[frci] * b->v[a->ci[frci]];
//					ssr += a->v[srci] * b->v[a->ci[srci]];
//					frci++;
//					srci++;
//				}
//				while (srci < sre) {
//					ssr += a->v[srci] * b->v[a->ci[srci]];
//					srci++;
//				}
//			} else {
//				while (srci < sre) {
//					sfr += a->v[frci] * b->v[a->ci[frci]];
//					ssr += a->v[srci] * b->v[a->ci[srci]];
//					frci++;
//					srci++;
//				}
//				while (frci < fre) {
//					sfr += a->v[frci] * b->v[a->ci[frci]];
//					frci++;
//				}
//			}
//
//			c->v[row_ptr] = sfr;
//			c->v[row_ptr + 1] = ssr;
//		}
//
//	}
//
//	if (a->info.w % 2 == 1) {
//		sfr = 0;
//		frs = a->rp[a->info.w - 1];
//		fre = a->rp[a->info.w];
//		frci = frs;
//
//		while (frci < fre) {
//			sfr += a->v[frci] * b->v[a->ci[frci]];
//			frci++;
//		}
//
//		c->v[a->info.w - 1] = sfr;
//	}
//
//	end_time = omp_get_wtime();
//	return end_time - start_time;
//
//}
//
//double csr_matrix_vector_mul_unrolled(csr_matrix_t *a, vector_t *b, vector_t *c) {
//
//	double start_time;
//	double end_time;
//	int row_ptr;
//	int frs; /* first row start */
//	int fre; /* first row end */
//	int sre; /* second row end */
//	int fri; /* first row items */
//	int sri; /* second row items*/
//	int frci; /* first row column items */
//	int srci; /* second row column items */
//	datatype_t sfr; /* sum - first row */
//	datatype_t ssr; /* sum - second row */
//
//	assert(a->info.w == b->size);
//	vector_init(c, b->size);
//	start_time = omp_get_wtime();
//
//	frs = a->rp[0];
//	for (row_ptr = 0; row_ptr < a->info.w - 1; row_ptr += 2) {
//		fre = a->rp[row_ptr + 1];
//		sre = a->rp[row_ptr + 2];
//		fri = fre - frs;
//		sri = sre - fre;
//		sfr = 0;
//		ssr = 0;
//		frci = frs;
//		srci = fre;
//
//		if (fri < sri) {
//			while (frci < fre) {
//				sfr += a->v[frci] * b->v[a->ci[frci]];
//				ssr += a->v[srci] * b->v[a->ci[srci]];
//				frci++;
//				srci++;
//			}
//			while (srci < sre) {
//				ssr += a->v[srci] * b->v[a->ci[srci]];
//				srci++;
//			}
//		} else {
//			while (srci < sre) {
//				sfr += a->v[frci] * b->v[a->ci[frci]];
//				ssr += a->v[srci] * b->v[a->ci[srci]];
//				frci++;
//				srci++;
//			}
//			while (frci < fre) {
//				sfr += a->v[frci] * b->v[a->ci[frci]];
//				frci++;
//			}
//		}
//
//		c->v[row_ptr] = sfr;
//		c->v[row_ptr + 1] = ssr;
//		frs = sre;
//	}
//	if (a->info.w % 2 == 1) {
//		sfr = 0;
//		frs = a->rp[a->info.w - 1];
//		fre = a->rp[a->info.w];
//		frci = frs;
//
//		while (frci < fre) {
//			sfr += a->v[frci] * b->v[a->ci[frci]];
//			frci++;
//		}
//
//		c->v[a->info.w - 1] = sfr;
//	}
//
//	end_time = omp_get_wtime();
//	return end_time - start_time;
//}
//
//double csr_matrix_vector_mul_parallel(csr_matrix_t *a, vector_t *b, vector_t *c) {
//
//	double start_time;
//	double end_time;
//	int r;
//	int ci;
//
//	assert(a->info.w == b->size);
//	vector_init(c, b->size);
//	start_time = omp_get_wtime();
//
//#pragma omp parallel for default(none) shared(a,b,c) private(r,ci) \?
//	__OMP_NUM_THREADS__
//	for (r = 0; r < a->info.h; r++) {
//		for (ci = a->rp[r]; ci < a->rp[r + 1]; ci++) {
//			c->v[r] += a->v[ci] * b->v[a->ci[ci]];
//		}
//	}
//
//	end_time = omp_get_wtime();
//	return end_time - start_time;
//}
//
//double csr_matrix_vector_mul(csr_matrix_t *a, vector_t *b, vector_t *c) {
//
//	int row_ptr;
//	int row_start;
//	int row_end;
//	int col_ind;
//	datatype_t sum;
//	double start_time;
//	double end_time;
//
//	assert(a->info.w == b->size);
//	vector_init(c, b->size);
//	start_time = omp_get_wtime();
//
//	row_start = a->rp[0];
//	for (row_ptr = 0; row_ptr < a->info.w; row_ptr++) {
//
//		row_end = a->rp[row_ptr + 1];
//		sum = 0;
//
//		for (col_ind = row_start; col_ind < row_end; col_ind++) {
//			sum += a->v[col_ind] * b->v[a->ci[col_ind]];
//		}
//
//		c->v[row_ptr] = sum;
//		row_start = row_end;
//	}
//
//	end_time = omp_get_wtime();
//	return end_time - start_time;
//}
//

//
//int csr_matrix_save(csr_matrix_t *csr_matrix, const char *filename) {
//
//	FILE *fp;
//	int i;
//
//	fp = fopen(filename, "w");
//
//	if (fp == NULL) {
//		fprintf(stderr, "Cannot open file %s for writing\n", filename);
//		return 0;
//	}
//
//	fprintf(fp, "%d %d %d\n", csr_matrix->info.h, csr_matrix->info.w,
//		csr_matrix->info.nnz);
//
//	for (i = 0; i < csr_matrix->info.nnz; i++) {
//		fprintf(fp, DATATYPE_FORMAT " %d\n", csr_matrix->v[i],
//			csr_matrix->ci[i]);
//	}
//
//	for (i = 0; i <= csr_matrix->info.h; i++) {
//		fprintf(fp, "%d\n", csr_matrix->rp[i]);
//	}
//
//	fclose(fp);
//	return 1;
//}
//
//int csr_matrix_load(csr_matrix_t *csr_matrix, const char *filename) {
//
//	FILE *fp;
//	int i;
//	int width;
//	int height;
//	int nnz;
//
//	fp = fopen(filename, "r");
//
//	if (fp == NULL) {
//		fprintf(stderr, "Cannot open file %s for reading\n", filename);
//		goto error;
//	}
//
//	if (fscanf(fp, "%d %d %d\n", &height, &width, &nnz) != 3) {
//		fprintf(stderr, "Cannot read header of csr file %s.\n", filename);
//		goto error;
//	}
//
//	csr_matrix_init(csr_matrix, nnz, width, height);
//
//	for (i = 0; i < nnz; i++) {
//		if (fscanf(fp, "%lg %d\n", &csr_matrix->v[i], &csr_matrix->ci[i])
//			!= 2) {
//			fprintf(stderr,
//				"Cannot read values and col_indicies of csr file %s for i=%d.\n",
//				filename, i);
//			goto error;
//		}
//	}
//
//	for (i = 0; i <= height; i++) {
//		if (fscanf(fp, "%d\n", &csr_matrix->rp[i]) != 1) {
//			fprintf(stderr, "Cannot read row of csr file %s for i=%d.\n",
//				filename, i);
//			goto error;
//		}
//	}
//
//	fclose(fp);
//
//	return 1;
//
//	error:
//
//	fclose(fp);
//	fprintf(stderr, "ERROR in csr_matrix_load\n");
//	csr_matrix_free(csr_matrix);
//
//	return 0;
//}
//
//void csr_matrix_free(csr_matrix_t *csr_matrix) {
//	free(csr_matrix->v);
//	csr_matrix->v = NULL;
//
//	free(csr_matrix->rp);
//	csr_matrix->rp = NULL;
//
//	free(csr_matrix->ci);
//	csr_matrix->ci = NULL;
//}
//
//int csr_matrix_copy(csr_matrix_t *a, csr_matrix_t *b) {
//
//	int i;
//
//	csr_matrix_init(b, a->info.nnz, a->info.w, a->info.h);
//
//	for (i = 0; i < a->info.nnz; i++) {
//		b->v[i] = a->v[i];
//		b->ci[i] = a->ci[i];
//	}
//
//	for (i = 0; i <= a->info.h; i++) {
//		b->rp[i] = a->rp[i];
//	}
//
//	return 1;
//}
//
//void csr_to_html(csr_matrix_t *a, char *filename) {
//
//	den_matrix_t dense_matrix;
//
//	csr_to_dense(a, &dense_matrix);
//	dense_to_html(&dense_matrix, filename);
//	dense_matrix_free(&dense_matrix);
//}
//
//void csr_matrix_print(csr_matrix_t *a) {
//	int i;
//
//	printf("printing CSR matrix\n");
//	printf("width=%d, height=%d, nnz=%d\n", a->info.w, a->info.h, a->info.nnz);
//
//	for (i = 0; i < a->info.h + 1; i++) {
//		printf("a->rp[%d]=%d\n", i, a->rp[i]);
//	}
//
//	for (i = 0; i < a->info.nnz; i++) {
//		printf("a->v[%d]=%lf, a->ci[%d]=%d\n", i, a->v[i], i, a->ci[i]);
//	}
//
//	printf("end of CSR matrix\n");
//}
//
//void csr_to_dense(csr_matrix_t *csr_matrix, den_matrix_t *dense_matrix) {
//
//	int col_ind;
//	int row_ptr;
//
//	den_matrix_init(dense_matrix, csr_matrix->info.w, csr_matrix->info.h);
//
//	for (row_ptr = 0; row_ptr < csr_matrix->info.h; row_ptr++) {
//		for (col_ind = csr_matrix->rp[row_ptr];
//			col_ind < csr_matrix->rp[row_ptr + 1]; col_ind++) {
//			dense_matrix->v[row_ptr][csr_matrix->ci[col_ind]] =
//				csr_matrix->v[col_ind];
//		}
//	}
//}
//
//int csr_matrix_load_mm(csr_matrix_t *csr_matrix, const char *filename) {
//
//	coo_matrix_t coo_matrix;
//
//	coo_matrix_load_mm(&coo_matrix, filename);
//	coo_to_csr(&coo_matrix, csr_matrix);
//	coo_matrix_free(&coo_matrix);
//
//	return 1;
//}
//
//int csr_matrix_generate(csr_matrix_t *csr_matrix, int width, int height,
//	int nnz) {
//
//	int row;
//	int col;
//	int space;
//	int items;
//	int index;
//	int first_in_row;
//
//	csr_matrix_init(csr_matrix, nnz, width, height);
//
//	space = width * height;
//	items = nnz;
//	index = 0;
//
//	for (row = 0; row < height; row++) {
//		first_in_row = 1;
//
//		for (col = 0; col < width; col++) {
//			if (rand_is_nonzero(space, items)) {
//				csr_matrix->v[index] = 2 + (100 * random_double());
//				csr_matrix->ci[index] = col;
//
//				if (first_in_row) {
//					first_in_row = 0;
//
//					csr_matrix->rp[row] = index;
//				}
//
//				index++;
//				items--;
//			}
//
//			space--;
//		}
//
// 		if (first_in_row) {
//			first_in_row = 0;
//
//			csr_matrix->rp[row] = index;
//		}
//	}
//
//	csr_matrix->rp[height] = index;
//
//	return 1;
//}
//
//time_record_t csr_matrix_mm_mul(const char * matrix_a, const char *matrix_b) {
//	csr_matrix_t a;
//	csr_matrix_t b;
//	den_matrix_t c;
//	time_record_t tr;
//
//	tr.load_a = csr_matrix_load_mm(&a, matrix_a);
//	tr.load_b = csr_matrix_load_mm(&b, matrix_b);
//
//	tr.multiplication = csr_matrix_matrix_mul(&a, &b, &c);
//
//	csr_matrix_free(&a);
//	csr_matrix_free(&b);
//	dense_matrix_free(&c);
//
//	return tr;
//}
