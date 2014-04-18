/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
 * https://github.com/nesro/sparse-matrices
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <math.h>

#include "utils.h"
#include "mm_load.h"
#include "den_matrix.h"
#include "coo_matrix.h"

/***************************************************************************/

static vm_vmt_t den_vmt = { /**/
(reset_t) NULL, /**/
(free_t) den_matrix_free, /**/
(mm_load_t) NULL,/**/
(mm_save_t) NULL, /**/
(print_t) den_matrix_print, /**/
(compare_t) den_compare, /**/
(distance_t) den_distance, /**/
(convert_t) NULL, /**/
(mul_t) mul, /**/
};

int den_compare(den_matrix_t *a, vm_t *b) {

	int i;
	int j;
	vm_t *tmp = NULL;

	if (a->_.w != b->w || a->_.h != b->h)
		return 2;

	if (b->type != DEN)
		tmp = b->f.convert(b, DEN);
	else
		tmp = b;

	for (i = 0; i < a->_.h; i++)
		for (j = 0; j < a->_.w; j++)
			if (a->v[i][j] != ((den_matrix_t *) tmp)->v[i][j])
				return 1;

	if (b != tmp && tmp != NULL)
		tmp->f.free(tmp);

	return 0;
}

double den_distance(den_matrix_t *a, den_matrix_t *b) {

	double ret = 0.;
	int i;
	int j;

	if (a->_.w != b->_.w || a->_.h != a->_.h)
		return INFINITY;

	for (i = 0; i < a->_.h; i++) {
		for (j = 0; j < a->_.w; j++) {
			ret += fabs(a->v[i][j] - b->v[i][j]);
		}
	}

	return ret;
}

int den_count_nnz(const den_matrix_t *a) {

	int nnz = 0;
	int i;
	int j;

	for (i = 0; i < a->_.h; i++)
		for (j = 0; j < a->_.w; j++)
			if (a->v[i][j] != ((datatype_t) 0))
				nnz++;

	return nnz;
}

/***************************************************************************/

void den_vm_init(den_matrix_t **den, va_list va) {

	int va_flag = 0;
	int width;
	int height;
	int zero;

	zero = va_get_int(va, 1, &va_flag);
	height = va_get_int(va, -1, &va_flag);
	width = va_get_int(va, -1, &va_flag);

//	printf("width=%d, height=%d, zero=%d\n", width, height, zero);

	den_matrix_init(den, width, height, zero);
}

void den_from_mm(den_matrix_t **den, const char *mm_filename,
		va_list va /* unused */) {

	mm_file_t *mm_file;
	int i;

	mm_file = mm_load(mm_filename);

	den_matrix_init(den, mm_file->width, mm_file->height, 1);

	for (i = 0; i < mm_file->nnz; i++) {
		(*den)->v[mm_file->data[i].row][mm_file->data[i].col] +=
				mm_file->data[i].value;
	}

	mm_free(mm_file);
}

void den_matrix_init(den_matrix_t **den, int width, int height, int zero) {

	int row;

	if (den != NULL && *den != NULL)
		den_matrix_free(*den);

	*den = malloc(sizeof(den_matrix_t));

	//_debug("w=%d h=%d z=%d\n", width, height, zero);

	(*den)->_.type = DEN;
	(*den)->_.f = den_vmt;
	(*den)->_.w = width;
	(*den)->_.h = height;

	(*den)->v = malloc(height * sizeof(datatype_t *));
	if ((*den)->v == NULL)
		die("dense matrix memory error");

	if (zero)
		(*den)->rows_block = calloc(width * height, sizeof(datatype_t));
	else
		(*den)->rows_block = malloc(width * height * sizeof(datatype_t));

	for (row = 0; row < height; row++)
		(*den)->v[row] = (*den)->rows_block + row * height;

	(*den)->strassen_block_treshold = 2;
}

void den_matrix_free(den_matrix_t *den_matrix) {
	free(den_matrix->v);
	free(den_matrix->rows_block);
	free(den_matrix);
}

/******************************************************************************/

void den_matrix_print(den_matrix_t *den) {

	int i;
	int j;

	vm_print(&den->_);

	for (i = 0; i < den->_.h; i++) {
		for (j = 0; j < den->_.w; j++) {
			printf(DPF ", ", den->v[i][j]);
		}

		printf("\n");
	}

	fflush(stdout);
}

/*double den_mul_naive(const den_matrix_t *a, const den_matrix_t *b,
 den_matrix_t **_c, char flag) {

 int row;
 int col;
 int i;
 datatype_t sum;
 double start_time;
 den_matrix_t *c;

 *
 * Because we are setting every element at once, we don't need
 * initialize matrix to zero. The matrix items are uninitialized
 * for now.
 *
 if (*_c == NULL) {
 den_matrix_init(_c, b->_.w, a->_.h, 0);
 } else if ((*_c)->_.w != b->_.w || (*_c)->_.h != a->_.h) {
 free(*_c);
 c = NULL;
 den_matrix_init(_c, b->_.w, a->_.h, 0);
 }

 c = *_c;

 start_time = omp_get_wtime();

 if (flag & UNROLLED) {
 for (row = 0; row < a->_.h; row++) {
 for (col = 0; col < b->_.w; col++) {
 sum = 0;

 for (i = 0; i < c->_.h; i += 8) {
 sum += a->v[row][i + 0] * b->v[i + 0][col];
 sum += a->v[row][i + 1] * b->v[i + 1][col];
 sum += a->v[row][i + 2] * b->v[i + 2][col];
 sum += a->v[row][i + 3] * b->v[i + 3][col];
 sum += a->v[row][i + 4] * b->v[i + 4][col];
 sum += a->v[row][i + 5] * b->v[i + 5][col];
 sum += a->v[row][i + 6] * b->v[i + 6][col];
 sum += a->v[row][i + 7] * b->v[i + 7][col];
 }
 for (i -= 8; i < c->_.h; i++) {
 sum += a->v[row][i] * b->v[i][col];
 }

 c->v[row][col] = sum;
 }
 }
 } else {
 for (row = 0; row < a->_.h; row++) {
 for (col = 0; col < b->_.w; col++) {
 sum = 0;

 for (i = 0; i < c->_.h; i++) {
 sum += a->v[row][i] * b->v[i][col];
 }

 c->v[row][col] = sum;
 }
 }
 }

 return (omp_get_wtime() - start_time);
 }*/

double den_mul_unrolled_parallel(const den_matrix_t *a, const den_matrix_t *b,
		den_matrix_t *c) {
	return 0.;
}

double den_mul_unrolled(const den_matrix_t *a, const den_matrix_t *b,
		den_matrix_t *c) {

	int row;
	int col;
	int i;
	datatype_t sum;
	double start_time;

	start_time = omp_get_wtime();

	for (row = 0; row < a->_.h; row++) {
		for (col = 0; col < b->_.w; col++) {
			sum = 0;

			for (i = 0; i < c->_.h - UNROLL; i += UNROLL) {
				sum += a->v[row][i + 0] * b->v[i + 0][col];
				sum += a->v[row][i + 1] * b->v[i + 1][col];
				sum += a->v[row][i + 2] * b->v[i + 2][col];
				sum += a->v[row][i + 3] * b->v[i + 3][col];
				sum += a->v[row][i + 4] * b->v[i + 4][col];
				sum += a->v[row][i + 5] * b->v[i + 5][col];
				sum += a->v[row][i + 6] * b->v[i + 6][col];
				sum += a->v[row][i + 7] * b->v[i + 7][col];
			}
			for (; i < c->_.h; i++) {
				sum += a->v[row][i] * b->v[i][col];
			}

			c->v[row][col] = sum;
		}
	}

	return (omp_get_wtime() - start_time);
}

double den_mul_parallel(const den_matrix_t *a, const den_matrix_t *b,
		den_matrix_t *c) {
	return 0.;
}

double den_mul_strassen_unrolled_parallel(const den_matrix_t *a,
		const den_matrix_t *b, den_matrix_t *c) {
	return 0.;
}

double den_mul_strassen_unrolled(const den_matrix_t *a, const den_matrix_t *b,
		den_matrix_t *c) {
	return 0.;
}

double den_mul_strassen_parallel(const den_matrix_t *a, const den_matrix_t *b,
		den_matrix_t *c) {
	return 0.;
}

double den_mul_naive(const den_matrix_t *a, const den_matrix_t *b,
		den_matrix_t *c) {

	int row;
	int col;
	int i;
	datatype_t sum;
	double start_time;

#if 0
	double zero=0;
	double nonzero=0;
#endif /* 0 */

	start_time = omp_get_wtime();

	for (row = a->_.h - 1; row >= 0; row--) {
		for (col = b->_.w - 1; col >= 0; col--) {
			sum = 0;

			for (i = c->_.h - 1; i >= 0; i--) {

#if 0
				if(a->v[row][i] == 0. || b->v[i][col] == 0.)
				zero++;
				else
				nonzero++;
#endif /* 0 */

				sum += a->v[row][i] * b->v[i][col];
			}

			c->v[row][col] = sum;
		}
	}

#if 0
	printf("zero=%lf nonzero=%lf\n", zero, nonzero);
#endif /* 0 */

	return (omp_get_wtime() - start_time);
}

/******************************************************************************/

static void den_mul_recursion_inner(const den_matrix_t *a,
		const den_matrix_t *b, den_matrix_t *c, int ax, int ay, int bx, int by,
		int cx, int cy, int s) {

	int row, col, i;

	if (s == 1) {
		//printf("c[%d][%d] += a[%d][%d] * b[%d][%d]\n", cy, cx, ay, ax, by, bx);
		c->v[cy][cx] += a->v[ay][ax] * b->v[by][bx];
		return;
	}

	for (row = 0; row < 2; row++) {
		for (col = 0; col < 2; col++) {
			for (i = 0; i < 2; i++) {
				//printf("row=%d,col=%d,i=%d\n", row, col, i);

				den_mul_recursion_inner(a, b, c, /**/
				ax + (i * (s / 2)), ay + (row * (s / 2)), /**/
				bx + (col * (s / 2)), by + (i * (s / 2)), /**/
				cx + (col * (s / 2)), cy + (row * (s / 2)), /**/
				s / 2);
			}
		}
	}

}

double den_mul_recursion(const den_matrix_t *a, const den_matrix_t *b,
		den_matrix_t *c) {
	den_mul_recursion_inner(a, b, c, 0, 0, 0, 0, 0, 0, a->_.w);
	return 0.;
}

/********************************************************** offset operations */

void den_offset_add(const den_matrix_t *a, const den_matrix_t *b,
		den_matrix_t *c, int ar, int ac, int br, int bc, int cr, int cc, int s) {

	int row;
	int col;

	for (row = 0; row < s; row++) {
		for (col = 0; col < s; col++) {
			c->v[cr + row][cc + col] = a->v[ar + row][ac + col]
					+ b->v[br + row][bc + col];
		}
	}
}

void den_offset_addto(const den_matrix_t *a, den_matrix_t *c, int ar, int ac,
		int cr, int cc, int s) {

	int row;
	int col;

	for (row = 0; row < s; row++) {
		for (col = 0; col < s; col++) {
			c->v[cr + row][cc + col] += a->v[ar + row][ac + col];
		}
	}
}

void den_offset_sub(const den_matrix_t *a, const den_matrix_t *b,
		den_matrix_t *c, int ar, int ac, int br, int bc, int cr, int cc, int s) {

	int row;
	int col;

	for (row = 0; row < s; row++) {
		for (col = 0; col < s; col++) {
			c->v[cr + row][cc + col] = a->v[ar + row][ac + col]
					- +b->v[br + row][bc + col];
		}
	}
}

/******************************************************************************/

/*

 An example of the Strassen algorithm in use:

 {{1,2},
 {3,4}}
 *
 {{5,6},
 {7,8}}
 =
 {{19,22},
 {43,50}}

 A11 = 1
 A12 = 2
 A21 = 3
 A22 = 4

 B11 = 5
 B12 = 6
 B21 = 7
 B22 = 8

 M1 = (A11 + A22) * (B11 + B22) = (1 + 4) * (5 + 8) = 5 * 13 = 65
 M2 = (A21 + A22) * B11 = (3 + 4) * 5 = 7 * 5 = 35
 M3 = A11 * (B12 - B22) = 1 * (6 - 8) = 1 * -2 = -2
 M4 = A22 * (B21 - B11) = 4 * (7 - 5) = 4 * 2 = 8
 M5 = (A11 + A12) * B22 = (1 + 2) * 8 = 3 * 8 = 24
 M6 = (A21 - A11) * (B11 + B12) = (3 - 1) * (5 + 6) = 2 * 11 = 22
 M7 = (A12 - A22) * (B21 + B22) = (2 - 4) * (7 + 8) = -2 * 15 = -30

 C11 = M1 + M4 - M5 + M7 = 65 + 8 - 24 + (-30) = 19
 C12 = M3 + M5 = -2 + 24 = 22
 C21 = M2 + M4 = 35 + 8 = 43
 C22 = M1 - M2 + M3 + M6 = 65 - 35 + (-2) + 22 = 50

 */

static void den_mul_strassen_inner(const den_matrix_t *a, const den_matrix_t *b,
		den_matrix_t *c, int ay, int ax, int by, int bx, int cy, int cx, int s) {

	int i;
	int h; /* half */
	den_matrix_t *m[9];

#if 1
	int block_size = a->strassen_block_treshold;

	if (s == block_size) {
		int j, k;
		double sum;

		for (j = block_size - 1; j >= 0; j--) {
			for (k = block_size - 1; k >= 0; k--) {
				sum = 0;

				for (i = block_size - 1; i >= 0; i--) {
					sum += a->v[ay + j][ax + i] * b->v[by + i][bx + k];
				}

				c->v[cy + j][cx + k] = sum;
			}
		}

		return;
	}
#else
	if (s == 1) {
		c->v[cy][cx] = a->v[ay][ax] * b->v[by][bx];
		return;
	}
#endif

	h = s / 2;
	for (i = 0; i < 9; i++) {
		m[i] = NULL;
		vm_create((vm_t **) &m[i], DEN, 1, h, h);
	}

	/*
	 (1 ) [0]M1 = (A11TL + A22BR) * (B11TL + B22BR)
	 (2 ) [1]M2 = (A21BL + A22BR) * B11TL
	 (3 ) [2]M3 = A11TL * (B12TR - B22BR)
	 (4 ) [3]M4 = A22BR * (B21BL - B11TL)
	 (5 ) [4]M5 = (A11TL + A12TR) * B22BR
	 (6 ) [5]M6 = (A21BL - A11TL) * (B11TL + B12TR)
	 (7 ) [6]M7 = (A12TR - A22BR) * (B21BL + B22BR)

	 (8 ) C11TL = M1 + M4 - M5 + M7
	 (9 ) C12TR = M3 + M5
	 (10) C21BL = M2 + M4
	 (11) C22BR = M1 - M2 + M3 + M6
	 */

	/*M1*/
	den_offset_add(a, a, m[7], ay, ax, ay + h, ax + h, 0, 0, h);
	den_offset_add(b, b, m[8], by, bx, by + h, bx + h, 0, 0, h);
	den_mul_strassen_inner(m[7], m[8], m[0], 0, 0, 0, 0, 0, 0, h);

	/*M2*/
	den_offset_add(a, a, m[7], ay + h, ax, ay + h, ax + h, 0, 0, h);
	den_mul_strassen_inner(m[7], b, m[1], 0, 0, bx, by, 0, 0, h);

	/*M3 */
	den_offset_sub(b, b, m[7], by, bx + h, by + h, bx + h, 0, 0, h);
	den_mul_strassen_inner(a, m[7], m[2], ay, ax, 0, 0, 0, 0, h);

	/* M4 */
	den_offset_sub(b, b, m[7], by + h, bx, by, bx, 0, 0, h);
	den_mul_strassen_inner(a, m[7], m[3], ay + h, ax + h, 0, 0, 0, 0, h);

	/* M5 */
	den_offset_add(a, a, m[7], ay, ax, ay, ax + h, 0, 0, h);
	den_mul_strassen_inner(m[7], b, m[4], 0, 0, by + h, bx + h, 0, 0, h);

	/* M6 */
	den_offset_sub(a, a, m[7], ay + h, ax, ay, ax, 0, 0, h);
	den_offset_add(b, b, m[8], by, bx, by, bx + h, 0, 0, h);
	den_mul_strassen_inner(m[7], m[8], m[5], 0, 0, 0, 0, 0, 0, h);

	/* M7 */
	den_offset_sub(a, a, m[7], ay, ax + h, ay + h, ax + h, 0, 0, h);
	den_offset_add(b, b, m[8], by + h, bx, by + h, bx + h, 0, 0, h);
	den_mul_strassen_inner(m[7], m[8], m[6], 0, 0, 0, 0, 0, 0, h);

	/* c1,1*/
	den_offset_add(m[0], m[3], m[7], 0, 0, 0, 0, 0, 0, h);
	den_offset_sub(m[7], m[4], m[7], 0, 0, 0, 0, 0, 0, h);
	den_offset_add(m[7], m[6], c, 0, 0, 0, 0, cy, cx, h);

	/* c1,2*/
	den_offset_add(m[2], m[4], c, 0, 0, 0, 0, cy, cx + h, h);

	/* c2,1*/
	den_offset_add(m[1], m[3], c, 0, 0, 0, 0, cy + h, cx, h);

	/* c2,2*/
	den_offset_sub(m[0], m[1], m[7], 0, 0, 0, 0, 0, 0, h);
	den_offset_add(m[7], m[2], m[7], 0, 0, 0, 0, 0, 0, h);
	den_offset_add(m[7], m[5], c, 0, 0, 0, 0, cy + h, cx + h, h);

	for (i = 0; i < 9; i++) {
		((vm_t *) m[i])->f.free(((vm_t *) m[i]));
	}
}

double den_mul_strassen(const den_matrix_t *a, const den_matrix_t *b,
		den_matrix_t *c) {
	den_mul_strassen_inner(a, b, c, 0, 0, 0, 0, 0, 0, a->_.h);
	return 0.;
}

/******************************************************************************/

double mul(const den_matrix_t *a, const den_matrix_t *b, den_matrix_t **c,
		char flag) {

	int init_zeros = 1;

	/*
	 * We have to initialize C matrix. For naive method we can compute every
	 * element alone. So we don't need to initialize matrix with zeros for
	 * naive algorithm. But for Strassen's we need to initialize matrix to
	 * zeros because we are just adding up to elements and thus we need start
	 * on zero.
	 */
	if (flag & NAIVE)
		init_zeros = 0;
	else if ((flag & STRASSEN) || (flag & RECURSIVE))
		init_zeros = 1;

	if (*c == NULL) {
		den_matrix_init(c, b->_.w, a->_.h, init_zeros);
	} else if ((*c)->_.w != b->_.w || (*c)->_.h != a->_.h) {
		free(*c);
		c = NULL;
		den_matrix_init(c, b->_.w, a->_.h, init_zeros);
	}

	if (flag & NAIVE)
		return den_mul_naive(a, b, *c);
	else if (flag & STRASSEN)
		return den_mul_strassen(a, b, *c);
	else if (flag & RECURSIVE)
		return den_mul_recursion(a, b, *c);

	fdie("The flag \"%x\" is not valid.\n", flag);

	/*NOTREACHABE*/
	return -1;
}

//
///**
// * Algorithm has been taken from:
// * https://fedcsis.org/proceedings/2013/pliks/19.pdf
// */
//void dense_matrix_matrix_mul(den_matrix_t *a, den_matrix_t *b, den_matrix_t *c) {
//
//	int row;
//	int col;
//	int i;
//	datatype_t sum;
//
//	/*
//	 * IxJ * JxK = 	IxK
//	 */
//	den_matrix_init(c, a->info.h, b->info.w);
//
//	for (row = 0; row < a->info.h; row++) {
//		for (col = 0; col < b->info.w; col++) {
//			sum = 0;
//
//			for (i = 0; i < c->info.h; i++) {
//				sum += a->v[row][i] * b->v[i][col];
//			}
//
//			c->v[row][col] = sum;
//		}
//	}
//}
//
//void dense_matrix_vector_mul(den_matrix_t *a, vector_t *b, vector_t *c) {
//
//	int row;
//	int i;
//	datatype_t sum;
//
//	assert(a->info.w == b->size);
//
//	vector_init(c, b->size);
//
//	for (row = 0; row < a->info.h; row++) {
//		sum = 0;
//
//		for (i = 0; i < b->size; i++) {
//			sum += a->v[row][i] * b->v[i];
//		}
//
//		c->v[row] = sum;
//	}
//}
//
//void dense_to_html(den_matrix_t *matrix, char *filename) {
//
//	int td_size = TD_SIZE;
//	FILE *fp;
//	int row;
//	int col;
//
//	fp = fopen(filename, "w");
//	fprintf(fp, "<html><head><style>\n"
//		"*{margin:0;padding:0}\n"
//		"table{border-spacing:0;border-collapse:collapse;}\n"
//		"td{width:%dpx;height:%dpx;border:1px solid black;padding:2px}\n"
//		"td.nonzero{background:#cc0000;}\n"
//		"td.nonzero:hover{background:#ff0000;}\n"
//		"td.zero{background:#cccccc;}\n"
//		"td.zero:hover{background:#ffffff;}\n"
//		"</style></head><body><table>\n", td_size, td_size);
//
//	for (row = 0; row < matrix->info.h; row++) {
//		fprintf(fp, "<tr>\n");
//
//		for (col = 0; col < matrix->info.w; col++) {
//			fprintf(fp,
//				"<td class=\"%s\" title=\"row=%2d,col=%2d,val=" DATATYPE_FORMAT
//				"\">%.0f</td>\n",
//				((matrix->v[row][col] == 0) ? "zero" : "nonzero"), row, col,
//				matrix->v[row][col], matrix->v[row][col]);
//		}
//
//		fprintf(fp, "</tr>\n");
//	}
//
//	fprintf(fp, "</table></body></html>\n");
//	fclose(fp);
//}
//
//int dense_matrix_load_mm(den_matrix_t *dense_matrix, const char *filename) {
//
//	coo_matrix_t coo_matrix;
//
//	coo_matrix_load_mm(&coo_matrix, filename);
//	coo_to_dense(&coo_matrix, dense_matrix);
//	coo_matrix_free(&coo_matrix);
//
//	return 1;
//}
//
//void dense_matrix_hash_save(den_matrix_t *matrix, const char *filename) {
//
//	FILE *fp;
//
//	double hash[16];
//	int i;
//	int row;
//	int col;
//
//	fp = fopen(filename, "w");
//
//	if (fp == NULL) {
//		fprintf(stderr, "cannot create file %s", filename);
//		perror("\n");
//	}
//
//	for (i = 0; i < 16; i++) {
//		hash[i] = 0;
//	}
//
//	i = 0;
//	for (row = 0; row < matrix->info.h; row++) {
//		for (col = 0; col < matrix->info.w; col++) {
//			hash[i] += row * col * matrix->v[row][col];
//
//			i++;
//
//			if (i % 16 == 0) {
//				i = 0;
//			}
//		}
//	}
//
//	for (i = 0; i < 16; i++) {
//		fprintf(fp, "%lf\n", hash[i]);
//	}
//
//	fclose(fp);
//}
//
//int dense_matrix_compare(den_matrix_t *a, den_matrix_t *b) {
//
//	int i;
//	int j;
//
//	if (a == NULL && b == NULL) {
//		return 1;
//	} else if (a == NULL || b == NULL) {
//		return 0;
//	}
//
//	if (a->info.w != b->info.w || a->info.h != b->info.h) {
//		return 0;
//	}
//
//	for (i = 0; i < a->info.h; i++) {
//		for (j = 0; j < a->info.w; j++) {
//			if (a->v[i][j] != b->v[i][j]) {
//				return 0;
//			}
//		}
//	}
//
//	return 1;
//}
//
