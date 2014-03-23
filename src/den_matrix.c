/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
 * https://github.com/nesro/sparse-matrices
 */

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

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
(compare_t) NULL, /**/
(convert_t) NULL, /**/
(mul_t) mul, /**/
};

/***************************************************************************/

void den_vm_init(den_matrix_t **den, va_list va) {

	int va_flag = 0;
	int width;
	int height;
	int zero;

	zero = va_get_int(va, 1, &va_flag);
	height = va_get_int(va, -1, &va_flag);
	width = va_get_int(va, -1, &va_flag);

	den_matrix_init(den, width, height, zero);
}

void den_from_mm(den_matrix_t **den, const char *mm_filename, va_list va) {

	mm_file_t *mm_file;
	int i;

	mm_file = mm_load(mm_filename);

	den_matrix_init(den, mm_file->width, mm_file->height, 1);

	for (i = 0; i < mm_file->items; i++) {
		(*den)->v[mm_file->data[i].col - 1][mm_file->data[i].row - 1] +=
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
}

void den_matrix_free(den_matrix_t *den_matrix) {
	free(den_matrix->v);
	free(den_matrix->rows_block);
	free(den_matrix);
}

void dense_matrix_free(den_matrix_t *dense_matrix) {
	free(dense_matrix->v);
	free(dense_matrix->rows_block);
	free(dense_matrix);
}

void den_matrix_print(den_matrix_t *den) {

	int i;
	int j;

	vm_print(&den->_);

	for (i = 0; i < den->_.h; i++) {
		for (j = 0; j < den->_.w; j++) {
			printf(DPF ", ", den->v[j][i]);
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

double den_mul_strassen(const den_matrix_t *a, const den_matrix_t *b,
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

	start_time = omp_get_wtime();

	for (row = a->_.h - 1; row >= 0; row--) {
		for (col = b->_.w - 1; col >= 0; col--) {
			sum = 0;

			for (i = c->_.h - 1; i >= 0; i--) {
				sum += a->v[row][i] * b->v[i][col];
			}

			c->v[row][col] = sum;
		}
	}

	return (omp_get_wtime() - start_time);
}

/******************************************************************************/

static void den_mul_recursion_inner(const den_matrix_t *a,
		const den_matrix_t *b, den_matrix_t *c, int ax, int ay, int bx, int by,
		int cx, int cy, int s) {

	int row, col, i;

	if (s == 1) {
		printf("c[%d][%d] += a[%d][%d] * b[%d][%d]\n", cy, cx, ay, ax, by, bx);
		c->v[cy][cx] += a->v[ay][ax] * b->v[by][bx];
		return;
	}

	for (row = 0; row < 2; row++) {
		for (col = 0; col < 2; col++) {
			for (i = 0; i < 2; i++) {
				printf("row=%d,col=%d,i=%d\n", row, col, i);

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

	die("The flag \"%x\" is not valid.\n", flag);

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
