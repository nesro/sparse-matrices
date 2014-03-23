/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
 * https://github.com/nesro/sparse-matrices
 */

#include <stdio.h>
#include <stdlib.h>

#include "dense_matrix.h"
#include "coo_matrix.h"

int dense_matrix_init(dense_matrix_t *dense_matrix, int width, int height) {

	int row;

	dense_matrix->info.nnz = -1;
	dense_matrix->info.w = width;
	dense_matrix->info.h = height;

	dense_matrix->v = malloc(dense_matrix->info.h * sizeof(datatype_t *));
	if (dense_matrix->v == NULL) {
		goto err_out;
	}

	dense_matrix->rows_block = calloc(
			dense_matrix->info.h * dense_matrix->info.w, sizeof(datatype_t));

	for (row = 0; row < dense_matrix->info.h; row++) {
		dense_matrix->v[row] = dense_matrix->rows_block
				+ row * dense_matrix->info.h;
	}

	return 1;

	err_out: /* LABEL */
	free(dense_matrix->v);
	free(dense_matrix->rows_block);

	return 0;
}

void dense_matrix_free(dense_matrix_t *dense_matrix) {
	free(dense_matrix->v);
	free(dense_matrix->rows_block);
}

/******************************************************************************/

static dense_mmm_recursion_inner(den_matrix_t *a, den_matrix_t *b,
		den_matrix_t *c, int ax, int ay, int bx, int by, int cx, int cy, int s) {

	int i, j;

	if (s == 1) {
		printf("c[%d][%d] += a[%d][%d] * b[%d][%d]\n", cy, cx, ay, ax, by, bx);
		return;
	}

	for (i = 0; i < 2; i++) {
		for (j = 0; j < 2; j++) {
			//printf("i=%d,j=%d\n", i, j);
			dense_mmm_recursion_inner(a, b, c, ax + (i * s / 2),
					ay + (j * s / 2), bx, by, cx, cy, s / 2);
		}
	}

}

void dense_matrix_matrix_mul_recursion(den_matrix_t *a, den_matrix_t *b,
		den_matrix_t *c) {

	dense_mmm_recursion_inner(a, b, c, 0, 0, 0, 0, 0, 0, a->_.w);
}

/******************************************************************************/

/**
 * Algorithm has been taken from:
 * https://fedcsis.org/proceedings/2013/pliks/19.pdf
 */
void dense_matrix_matrix_mul(dense_matrix_t *a, dense_matrix_t *b,
		dense_matrix_t *c) {

	int row;
	int col;
	int i;
	datatype_t sum;

	/*
	 * IxJ * JxK = 	IxK
	 */
	dense_matrix_init(c, a->info.h, b->info.w);

	for (row = 0; row < a->info.h; row++) {
		for (col = 0; col < b->info.w; col++) {
			sum = 0;

			for (i = 0; i < c->info.h; i++) {
				sum += a->v[row][i] * b->v[i][col];
			}

			c->v[row][col] = sum;
		}
	}
}

void dense_matrix_vector_mul(dense_matrix_t *a, vector_t *b, vector_t *c) {

	int row;
	int i;
	datatype_t sum;

	assert(a->info.w == b->size);

	vector_init(c, b->size);

	for (row = 0; row < a->info.h; row++) {
		sum = 0;

		for (i = 0; i < b->size; i++) {
			sum += a->v[row][i] * b->v[i];
		}

		c->v[row] = sum;
	}
}

void dense_to_html(dense_matrix_t *matrix, char *filename) {

	int td_size = TD_SIZE;
	FILE *fp;
	int row;
	int col;

	fp = fopen(filename, "w");
	fprintf(fp, "<html><head><style>\n"
			"*{margin:0;padding:0}\n"
			"table{border-spacing:0;border-collapse:collapse;}\n"
			"td{width:%dpx;height:%dpx;border:1px solid black;padding:2px}\n"
			"td.nonzero{background:#cc0000;}\n"
			"td.nonzero:hover{background:#ff0000;}\n"
			"td.zero{background:#cccccc;}\n"
			"td.zero:hover{background:#ffffff;}\n"
			"</style></head><body><table>\n", td_size, td_size);

	for (row = 0; row < matrix->info.h; row++) {
		fprintf(fp, "<tr>\n");

		for (col = 0; col < matrix->info.w; col++) {
			fprintf(fp,
					"<td class=\"%s\" title=\"row=%2d,col=%2d,val=" DATATYPE_FORMAT
					"\">%.0f</td>\n",
					((matrix->v[row][col] == 0) ? "zero" : "nonzero"), row, col,
					matrix->v[row][col], matrix->v[row][col]);
		}

		fprintf(fp, "</tr>\n");
	}

	fprintf(fp, "</table></body></html>\n");
	fclose(fp);
}

int dense_matrix_load_mm(dense_matrix_t *dense_matrix, const char *filename) {

	coo_matrix_t coo_matrix;

	coo_matrix_load_mm(&coo_matrix, filename);
	coo_to_dense(&coo_matrix, dense_matrix);
	coo_matrix_free(&coo_matrix);

	return 1;
}

void dense_matrix_hash_save(dense_matrix_t *matrix, const char *filename) {

	FILE *fp;

	double hash[16];
	int i;
	int row;
	int col;

	fp = fopen(filename, "w");

	if (fp == NULL) {
		fprintf(stderr, "cannot create file %s", filename);
		perror("\n");
	}

	for (i = 0; i < 16; i++) {
		hash[i] = 0;
	}

	i = 0;
	for (row = 0; row < matrix->info.h; row++) {
		for (col = 0; col < matrix->info.w; col++) {
			hash[i] += row * col * matrix->v[row][col];

			i++;

			if (i % 16 == 0) {
				i = 0;
			}
		}
	}

	for (i = 0; i < 16; i++) {
		fprintf(fp, "%lf\n", hash[i]);
	}

	fclose(fp);
}

int dense_matrix_compare(dense_matrix_t *a, dense_matrix_t *b) {

	int i;
	int j;

	if (a == NULL && b == NULL) {
		return 1;
	} else if (a == NULL || b == NULL) {
		return 0;
	}

	if (a->info.w != b->info.w || a->info.h != b->info.h) {
		return 0;
	}

	for (i = 0; i < a->info.h; i++) {
		for (j = 0; j < a->info.w; j++) {
			if (a->v[i][j] != b->v[i][j]) {
				return 0;
			}
		}
	}

	return 1;
}

