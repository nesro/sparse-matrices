/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
 * https://github.com/nesro/sparse-matrices
 */

#include <stdio.h>
#include <stdlib.h>

#include "utils.h"
#include "mmio.h"
#include "coo_matrix.h"

int coo_matrix_init(coo_matrix_t *coo_matrix, int nnz, int width, int height) {

	int i;

	coo_matrix->info.nnz = nnz;
	coo_matrix->info.w = width;
	coo_matrix->info.h = height;

	coo_matrix->values = malloc(coo_matrix->info.nnz * sizeof(datatype_t));
	if (coo_matrix->values == NULL) {
		perror("coo_matrix->values");
		exit(1);
	}

	coo_matrix->col = malloc(coo_matrix->info.nnz * sizeof(int));
	if (coo_matrix->values == NULL) {
		perror("coo_matrix->col");
		exit(1);
	}

	coo_matrix->row = malloc(coo_matrix->info.nnz * sizeof(int));
	if (coo_matrix->values == NULL) {
		perror("coo_matrix->row");
		exit(1);
	}

	for (i = 0; i < coo_matrix->info.nnz; i++) {
		coo_matrix->col[i] = 0;
		coo_matrix->row[i] = 0;
		coo_matrix->values[i] = 0;
	}

	return 1;
}

void coo_matrix_free(coo_matrix_t *coo_matrix) {
	free(coo_matrix->values);
	free(coo_matrix->col);
	free(coo_matrix->row);
}

int coo_matrix_load_mm(coo_matrix_t *coo_matrix, const char *filename) {

	MM_typecode matcode;
	FILE *f;
	int i;
	int M;
	int N;
	int nz;

	if ((f = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "File %s doesn't exists. Exiting.\n", filename);
		exit(1);
	}

	if (mm_read_banner(f, &matcode) != 0) {
		printf("Could not process Matrix Market banner.\n");
		exit(1);
	}

	if (mm_is_complex(matcode) && mm_is_matrix(matcode) && mm_is_sparse(matcode)) {
		printf("Sorry, this application does not support ");
		printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
		exit(1);
	}

	if (mm_read_mtx_crd_size(f, &M, &N, &nz) != 0) {
		perror("mm_read_mtx_crd_size");
		exit(1);
	}

	coo_matrix_init(coo_matrix, nz, N, M); /* N=width, M=height */

	for (i = 0; i < nz; i++) {

		if (LOAD_VALUE) {
			int success = -1;

			/* TODO: make better error */
			if ((success = fscanf(f, "%d %d %lg\n", &coo_matrix->row[i],
				&coo_matrix->col[i], &coo_matrix->values[i])) != 2) {
				/*printf("ounou: %d\n", success);*/

			}
		} else {
			if (fscanf(f, "%d %d \n", &coo_matrix->row[i], &coo_matrix->col[i])
				!= 2) {
				exit(1);
			}
			coo_matrix->values[i] = 1;
		}

		/*EIA adjust from 1-based to 0-based because of Matrix Market EIA*/
		coo_matrix->row[i]--;
		coo_matrix->col[i]--;
	}

	if (f != stdin)
		fclose(f);

	return 1;
}

void sort(int *col_idx, double *a, int start, int end) {

	int i;
	int j;
	int it;
	double dt;

	for (i = end - 1; i > start; i--) {
		for (j = start; j < i; j++) {
			if (col_idx[j] > col_idx[j + 1]) {
				if (a) {
					dt = a[j];
					a[j] = a[j + 1];
					a[j + 1] = dt;
				}

				it = col_idx[j];
				col_idx[j] = col_idx[j + 1];
				col_idx[j + 1] = it;
			}
		}
	}
}

//void coo_to_csr(coo_matrix_t *coo_matrix, csr_t *csr_matrix) {
//
//	int row;
//	int col;
//
//	csr_matrix_init(csr_matrix, coo_matrix->info.nnz, coo_matrix->info.w,
//		coo_matrix->info.h);
//
//	/*EIA reset row_pointers EIA*/
//	for (row = 0; row < csr_matrix->info.h; row++) {
//		csr_matrix->rp[row] = 0;
//	}
//
//	/*EIA pocet prvku v row EIA*/
//	for (row = 0; row < csr_matrix->info.nnz; row++) {
//		csr_matrix->rp[coo_matrix->row[row] + 1]++;
//	}
//
//	/*EIA posunuti prvku v poli EIA*/
//	for (row = 0; row < csr_matrix->info.h; row++) {
//		csr_matrix->rp[row + 1] += csr_matrix->rp[row];
//	}
//
//	for (col = 0; col < coo_matrix->info.nnz; col++) {
//		row = csr_matrix->rp[coo_matrix->row[col]];
//		csr_matrix->v[row] = coo_matrix->values[col];
//		csr_matrix->ci[row] = coo_matrix->col[col];
//		csr_matrix->rp[coo_matrix->row[col]]++;
//	}
//
//	/* shift back row_start */
//	for (row = csr_matrix->info.h; row > 0; row--)
//		csr_matrix->rp[row] = csr_matrix->rp[row - 1];
//
//	csr_matrix->rp[0] = 0;
//
//	/* FIXME: is this necessary? */
//	/*for (row = 0; row < csr_matrix->info.height; row++) {
//	 int i;
//	 int j;
//	 int it;
//	 double dt;
//
//	 for (i = csr_matrix->row_pointers[row + 1] - 1;
//	 i > csr_matrix->row_pointers[row]; i--) {
//	 for (j = csr_matrix->row_pointers[row]; j < i; j++) {
//	 if (csr_matrix->col_indicies[j]
//	 > csr_matrix->col_indicies[j + 1]) {
//	 if (csr_matrix->values) {
//	 dt = csr_matrix->values[j];
//	 csr_matrix->values[j] = csr_matrix->values[j + 1];
//	 csr_matrix->values[j + 1] = dt;
//	 }
//
//	 it = csr_matrix->col_indicies[j];
//	 csr_matrix->col_indicies[j] =
//	 csr_matrix->col_indicies[j + 1];
//	 csr_matrix->col_indicies[j + 1] = it;
//	 }
//	 }
//	 }
//
//	 }*/
//}

void coo_to_vector(coo_matrix_t *coo_matrix, vector_t *vector) {

	int el;

	vector_init(vector, coo_matrix->info.h);

	for (el = 0; el < coo_matrix->info.nnz; el++) {
		vector->v[coo_matrix->row[el]] = coo_matrix->values[el];
	}
}

void coo_to_dense(coo_matrix_t *coo_matrix, den_matrix_t *dense_matrix) {

	int i;

//	FIXME:
//	den_matrix_init(&dense_matrix, (den_matrix_init_t ) { coo_matrix->info.w,
//			coo_matrix->info.h, 1 });

	for (i = 0; i < coo_matrix->info.nnz; i++) {
		dense_matrix->v[coo_matrix->row[i]][coo_matrix->col[i]] =
			coo_matrix->values[i];
	}
}
