/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
 * https://github.com/nesro/sparse-matrices
 */

#include <stdio.h>
#include <stdlib.h>

#include "generator.h"

#include "den_matrix.h"

void generate_dense(den_matrix_t *dense, int n, int density) {

	int i;
	int j;

//	FIXME:
//	den_matrix_init(&dense, (den_matrix_init_t ) { n, n, 0 });

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			dense->v[i][j] = possibility(density / 100) == 1 ? 1 : 0;

}

//void generate_mm(const char* filename, int width, int height, double density,
//	int srand_seed) {
//
//	FILE *f;
//	int nnz;
//	int i;
//	int last_x;
//	int last_y;
//
//	if (srand < 0) {
//		fprintf(stderr, "ERROR: srand_seed must be >= 0.", filename);
//		exit(1);
//	}
//
//	srand(srand_seed);
//
//	f = fopen(filename, "w");
//	if (f == NULL) {
//		fprintf(stderr, "ERROR: Cannot open the file %s in the w mode.",
//			filename);
//		exit(1);
//	}
//
//	nnz = (width * height) / density;
//
//	fprintf(f, "%%MatrixMarket matrix coordinate real general\n");
//	fprintf(f, "%d %d %d\n", width, height, nnz);
//
//	last_x = 0;
//	last_y = 0;
//
//	for (i = 0; i < nnz; i++) {
//
//	}
//
//	fclose(f);
//}
