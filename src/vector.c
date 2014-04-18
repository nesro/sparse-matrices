/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
 * https://github.com/nesro/sparse-matrices
 */

#include <stdio.h>
#include <stdlib.h>

#include "utils.h"
#include "vector.h"
#include "coo_matrix.h"

error_t vector_init(vector_t *vector, int size) {

	vector->size = size;
	vector->v = calloc(vector->size, sizeof(datatype_t));
	if (vector->v == NULL) {
		return ERROR_NOMEMORY;
	}

	return ERROR_NONE;
}

void vector_free(vector_t *vector) {
	free(vector->v);
}

int vector_load_mm(vector_t *vector, const char *filename) {

//	coo_matrix_t *coo_matrix = NULL;
//
//	coo_from_mm(&coo_matrix, filename, 0);
//
//	if (coo_matrix->_.w != 1) {
//		fprintf(stderr,
//			"vector_load_mm failed: loaded matrix %s is not width 1 element",
//			filename);
//		coo_matrix_free(&coo_matrix);
//		return 0;
//	}
//
//	coo_(&coo_matrix, vector);
//	coo_free(coo_matrix);

	return 1;
}

void vector_hash_save(vector_t *vector, const char *filename) {

	FILE *fp;

	double hash[16];
	int i;
	int row;

	fp = fopen(filename, "w");

	for (i = 0; i < 16; i++) {
		hash[i] = 0;
	}

	i = 0;
	for (row = 0; row < vector->size; row++) {
		hash[i] += row * vector->v[row];

		i++;

		if (i % 16 == 0) {
			i = 0;
		}
	}

	for (i = 0; i < 16; i++) {
		fprintf(fp, "%lf\n", hash[i]);
	}

	fclose(fp);
}

int vector_compare(vector_t *a, vector_t *b) {

	int i;

	if (a->size != b->size) {
		return -1;
	}

	for (i = 0; i < a->size; i++) {
		if (a->v[i] != b->v[i]) {
			return 0;
		}
	}

	return 1;
}

error_t vector_generate(vector_t *vector, int size) {

	error_t error;
	int i;

	srand(0);

	error = vector_init(vector, size);
	if (error != ERROR_NONE) {
		return error;
	}

	for (i=0;i<size;i++) {
		vector->v[i] = 2 + random_double() * 100;
	}

	return ERROR_NONE;
}



