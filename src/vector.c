/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
 * https://github.com/nesro/sparse-matrices
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>

#include "utils.h"
#include "mm_load.h"
#include "vector.h"
#include "coo_matrix.h"

static vm_vmt_t vec_vmt = { /**/
(reset_t) NULL, /**/
(free_t) vec_free, /**/
(mm_load_t) NULL,/**/
(mm_save_t) vec_mm_save, /**/
(print_t) NULL, /**/
(compare_t) vec_compare, /**/
(distance_t) vec_distance, /**/
(convert_t) NULL, /**/
(mul_t) NULL, /**/
};

void vec_mm_save(vec_t *vec, const char *output) {

	FILE *f = NULL;
	int f_stdout = 0;
	int i;

	if (strcmp(output, "stdout") == 0) {
		f = stdout;
		f_stdout = 1;
	} else {
		f = fopen(output, "w");
	}

	/* FIXME: w h ? */
	fprintf(f, "%%MatrixMarket matrix coordinate real general\n"
			"%d %d %d\n", 1, vec->_.h, vec->size);

	for (i = 0; i < vec->size; i++)
		fprintf(f, "%d %d %lf\n", i + 1, 1, vec->v[i]);

	if (!f_stdout)
		fclose(f);
}

void vec_vm_init(vec_t **vec, va_list va) {

	int va_flag = 0;
	int size;

	size = va_get_int(va, 1, &va_flag);

	vec_init(vec, size);
}

void vec_init(vec_t **vec, int size) {

	*vec = malloc(sizeof(vec_t));
	assert(*vec != NULL);
	(*vec)->_.object_size += sizeof(vec_t);

	(*vec)->_.type = VEC;
	(*vec)->_.f = vec_vmt;
	(*vec)->_.h = size;

	/*
	 * Because the .mtx file is sparse, it would be possible to have some
	 * uninitialized elements.
	 */
	(*vec)->v = calloc(size, sizeof(datatype_t));
	assert((*vec)->v != NULL);
	(*vec)->_.object_size += size * sizeof(datatype_t);
}

void vec_free(vec_t *vec) {

	free(vec->v);
	free(vec);
}

void vec_from_mm(vec_t **vec, const char *file, va_list va /* unused */) {
	vec_load_mm(vec, file);
}

double vec_load_mm(vec_t **vec, const char *filename) {

	double start_time;
	double end_time;
	mm_file_t *mm_file;
	int i;

	mm_file = mm_load(filename, 0);

	/*
	 * Start timer after mm_file is loaded. We don't want to IO operations
	 * in our timer.
	 */
	start_time = omp_get_wtime();

	assert(mm_file->width == 1);

	vec_init(vec, mm_file->height);

	for (i = 0; i < mm_file->nnz; i++)
		(*vec)->v[mm_file->data[i].row] = mm_file->data[i].value;

	end_time = omp_get_wtime();
	mm_free(mm_file);
	return end_time - start_time;
}

/******************************************************************************/

/*
 * Compare two vectors.
 * If they're same, the return value is ZERO, to be in par with
 * the vec_distance function.
 */
int vec_compare(vec_t *a, vec_t *b) {

	int i;

	if (a->_.h != b->_.h)
		return 1;

	for (i = 0; i < a->_.h; i++)
		if (a->v[i] != a->v[i])
			return 1;

	return 0;
}

/*
 * Distance of two vectors.
 */
int vec_distance(vec_t *a, vec_t *b) {

	int i;
	int distance;

	if (a->size != b->size)
		return -1;

	distance = 0;
	for (i = 0; i < a->size; i++)
		if (a->v[i] != a->v[i])
			distance++;

	return distance;
}

/******************************************************************************/

//error_t vector_init(vec_t *vector, int size) {
//
//	vector->size = size;
//	vector->v = calloc(vector->size, sizeof(datatype_t));
//	if (vector->v == NULL) {
//		return ERROR_NOMEMORY;
//	}
//
//	return ERROR_NONE;
//}
//
//void vector_free(vec_t *vector) {
//	free(vector->v);
//}
//
//int vector_load_mm(vec_t *vector, const char *filename) {
//
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
//
//	return 1;
//}
//
//void vector_hash_save(vec_t *vector, const char *filename) {
//
//	FILE *fp;
//
//	double hash[16];
//	int i;
//	int row;
//
//	fp = fopen(filename, "w");
//
//	for (i = 0; i < 16; i++) {
//		hash[i] = 0;
//	}
//
//	i = 0;
//	for (row = 0; row < vector->size; row++) {
//		hash[i] += row * vector->v[row];
//
//		i++;
//
//		if (i % 16 == 0) {
//			i = 0;
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
//int vector_compare(vec_t *a, vec_t *b) {
//
//	int i;
//
//	if (a->size != b->size) {
//		return -1;
//	}
//
//	for (i = 0; i < a->size; i++) {
//		if (a->v[i] != b->v[i]) {
//			return 0;
//		}
//	}
//
//	return 1;
//}
//
//error_t vector_generate(vec_t *vector, int size) {
//
//	error_t error;
//	int i;
//
//	srand(0);
//
//	error = vector_init(vector, size);
//	if (error != ERROR_NONE) {
//		return error;
//	}
//
//	for (i = 0; i < size; i++) {
//		vector->v[i] = 2 + random_double() * 100;
//	}
//
//	return ERROR_NONE;
//}
