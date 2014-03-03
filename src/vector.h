/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
 * https://github.com/nesro/sparse-matrices
 */

#ifndef VECTOR_H_
#define VECTOR_H_

#include "utils.h"
#include "virtual_matrix.h"

typedef struct vector {
	int size;
	datatype_t *v;
} vector_t;

error_t vector_init(vector_t *vector, int size);
void vector_free(vector_t *vector);

int vector_load_mm(vector_t *vector, const char *filename);
void vector_hash_save(vector_t *vector, const char *filename);

int vector_compare(vector_t *a, vector_t *b);

error_t vector_generate(vector_t *vector, int size);

#endif /* VECTOR_H_ */
