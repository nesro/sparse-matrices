/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013,2014
 * https://github.com/nesro/sparse-matrices
 */

#ifndef VECTOR_H_
#define VECTOR_H_

#include "utils.h"
#include "virtual_matrix.h"

typedef struct vector {
	vm_t _;
	int size;
	datatype_t *v;
} vec_t;

void vec_vm_init(vec_t **vec, va_list va);
void vec_init(vec_t **vec, int size);
void vec_free(vec_t *vec);
void vec_from_mm(vec_t **vec, const char *file, va_list va);
double vec_load_mm(vec_t **vec, const char *filename, int min_size);
int vec_distance(vec_t *a, vec_t *b);
int vec_compare(vec_t *a, vec_t *b);
void vec_mm_save(vec_t *vec, const char *output);

//error_t vector_init(vec_t *vector, int size);
//void vector_free(vec_t *vector);
//int vector_load_mm(vec_t *vector, const char *filename);
//void vector_hash_save(vec_t *vector, const char *filename);
//int vector_compare(vec_t *a, vec_t *b);
//error_t vector_generate(vec_t *vector, int size);

#endif /* VECTOR_H_ */
