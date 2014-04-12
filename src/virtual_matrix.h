/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
 * https://github.com/nesro/sparse-matrices
 */

#include "vector.h"
#include "utils.h"

#ifndef VIRTUAL_MATRIX_H_
#define VIRTUAL_MATRIX_H_

#define DPF "%lf" /* datatype_t printf format */
typedef double datatype_t;

#define UNROLL 8

/***************************************************************************/

typedef enum vm_type {
	UNKNOWN, /**/
	VEC, /* vector */
	DEN, /* dense */
	COO, /* coordinate */
	CSR, /* compressed sparse rows */
	BSR, /* block sparse row*/
	QDT, /* quadtree */
	KAT, /* k-ary tree matrix */
}vm_type_t;

enum {
	NAIVE = 0x01, /**/
	UNROLLED = 0x02, /**/
	PARALLEL = 0x04, /**/
	STRASSEN = 0x08, /**/
	RECURSIVE = 0x10, /**/
};

/***************************************************************************/

typedef struct vm vm_t;

/***************************************************************************/

typedef void (*reset_t)(vm_t *);
typedef void (*free_t)(vm_t *);

typedef void (*mm_load_t)(vm_t *, const char *filename);
typedef void (*mm_save_t)(vm_t *, const char *filename);

typedef void (*print_t)(vm_t *);
typedef int (*compare_t)(vm_t *, vm_t *);
typedef double (*distance_t)(vm_t *, vm_t *);

typedef vm_t *(*convert_t)(vm_t *, vm_type_t);
typedef double (*mul_t)(const vm_t *, const vm_t *, vm_t **, char);

/***************************************************************************/

typedef struct vm_vmt {
	reset_t reset;
	free_t free;
	mm_load_t mm_load;
	mm_save_t mm_save;
	print_t print;
	compare_t compare;
	distance_t distance;
	convert_t convert;
	mul_t mul;
}vm_vmt_t;

struct vm {
	vm_type_t type;
	vm_vmt_t f;
	int w;
	int h;
	int nnz;
	char filename[30];
};

typedef enum action {
	MUL, /**/
}action_t;

typedef enum setting_tokens {
	MATRIX_A, /**/
	MATRIX_B, /**/
	LEAF_SIZE, /**/
}setting_tokens_t;

typedef union setting_values {
	int i;
	float f;
	double d;
	char *s;
}setting_values_t;

/***************************************************************************/

void vm_create(vm_t **, vm_type_t, ...);
void vm_load_mm(vm_t **, vm_type_t, const char *, ...);
void vm_print(vm_t *);

void vm_exec(action_t action, vm_type_t type_a, vm_type_t type_b,
	const char *file_a, const char *file_b, vm_t **c, char flag, ...);

/***************************************************************************/

#endif /* VIRTUAL_MATRIX_H_ */
