/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
 * https://github.com/nesro/sparse-matrices
 */

#include "stdlib.h"
#include "vector.h"
#include "utils.h"

#ifndef VIRTUAL_MATRIX_H_
#define VIRTUAL_MATRIX_H_

#ifdef _DIFF_TRESHOLD
#define DIFF_TRESHOLD _DIFF_TRESHOLD
#else /* _DIFF_TRESHOLD */
#define DIFF_TRESHOLD 0.001
#endif /* _DIFF_TRESHOLD */

#ifdef _PRECISION
#	if _PRECISION == 1
#		define DPF "%f"
		typedef float datatype_t;
#	elif _PRECISION == 2
#		define DPF "%lf"
		typedef double datatype_t;
#	elif _PRECISION == 3
#		define DPF "%Lf"
		typedef long double datatype_t;
#	else
#		define DPF "%lf"
		typedef double datatype_t;
#	endif
#else /* _PRESICION */
#		define DPF "%lf"
		typedef double datatype_t;
#endif /* _PRESICION */


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
	size_t object_size;
};

/***************************************************************************/

void vm_create(vm_t **, vm_type_t, ...);
void vm_load_mm(vm_t **, vm_type_t, const char *, ...);
void vm_print(vm_t *);

int vm_has_blocks(vm_type_t type);

/***************************************************************************/

#endif /* VIRTUAL_MATRIX_H_ */
