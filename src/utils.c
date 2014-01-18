/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
 * https://github.com/nesro/sparse-matrices
 */

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"


double random_double(void) {
	return ((double) rand() / (double) RAND_MAX);
}

int possibility(double n) {
	if (n == 1)
		return 1;

	return (random_double() < n);
}

int rand_is_nonzero(int space, int items) {
	assert(space > 0);
	assert(space >= items);

	return possibility((double) items / (double) space);
}

int is_power_of_two(unsigned int x) {
	assert(x > 1);
	return ((x != 0) && ((x & (~x + 1)) == x));
}

/******************************************************************************/

#define MEMORY_USAGE 1

#if MEMORY_USAGE == 1

#include <math.h>

static int allocated_g = 0;

#define malloc(a) _malloc(a)
#define calloc(a,b) _calloc((a), (b))
#define realloc(a,b) _realloc((a),(b),(sqrt(b)))

void *_malloc(size_t size) {
	void *mem_ptr;
	allocated_g += size;
	mem_ptr = malloc(size);
	if (mem_ptr == NULL) {
		perror("malloc");
		exit(EXIT_FAILURE);
	}
	return mem_ptr;
}

void *_calloc(size_t count, size_t element_size) {
	void *mem_ptr;
	allocated_g += count * element_size;
	mem_ptr = calloc(count, element_size);
	if (mem_ptr == NULL) {
		perror("calloc");
		exit(EXIT_FAILURE);
	}
	return mem_ptr;
}

void *_realloc(void *mem_ptr, size_t new_size, size_t old_size) {
	void *new_mem_ptr;
	allocated_g += new_size - old_size;
	new_mem_ptr = realloc(mem_ptr, new_size);
	if (new_mem_ptr == NULL) {
		perror("realloc");
		exit(EXIT_FAILURE);
	}
	return new_mem_ptr;
}

#endif
