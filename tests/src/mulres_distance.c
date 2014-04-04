/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
 * https://github.com/nesro/sparse-matrices
 */

#include <stdio.h>
#include <stdlib.h>

#define STRASSEN_STABILITY_TEST

#include "../../src/virtual_matrix.h"
#include "../../src/den_matrix.h"

int main(int argc, char *argv[]) {

	vm_t *a = NULL;
	vm_t *b = NULL;
	vm_t *c_def = NULL;
	vm_t *c_str = NULL;

	vm_load_mm(&a, DEN, argv[1]);
	vm_load_mm(&b, DEN, argv[2]);

	/* todo: check for ^2 */
	((den_matrix_t*) a)->strassen_block_treshold = atoi(argv[3]);

	a->f.mul(a, b, &c_def, NAIVE);
	a->f.mul(a, b, &c_str, STRASSEN);

//	printf("---- naive:\n");
//			c_def->f.print(c_def);
//			printf("---- recursive:\n");
//			c_str->f.print(c_str);

	printf("%lf\n", c_def->f.distance(c_def, c_str));

	a->f.free(a);
	b->f.free(b);
	c_def->f.free(c_def);
	c_str->f.free(c_str);

	return EXIT_SUCCESS;
}
