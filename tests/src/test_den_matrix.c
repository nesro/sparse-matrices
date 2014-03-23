#include <stdio.h>
#include <stdlib.h>

#include "../../src/virtual_matrix.h"
#include "../../src/den_matrix.h"

#include "../../cassertion/cassertion.h"

static void rec() {
	vm_t *a = NULL;
	vm_t *b = NULL;
	vm_t *c1 = NULL;
	vm_t *c2 = NULL;

	vm_load_mm(&a, DEN, "./matrices/2x2_4nz_01.mtx");
	vm_load_mm(&b, DEN, "./matrices/2x2_4nz_02.mtx");

	a->f.mul(a, b, &c1, NAIVE);
	a->f.mul(a, b, &c2, RECURSIVE);

	printf("---- naive:\n");
	c1->f.print(c1);
	printf("---- recursive:\n");
	c2->f.print(c2);

	a->f.free(a);
	b->f.free(b);
	c1->f.free(c1);
	c2->f.free(c2);
}

int main(int argc, char *argv[]) {

	CASSERTION_INIT(argc, argv);

	rec();

	CASSERTION_RESULTS();

	return EXIT_SUCCESS;
}
