
#include <stdlib.h>
#include "../src/virtual_matrix.h"


void ranom_quadtree() {

	vm_t *a;
	vm_t *b;
	vm_t *c;

	vm_load_mm(&a, QDT, "file.mtx", 12);
	vm_load_mm(&b, QDT, "file.mtx", 12);
	a->f.mul(a, b, &c);
}

int main(int argc, char *argv[]) {

	return EXIT_SUCCESS;
}
