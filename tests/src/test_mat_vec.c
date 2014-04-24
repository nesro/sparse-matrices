/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
 * https://github.com/nesro/sparse-matrices
 */

#include <stdio.h>
#include <stdlib.h>

#include "../../src/virtual_matrix.h"
#include "../../src/den_matrix.h"
#include "../../src/coo_matrix.h"
#include "../../src/vector.h"

#include "../../cassertion/cassertion.h"
#include "test_utils.h"

static void run() {

	vm_t *a_den = NULL;
	vm_t *a_coo = NULL;
	vm_t *a_csr = NULL;
	vm_t *b_vec = NULL;
	vm_t *c_vec_den = NULL;
	vm_t *c_vec_coo = NULL;
	vm_t *c_vec_csr = NULL;
	const test_matrices_pair_t *tp;

	while ((tp = foreach_pair(mat_vec_pairs)) != NULL) {
		printf("%s - %s\n", tp->a.path, tp->b.path);

		vm_load_mm(&a_den, DEN, tp->a.path);
		vm_load_mm(&a_coo, COO, tp->a.path);
		vm_load_mm(&a_csr, CSR, tp->a.path);

		vm_load_mm(&b_vec, VEC, tp->b.path);
		a_den->f.mul(a_den, b_vec, &c_vec_den, NAIVE);

		CASSERTION_TIME();
		a_coo->f.mul(a_coo, b_vec, &c_vec_coo, NAIVE);
		CASSERTION(c_vec_den->f.compare(c_vec_den, c_vec_coo) == 0,
				"a=%s,b=%s, coo", tp->a.path, tp->b.path);

		CASSERTION_TIME();
		a_csr->f.mul(a_csr, b_vec, &c_vec_csr, NAIVE);
		CASSERTION(c_vec_den->f.compare(c_vec_den, c_vec_csr) == 0,
				"a=%s,b=%s, csr", tp->a.path, tp->b.path);

#if 0
		printf("---- dense:\n");
		dc->f.print(dc);
		printf("---- coo:\n");
		cc->f.print(cc);
#endif

		a_den->f.free(a_den);
		a_coo->f.free(a_coo);
		a_csr->f.free(a_csr);
		b_vec->f.free(b_vec);
		c_vec_den->f.free(c_vec_den);
		c_vec_coo->f.free(c_vec_coo);
		c_vec_csr->f.free(c_vec_csr);

		a_den = NULL;
		a_coo = NULL;
		a_csr = NULL;
		b_vec = NULL;
		c_vec_den = NULL;
		c_vec_coo = NULL;
	}
}

int main(int argc, char *argv[]) {

	CASSERTION_INIT(argc, argv);

	run();

	CASSERTION_RESULTS();

	return EXIT_SUCCESS;
}
