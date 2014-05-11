/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
 * https://github.com/nesro/sparse-matrices
 */

#include <stdio.h>
#include <stdlib.h>

#include "../../src/virtual_matrix.h"
#include "../../src/coo_matrix.h"

#include "../../cassertion/cassertion.h"
#include "test_utils.h"

/* COO * COO is not implemented */
void run() {

	vm_t *da = NULL;
	vm_t *db = NULL;
	vm_t *dc = NULL;
	vm_t *ca = NULL;
	vm_t *cb = NULL;
	vm_t *cc = NULL;

	const test_matrices_pair_t *tp;

	while ((tp = foreach_pair(tm_pairs)) != NULL) {
		printf("%s - %s\n", tp->a.path, tp->b.path);

		vm_load_mm(&da, DEN, tp->a.path);
		vm_load_mm(&db, DEN, tp->b.path);

		vm_load_mm(&ca, COO, tp->a.path);
		vm_load_mm(&cb, COO, tp->b.path);

		da->f.mul(da, db, &dc, NAIVE);

		CASSERTION_TIME();
		ca->f.mul(ca, cb, &cc, NAIVE);

		CASSERTION(dc->f.compare(dc, cc) == 0, "a=%s,b=%s, coo", tp->a.path,
				tp->b.path);

#if 0
		printf("---- dense:\n");
		dc->f.print(dc);
		printf("---- coo:\n");
		cc->f.print(cc);
#endif

		da->f.free(da);
		db->f.free(db);
		dc->f.free(dc);
		ca->f.free(ca);
		cb->f.free(cb);
		cc->f.free(cc);

		da = NULL;
		db = NULL;
		dc = NULL;
		ca = NULL;
		cb = NULL;
		cc = NULL;
	}
}

int main(int argc, char *argv[]) {

	CASSERTION_INIT(argc, argv);

	//run();

	CASSERTION_RESULTS();

	return EXIT_SUCCESS;
}
