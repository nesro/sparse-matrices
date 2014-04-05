/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
 * https://github.com/nesro/sparse-matrices
 */

#include <stdio.h>
#include <stdlib.h>

#include "../../src/virtual_matrix.h"
#include "../../src/qdt_matrix.h"

#include "../../cassertion/cassertion.h"
#include "test_utils.h"

static void run() {

	vm_t *da = NULL;
	vm_t *db = NULL;
	vm_t *dc = NULL;
	vm_t *qa = NULL;
	vm_t *qb = NULL;
	vm_t *qc = NULL;

	const test_matrix_t *tm;
	const test_matrices_pair_t *tp;

	while ((tm = foreach_matrix(tm_small)) != NULL) {
		printf("%s\n", tm->path);
	}

	while ((tp = foreach_pair(tm_pairs)) != NULL) {
		printf("%s - %s\n", tp->a.path, tp->b.path);

		vm_load_mm(&da, DEN, tp->a.path);
		vm_load_mm(&db, DEN, tp->b.path);
		vm_load_mm(&qa, QDT, tp->a.path, 1);
		vm_load_mm(&qb, QDT, tp->b.path, 1);

		da->f.mul(da, db, &dc, NAIVE);

		CASSERTION_TIME();
		qa->f.mul(qa, qb, &qc, NAIVE);
		CASSERTION(dc->f.compare(dc, qc) == 0, "a=%s,b=%s", tp->a.path,
				tp->b.path);

#if 0
		printf("---- a:\n");
		a->f.print(a);
		printf("---- b:\n");
		b->f.print(b);
		printf("---- naive:\n");
		c_def->f.print(c_def);
		printf("---- recursive:\n");
		c_str->f.print(c_str);
#endif

		da->f.free(da);
		db->f.free(db);
		dc->f.free(dc);
		qa->f.free(qa);
		qb->f.free(qb);
		qc->f.free(qc);
		da = NULL;
		db = NULL;
		dc = NULL;
		qa = NULL;
		qb = NULL;
		qc = NULL;
	}
}

int main(int argc, char *argv[]) {

	CASSERTION_INIT(argc, argv);

	run();

	CASSERTION_RESULTS();

	return EXIT_SUCCESS;
}
