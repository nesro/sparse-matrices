/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
 * https://github.com/nesro/sparse-matrices
 */


/*

debug cheatsheet:

make DEBUG=1 tests && gdb ./tests/bin/test_bsr_matrix
break test_bsr_matrix.c:33
run
print ((bsr_t *)ba)->v[0]

*/

#include <stdio.h>
#include <stdlib.h>

#include "../../src/virtual_matrix.h"
#include "../../src/bsr_matrix.h"

#include "../../cassertion/cassertion.h"
#include "test_utils.h"

static void run() {

	vm_t *da = NULL;
	vm_t *db = NULL;
//	vm_t *dc = NULL;
	vm_t *ba = NULL;
	vm_t *bb = NULL;
//	vm_t *cc = NULL;

	const test_matrices_pair_t *tp;

	while ((tp = foreach_pair(bsr_pairs)) != NULL) {
		printf("%s - %s\n", tp->a.path, tp->b.path);

		vm_load_mm(&da, DEN, tp->a.path);
		vm_load_mm(&db, DEN, tp->b.path);

		vm_load_mm(&ba, BSR, tp->a.path, 2);

		printf("=============================================================\n");

		vm_load_mm(&bb, BSR, tp->b.path, 2);

//		da->f.mul(da, db, &dc, NAIVE);
//
//		CASSERTION_TIME();
//		ba->f.mul(ba, bb, &cc, NAIVE);
//
//		CASSERTION(dc->f.compare(dc, cc) == 0, "a=%s,b=%s, bsr", tp->a.path,
//				tp->b.path);
//
//#if 0
//		printf("---- dense:\n");
//		dc->f.print(dc);
//		printf("---- bsr:\n");
//		cc->f.print(cc);
//#endif

		da->f.free(da);
		db->f.free(db);
//		dc->f.free(dc);
		ba->f.free(ba);
		bb->f.free(bb);
//		cc->f.free(cc);

		da = NULL;
		db = NULL;
//		dc = NULL;
		ba = NULL;
		bb = NULL;
//		cc = NULL;

		break;
	}
}

int main(int argc, char *argv[]) {

	CASSERTION_INIT(argc, argv);

	run();

	CASSERTION_RESULTS();

	return EXIT_SUCCESS;
}
