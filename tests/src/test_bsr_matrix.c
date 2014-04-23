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
	vm_t *dc = NULL;
	vm_t *ba = NULL;
	vm_t *bb = NULL;
	vm_t *bc = NULL;
	int bs; /* block size */

#if BSR_DEBUG == 1
	vm_t *bsr_den = NULL;
#endif

	const test_matrices_pair_t *tp;

	while ((tp = foreach_pair(bsr_pairs)) != NULL) {
		printf("%s - %s\n", tp->a.path, tp->b.path);

		for (bs = 1; bs <= 128; bs *= 2) {

			vm_load_mm(&da, DEN, tp->a.path);
			vm_load_mm(&db, DEN, tp->b.path);

			vm_load_mm(&ba, BSR, tp->a.path, bs);

#if BSR_DEBUG == 1
			printf("=============================================================\n");
			bsr_den = ba->f.convert(ba, DEN);
			bsr_den->f.print(bsr_den);
			printf("=============================================================\n");
#endif

			vm_load_mm(&bb, BSR, tp->b.path, bs);

			da->f.mul(da, db, &dc, NAIVE);

			CASSERTION_TIME();
			ba->f.mul(ba, bb, &bc, NAIVE);

			CASSERTION(dc->f.compare(dc, bc) == 0, "a=%s,b=%s,bs=%d bsr",
					tp->a.path, tp->b.path, bs);

#if BSR_DEBUG == 1
			printf("---- dense:\n");
			dc->f.print(dc);
			printf("---- bsr:\n");
			bc->f.print(bc);
#endif

			da->f.free(da);
			db->f.free(db);
			dc->f.free(dc);
			ba->f.free(ba);
			bb->f.free(bb);
			bc->f.free(bc);

			da = NULL;
			db = NULL;
			dc = NULL;
			ba = NULL;
			bb = NULL;
			bc = NULL;

#if BSR_MATRIX == 1
			bsr_den->f.free(bsr_den);
			bsr_den = NULL;
#endif
		}
	}
}

int main(int argc, char *argv[]) {

	CASSERTION_INIT(argc, argv);

	run();

	CASSERTION_RESULTS();

	return EXIT_SUCCESS;
}
