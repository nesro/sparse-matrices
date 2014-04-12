/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
 * https://github.com/nesro/sparse-matrices
 */

#include <stdio.h>
#include <stdlib.h>

#include "../../src/virtual_matrix.h"
#include "../../src/kat_matrix.h"

#include "../../cassertion/cassertion.h"
#include "test_utils.h"

static void load() {

	vm_t *kat = NULL;
	const test_matrix_t *tm;

	while ((tm = foreach_matrix(kat_tm)) != NULL) {
		printf("%s\n", tm->path);
		vm_load_mm(&kat, KAT, tm->path, 4);
		kat->f.free(kat);
		kat = NULL;

		CASSERTION(1, "I survived: %s", tm->path);
	}
}

/*static void run() {

 vm_t *da = NULL;
 vm_t *db = NULL;
 vm_t *dc = NULL;
 vm_t *ka = NULL;
 vm_t *kb = NULL;
 vm_t *kc = NULL;

 const test_matrix_t *tm;
 const test_matrices_pair_t *tp;

 //	while ((tm = foreach_matrix(kat_tm_small)) != NULL) {
 //		printf("%s\n", tm->path);
 //	}

 while ((tp = foreach_pair(kat_tm_pairs)) != NULL) {
 printf("%s - %s\n", tp->a.path, tp->b.path);

 vm_load_mm(&da, DEN, tp->a.path);
 vm_load_mm(&db, DEN, tp->b.path);
 vm_load_mm(&da, KAT, tp->a.path, 1);
 vm_load_mm(&db, KAT, tp->b.path, 1);

 da->f.mul(da, db, &dc, NAIVE);

 CASSERTION_TIME();
 ka->f.mul(ka, kb, &kc, NAIVE);
 CASSERTION(dc->f.compare(dc, kc) == 0, "a=%s,b=%s", tp->a.path,
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
 ka->f.free(ka);
 kb->f.free(kb);
 kc->f.free(kc);
 da = NULL;
 db = NULL;
 dc = NULL;
 ka = NULL;
 kb = NULL;
 kc = NULL;
 }
 }*/

int main(int argc, char *argv[]) {

	CASSERTION_INIT(argc, argv);

	load();

	//run();

	CASSERTION_RESULTS();

	return EXIT_SUCCESS;
}

