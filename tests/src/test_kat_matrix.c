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

void load() {

	vm_t *kat = NULL;
	const test_matrix_t *tm;

	while ((tm = foreach_matrix(kat_tm)) != NULL) {
		printf("%s\n", tm->path);
		vm_load_mm(&kat, KAT, tm->path, 128);


#if KAT_DEBUG == 1
		vm_t *den = NULL;
		den = kat->f.convert(kat, DEN);
		printf(
				"=============================================================\n");
		den->f.print(den);
		printf(
				"=============================================================\n");
		den->f.free(den);
		den = NULL;
#endif

		kat->f.free(kat);
		kat = NULL;

		CASSERTION(1, "I survived: %s", tm->path);
	}
}

void run() {

	vm_t *da = NULL;
	vm_t *db = NULL;
	vm_t *dc = NULL;
	vm_t *ka = NULL;
	vm_t *kb = NULL;
	vm_t *kc = NULL;

	const test_matrices_pair_t *tp;
	double time;

	int sms; /* submatrix size */

	//	while ((tm = foreach_matrix(kat_tm_small)) != NULL) {
	//		printf("%s\n", tm->path);
	//	}

	while ((tp = foreach_pair(kat_tm_pairs)) != NULL) {

		vm_load_mm(&da, DEN, tp->a.path);
		vm_load_mm(&db, DEN, tp->b.path);
		da->f.mul(da, db, &dc, NAIVE);

		for (sms = 64; sms <= 1024; sms *= 2) {

			CASSERTION_DONTRUN(sms > tp->a.height,
					"sms=%d > tp->a.height=%d\n", sms, tp->a.height);
			if (cassertion.dontrun)
				continue;

			CASSERTION_MSG("begin a=%s,b=%s,sms=%d,kat_n=%d\n", tp->a.path,
					tp->b.path, sms, KAT_N);

			vm_load_mm(&ka, KAT, tp->a.path, sms);
			vm_load_mm(&kb, KAT, tp->b.path, sms);

			CASSERTION_TIME();
			time = ka->f.mul(ka, kb, &kc, NAIVE);
			CASSERTION(dc->f.compare(dc, kc) == 0, "a=%s,b=%s,sms=%d,time=%lf",
					tp->a.path, tp->b.path, sms, time);

#ifdef KAT_DEBUG
#if KAT_DEBUG
			printf("---- den:\n");
			dc->f.print(dc);
			printf("----  kat:\n");
			kc->f.print(kc);
#endif /* if KAT_DEBUG */
#endif /* ifdef KAT_DEBUG */

			ka->f.free(ka);
			kb->f.free(kb);
			kc->f.free(kc);
			ka = NULL;
			kb = NULL;
			kc = NULL;
		}

		da->f.free(da);
		db->f.free(db);
		dc->f.free(dc);
		da = NULL;
		db = NULL;
		dc = NULL;
	}
}

int main(int argc, char *argv[]) {

	CASSERTION_INIT(argc, argv);

//	load();

	run();

	CASSERTION_RESULTS();

	return EXIT_SUCCESS;
}

