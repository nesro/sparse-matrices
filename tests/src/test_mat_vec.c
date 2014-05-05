/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
 * https://github.com/nesro/sparse-matrices
 */

#include <stdio.h>
#include <stdlib.h>

#include "../../src/virtual_matrix.h"
#include "../../src/den_matrix.h"
#include "../../src/coo_matrix.h"
#include "../../src/csr_matrix.h"
#include "../../src/kat_matrix.h"
#include "../../src/vector.h"

#include "../../cassertion/cassertion.h"
#include "test_utils.h"

void mat_vec() {

	vm_t *a_den = NULL;
	vm_t *a_bsr = NULL;
	vm_t *a_coo = NULL;
	vm_t *a_csr = NULL;
	vm_t *a_kat = NULL;
	vm_t *b_vec = NULL;
	vm_t *c_vec_den = NULL;
	vm_t *c_vec_bsr = NULL;
	vm_t *c_vec_coo = NULL;
	vm_t *c_vec_csr = NULL;
	vm_t *c_vec_kat = NULL;
	const test_matrices_pair_t *tp;

	while ((tp = foreach_pair(mat_vec_pairs)) != NULL) {
		printf("%s - %s\n", tp->a.path, tp->b.path);

		vm_load_mm(&a_den, DEN, tp->a.path);
		vm_load_mm(&a_bsr, BSR, tp->a.path, 8);
		vm_load_mm(&a_coo, COO, tp->a.path);
		vm_load_mm(&a_csr, CSR, tp->a.path);
		vm_load_mm(&a_kat, KAT, tp->a.path, 8);

		vm_load_mm(&b_vec, VEC, tp->b.path, VA_END);
		a_den->f.mul(a_den, b_vec, &c_vec_den, NAIVE);

		CASSERTION_TIME();
		a_bsr->f.mul(a_bsr, b_vec, &c_vec_bsr, NAIVE);
		CASSERTION(c_vec_den->f.compare(c_vec_den, c_vec_bsr) == 0,
				"a=%s,b=%s, coo", tp->a.path, tp->b.path);

		CASSERTION_TIME();
		a_coo->f.mul(a_coo, b_vec, &c_vec_coo, NAIVE);
		CASSERTION(c_vec_den->f.compare(c_vec_den, c_vec_coo) == 0,
				"a=%s,b=%s, coo", tp->a.path, tp->b.path);

		CASSERTION_TIME();
		a_csr->f.mul(a_csr, b_vec, &c_vec_csr, NAIVE);
		CASSERTION(c_vec_den->f.compare(c_vec_den, c_vec_csr) == 0,
				"a=%s,b=%s, csr", tp->a.path, tp->b.path);

		CASSERTION_TIME();
		a_kat->f.mul(a_kat, b_vec, &c_vec_kat, NAIVE);
		CASSERTION(c_vec_den->f.compare(c_vec_den, c_vec_kat) == 0,
				"a=%s,b=%s, kat sms=%d", tp->a.path, tp->b.path, 2);

#if 0
		printf("---- dense:\n");
		dc->f.print(dc);
		printf("---- coo:\n");
		cc->f.print(cc);
#endif

		a_den->f.free(a_den);
		a_bsr->f.free(a_bsr);
		a_coo->f.free(a_coo);
		a_csr->f.free(a_csr);
		a_kat->f.free(a_kat);
		b_vec->f.free(b_vec);
		c_vec_den->f.free(c_vec_den);
		c_vec_bsr->f.free(c_vec_bsr);
		c_vec_coo->f.free(c_vec_coo);
		c_vec_csr->f.free(c_vec_csr);
		c_vec_kat->f.free(c_vec_kat);

		a_den = NULL;
		a_bsr = NULL;
		a_coo = NULL;
		a_csr = NULL;
		a_kat = NULL;
		b_vec = NULL;
		c_vec_den = NULL;
		c_vec_coo = NULL;
	}
}

/******************************************************************************/

vm_type_t types[] = { CSR, BSR, KAT, /**/ COO };
int types_size = 3;

void mat_mat() {

	int i;
	vm_t *den_a = NULL;
	vm_t *den_b = NULL;
	vm_t *den_c = NULL;
	vm_t *spa_a = NULL;
	vm_t *spa_b = NULL;
	vm_t *spa_c = NULL;
	const test_matrices_pair_t *tp;
	double time;
	int sms_start = 1;
//	int sms_stop = 32;
	int sms = sms_start;

//	while ((tp = foreach_pair(mat_mat_pairs)) != NULL) {
	while ((tp = foreach_pair(sparse_collection_pairs)) != NULL) {

		CASSERTION_MSG("begin dense a=%s,b=%s\n", tp->a.path, tp->b.path);

		vm_load_mm(&den_a, DEN, tp->a.path);
		vm_load_mm(&den_b, DEN, tp->b.path);
		den_a->f.mul(den_a, den_b, &den_c, NAIVE);

		for (i = 0; i < types_size; i++) {

			if (vm_has_blocks(types[i])) {
				block_loop: /**/
				CASSERTION_MSG("begin format %d witch sms=%d\n", types[i], sms);
				vm_load_mm(&spa_a, types[i], tp->a.path, sms);

				if (types[i] == KAT) {
					printf("kat_n %d\n", KAT_N);
					printf("kat_sm_size %d\n",
							((kat_matrix_t*) spa_a)->sm_size);
					printf("kat_a_inner %d\n",
							((kat_matrix_t*) spa_a)->nodes_inner);
					printf("kat_a_dense %d\n",
							((kat_matrix_t*) spa_a)->nodes_den);
					printf("kat_a_csr %d\n",
							((kat_matrix_t*) spa_a)->nodes_csr);
				}

				vm_load_mm(&spa_b, types[i], tp->b.path, sms);
			} else {
				CASSERTION_MSG("begin format %d\n", types[i]);
				vm_load_mm(&spa_a, types[i], tp->a.path);
				vm_load_mm(&spa_b, types[i], tp->b.path);
			}

			CASSERTION_TIME();

			time = spa_a->f.mul(spa_a, spa_b, &spa_c, NAIVE);

			CASSERTION(den_c->f.compare(den_c, spa_c) == 0,
					"a=%s,b=%s,sms=%d,time=%lf,format=%d", tp->a.path,
					tp->b.path, sms, time, types[i]);

#if 0
//			printf("---- dense bef:\n");
//			den_c->f.print(den_a);
			printf("---- dense:\n");
			den_c->f.print(den_c);
			printf("---- cspa:\n");
			spa_c->f.print(spa_c);
#endif

			spa_a->f.free(spa_a);
			spa_b->f.free(spa_b);
			spa_c->f.free(spa_c);
			spa_a = NULL;
			spa_b = NULL;
			spa_c = NULL;

			if (vm_has_blocks(types[i])) {
				sms *= 2;
#if 0
				if (sms <= sms_stop)
#else
				if (sms <= tp->a.height)
#endif
					goto block_loop;
				else
					sms = sms_start;
			}
		}

		den_a->f.free(den_a);
		den_b->f.free(den_b);
		den_c->f.free(den_c);
		den_a = NULL;
		den_b = NULL;
		den_c = NULL;
	}
}

int main(int argc, char *argv[]) {

	CASSERTION_INIT(argc, argv);

//	mat_vec();
	mat_mat();

	CASSERTION_RESULTS();

	return EXIT_SUCCESS;
}
