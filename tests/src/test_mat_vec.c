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


vm_type_t types[] = { KAT, CSR, BSR };
int types_size = 3;

/******************************************************************************/

void mat_vec() {

	int i;
	vm_t *den_a = NULL;
	vm_t *vec_b = NULL;
	vm_t *den_c = NULL;
	vm_t *spa_a = NULL;
	vm_t *spa_c = NULL;
	const test_matrices_pair_t *tp;
	double time;
	int sms_start = 1;
//	int sms_stop = 32;
	int sms = sms_start;

//	while ((tp = foreach_pair(mat_mat_pairs)) != NULL) {
	while ((tp = foreach_pair(sparse_collection_vectors_pairs)) != NULL) {

		CASSERTION_MSG("begin dense a=%s,b=%s\n", tp->a.path, tp->b.path);

		vm_load_mm(&den_a, DEN, tp->a.path);
		vm_load_mm(&vec_b, VEC, tp->b.path, den_a->w);
		den_a->f.mul(den_a, vec_b, &den_c, NAIVE);

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

			} else {
				CASSERTION_MSG("begin format %d\n", types[i]);
				vm_load_mm(&spa_a, types[i], tp->a.path);
			}

			CASSERTION_TIME();

			time = spa_a->f.mul(spa_a, vec_b, &spa_c, NAIVE);

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
			spa_c->f.free(spa_c);
			spa_a = NULL;
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
		vec_b->f.free(vec_b);
		den_c->f.free(den_c);
		den_a = NULL;
		vec_b = NULL;
		den_c = NULL;
	}
}

/******************************************************************************/

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

	mat_vec();
	mat_mat();

	CASSERTION_RESULTS();

	return EXIT_SUCCESS;
}
