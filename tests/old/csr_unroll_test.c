#include <stdio.h>
#include <stdlib.h>

#include "test_framework.h"

#include "../src/csr_matrix.h"
#include "../src/coo_matrix.h"

static void bacic2x2(void) {

	csr_t csr_a;
	csr_t csr_b;
	den_matrix_t dense_a;
	den_matrix_t dense_b;
	den_matrix_t dense_c;
	den_matrix_t dense_d;

	int html = 1;
	int i;

	char *matrices[][2] = { /**/
	/* 4nnz * 4nnz */
	{ "2x2_4nz_01", "2x2_4nz_02" }, /**/
	/**/
	/* 4nnz * 3nnz */
	{ "2x2_4nz_01", "2x2_3nz_01" }, /**/
	{ "2x2_4nz_01", "2x2_3nz_02" }, /**/
	{ "2x2_4nz_01", "2x2_3nz_03" }, /**/
	{ "2x2_4nz_01", "2x2_3nz_04" }, /**/
	/**/
	/* 3nnz * 4nnz */
	{ "2x2_3nz_01", "2x2_4nz_01" }, /**/
	{ "2x2_3nz_02", "2x2_4nz_01" }, /**/
	{ "2x2_3nz_03", "2x2_4nz_01" }, /**/
	{ "2x2_3nz_04", "2x2_4nz_01" }, /**/
	/**/
	/* 3nnz * 3nnz */
	{ "2x2_3nz_01", "2x2_3nz_01" }, /**/
	{ "2x2_3nz_01", "2x2_3nz_02" }, /**/
	{ "2x2_3nz_01", "2x2_3nz_03" }, /**/
	{ "2x2_3nz_01", "2x2_3nz_04" }, /**/
	{ "2x2_3nz_02", "2x2_3nz_01" }, /**/
	{ "2x2_3nz_02", "2x2_3nz_02" }, /**/
	{ "2x2_3nz_02", "2x2_3nz_03" }, /**/
	{ "2x2_3nz_02", "2x2_3nz_04" }, /**/
	{ "2x2_3nz_03", "2x2_3nz_01" }, /**/
	{ "2x2_3nz_03", "2x2_3nz_02" }, /**/
	{ "2x2_3nz_03", "2x2_3nz_03" }, /**/
	{ "2x2_3nz_03", "2x2_3nz_04" }, /**/
	{ "2x2_3nz_04", "2x2_3nz_01" }, /**/
	{ "2x2_3nz_04", "2x2_3nz_02" }, /**/
	{ "2x2_3nz_04", "2x2_3nz_03" }, /**/
	{ "2x2_3nz_04", "2x2_3nz_04" }, /**/
	/**/
	/* 4nnz * 2nnz */
	{ "2x2_4nz_01", "2x2_2nz_01" }, /**/
	{ "2x2_4nz_01", "2x2_2nz_02" }, /**/
	{ "2x2_4nz_01", "2x2_2nz_03" }, /**/
	{ "2x2_4nz_01", "2x2_2nz_04" }, /**/
	{ "2x2_4nz_01", "2x2_2nz_05" }, /**/
	{ "2x2_4nz_01", "2x2_2nz_06" }, /**/
	/**/
	/* 3nnz * 2nnz */
	{ "2x2_3nz_01", "2x2_2nz_01" }, /**/
	{ "2x2_3nz_01", "2x2_2nz_02" }, /**/
	{ "2x2_3nz_01", "2x2_2nz_03" }, /**/
	{ "2x2_3nz_01", "2x2_2nz_04" }, /**/
	{ "2x2_3nz_01", "2x2_2nz_05" }, /**/
	{ "2x2_3nz_01", "2x2_2nz_06" }, /**/
	{ "2x2_3nz_02", "2x2_2nz_01" }, /**/
	{ "2x2_3nz_02", "2x2_2nz_02" }, /**/
	{ "2x2_3nz_02", "2x2_2nz_03" }, /**/
	{ "2x2_3nz_02", "2x2_2nz_04" }, /**/
	{ "2x2_3nz_02", "2x2_2nz_05" }, /**/
	{ "2x2_3nz_02", "2x2_2nz_06" }, /**/
	{ "2x2_3nz_03", "2x2_2nz_01" }, /**/
	{ "2x2_3nz_03", "2x2_2nz_02" }, /**/
	{ "2x2_3nz_03", "2x2_2nz_03" }, /**/
	{ "2x2_3nz_03", "2x2_2nz_04" }, /**/
	{ "2x2_3nz_03", "2x2_2nz_05" }, /**/
	{ "2x2_3nz_03", "2x2_2nz_06" }, /**/
	{ "2x2_3nz_04", "2x2_2nz_01" }, /**/
	{ "2x2_3nz_04", "2x2_2nz_02" }, /**/
	{ "2x2_3nz_04", "2x2_2nz_03" }, /**/
	{ "2x2_3nz_04", "2x2_2nz_04" }, /**/
	{ "2x2_3nz_04", "2x2_2nz_05" }, /**/
	{ "2x2_3nz_04", "2x2_2nz_06" }, /**/
	/**/
	/* 2nnz * 2nnz */
	{ "2x2_2nz_01", "2x2_2nz_01" }, /**/
	{ "2x2_2nz_02", "2x2_2nz_01" }, /**/
	{ "2x2_2nz_03", "2x2_2nz_01" }, /**/
	{ "2x2_2nz_04", "2x2_2nz_01" }, /**/
	{ "2x2_2nz_05", "2x2_2nz_01" }, /**/
	{ "2x2_2nz_06", "2x2_2nz_01" }, /**/
	{ "2x2_2nz_01", "2x2_2nz_02" }, /**/
	{ "2x2_2nz_02", "2x2_2nz_02" }, /**/
	{ "2x2_2nz_03", "2x2_2nz_02" }, /**/
	{ "2x2_2nz_04", "2x2_2nz_02" }, /**/
	{ "2x2_2nz_05", "2x2_2nz_02" }, /**/
	{ "2x2_2nz_06", "2x2_2nz_02" }, /**/
	{ "2x2_2nz_01", "2x2_2nz_03" }, /**/
	{ "2x2_2nz_02", "2x2_2nz_03" }, /**/
	{ "2x2_2nz_03", "2x2_2nz_03" }, /**/
	{ "2x2_2nz_04", "2x2_2nz_03" }, /**/
	{ "2x2_2nz_05", "2x2_2nz_03" }, /**/
	{ "2x2_2nz_06", "2x2_2nz_03" }, /**/
	{ "2x2_2nz_01", "2x2_2nz_04" }, /**/
	{ "2x2_2nz_02", "2x2_2nz_04" }, /**/
	{ "2x2_2nz_03", "2x2_2nz_04" }, /**/
	{ "2x2_2nz_04", "2x2_2nz_04" }, /**/
	{ "2x2_2nz_05", "2x2_2nz_04" }, /**/
	{ "2x2_2nz_06", "2x2_2nz_04" }, /**/
	{ "2x2_2nz_01", "2x2_2nz_05" }, /**/
	{ "2x2_2nz_02", "2x2_2nz_05" }, /**/
	{ "2x2_2nz_03", "2x2_2nz_05" }, /**/
	{ "2x2_2nz_04", "2x2_2nz_05" }, /**/
	{ "2x2_2nz_05", "2x2_2nz_05" }, /**/
	{ "2x2_2nz_06", "2x2_2nz_05" }, /**/
	{ "2x2_2nz_01", "2x2_2nz_06" }, /**/
	{ "2x2_2nz_02", "2x2_2nz_06" }, /**/
	{ "2x2_2nz_03", "2x2_2nz_06" }, /**/
	{ "2x2_2nz_04", "2x2_2nz_06" }, /**/
	{ "2x2_2nz_05", "2x2_2nz_06" }, /**/
	{ "2x2_2nz_06", "2x2_2nz_06" }, /**/
	/**/
	/* 2nnz * 4nnz */
	{ "2x2_2nz_01", "2x2_4nz_01" }, /**/
	{ "2x2_2nz_02", "2x2_4nz_01" }, /**/
	{ "2x2_2nz_03", "2x2_4nz_01" }, /**/
	{ "2x2_2nz_04", "2x2_4nz_01" }, /**/
	{ "2x2_2nz_05", "2x2_4nz_01" }, /**/
	{ "2x2_2nz_06", "2x2_4nz_01" }, /**/
	{ "2x2_2nz_01", "2x2_4nz_02" }, /**/
	{ "2x2_2nz_02", "2x2_4nz_02" }, /**/
	{ "2x2_2nz_03", "2x2_4nz_02" }, /**/
	{ "2x2_2nz_04", "2x2_4nz_02" }, /**/
	{ "2x2_2nz_05", "2x2_4nz_02" }, /**/
	{ "2x2_2nz_06", "2x2_4nz_02" }, /**/
	/**/
	/* 2nnz * 3nnz */
	{ "2x2_2nz_01", "2x2_3nz_01" }, /**/
	{ "2x2_2nz_02", "2x2_3nz_01" }, /**/
	{ "2x2_2nz_03", "2x2_3nz_01" }, /**/
	{ "2x2_2nz_04", "2x2_3nz_01" }, /**/
	{ "2x2_2nz_05", "2x2_3nz_01" }, /**/
	{ "2x2_2nz_06", "2x2_3nz_01" }, /**/
	{ "2x2_2nz_01", "2x2_3nz_02" }, /**/
	{ "2x2_2nz_02", "2x2_3nz_02" }, /**/
	{ "2x2_2nz_03", "2x2_3nz_02" }, /**/
	{ "2x2_2nz_04", "2x2_3nz_02" }, /**/
	{ "2x2_2nz_05", "2x2_3nz_02" }, /**/
	{ "2x2_2nz_06", "2x2_3nz_02" }, /**/
	{ "2x2_2nz_01", "2x2_3nz_02" }, /**/
	{ "2x2_2nz_02", "2x2_3nz_03" }, /**/
	{ "2x2_2nz_03", "2x2_3nz_03" }, /**/
	{ "2x2_2nz_04", "2x2_3nz_03" }, /**/
	{ "2x2_2nz_05", "2x2_3nz_03" }, /**/
	{ "2x2_2nz_06", "2x2_3nz_03" }, /**/
	{ "2x2_2nz_01", "2x2_3nz_02" }, /**/
	{ "2x2_2nz_02", "2x2_3nz_04" }, /**/
	{ "2x2_2nz_03", "2x2_3nz_04" }, /**/
	{ "2x2_2nz_04", "2x2_3nz_04" }, /**/
	{ "2x2_2nz_05", "2x2_3nz_04" }, /**/
	{ "2x2_2nz_06", "2x2_3nz_04" }, /**/

	{ NULL, NULL } /**/
	};

	for (i = 0; i < 9999; i++) {

		if (matrices[i][0] == NULL) {
			break;
		}

		snprintf(str_buf, STR_BUF_LEN, "%s/%s.mtx", TEST_MATRICES_DIR,
			matrices[i][0]);
		csr_matrix_load_mm(&csr_a, str_buf);

		snprintf(str_buf, STR_BUF_LEN, "%s/%s.mtx", TEST_MATRICES_DIR,
			matrices[i][1]);
		csr_matrix_load_mm(&csr_b, str_buf);

		if (html) {
			csr_to_html(&csr_a, "orig_a.html");
			csr_to_html(&csr_b, "orig_b.html");
		}

		csr_matrix_matrix_mul(&csr_a, &csr_b, &dense_a);
		csr_matrix_matrix_mul_unrolled(&csr_a, &csr_b, &dense_b);
		csr_matrix_matrix_mul_parallel(&csr_a, &csr_b, &dense_c);
		csr_matrix_matrix_mul_unrolled_parallel(&csr_a, &csr_b, &dense_d);

		csr_matrix_free(&csr_a);
		csr_matrix_free(&csr_b);

		if (html) {
			dense_to_html(&dense_a, "csrmmm.html");
			dense_to_html(&dense_b, "csrmmm_unroll.html");
			dense_to_html(&dense_c, "csrmmm_parallel.html");
			dense_to_html(&dense_d, "csrmmm_unrolled_parallel.html");
		}

		snprintf(str_buf, STR_BUF_LEN, "%s.mtx vs %s.mtx normal-unrolled",
			matrices[i][0], matrices[i][1]);
		TEST_ASSERT(dense_matrix_compare(&dense_a, &dense_b) == 1, str_buf);

		snprintf(str_buf, STR_BUF_LEN, "%s.mtx vs %s.mtx normal-parallel",
			matrices[i][0], matrices[i][1]);
		TEST_ASSERT(dense_matrix_compare(&dense_a, &dense_c) == 1, str_buf);

		snprintf(str_buf, STR_BUF_LEN,
			"%s.mtx vs %s.mtx normal-unrolled-parallel", matrices[i][0],
			matrices[i][1]);
		TEST_ASSERT(dense_matrix_compare(&dense_a, &dense_d) == 1, str_buf);

		dense_matrix_free(&dense_a);
		dense_matrix_free(&dense_b);
		dense_matrix_free(&dense_c);
		dense_matrix_free(&dense_d);
	}
}

void csr_mmm_test(void) {
	csr_t csr_a = { { 0 } };
	csr_t csr_b = { { 0 } };

	den_matrix_t dense_c = { { 0 } };
	den_matrix_t dense_c_unrolled = { { 0 } };
	den_matrix_t dense_c_parallel = { { 0 } };
	den_matrix_t dense_c_unrolled_parallel = { { 0 } };

	int i, j;

	for (i = 1; i <= 10; i++) {
		for (j = 1; j <= i * i; j++) {

			csr_matrix_generate(&csr_a, i, i, j);
			csr_matrix_generate(&csr_b, i, i, j);

			csr_matrix_matrix_mul(&csr_a, &csr_b, &dense_c);

			csr_matrix_matrix_mul_unrolled(&csr_a, &csr_b, &dense_c_unrolled);

			csr_matrix_matrix_mul_parallel(&csr_a, &csr_b, &dense_c_parallel);

			csr_matrix_matrix_mul_unrolled_parallel(&csr_a, &csr_b,
				&dense_c_unrolled_parallel);

			snprintf(str_buf, STR_BUF_LEN, "i=%d j=%d ", i, j);
			TEST_ASSERT(dense_matrix_compare(&dense_c, &dense_c_unrolled) == 1,
				str_buf);

			snprintf(str_buf, STR_BUF_LEN, "i=%d j=%d ", i, j);
			TEST_ASSERT(dense_matrix_compare(&dense_c, &dense_c_parallel) == 1,
				str_buf);

			snprintf(str_buf, STR_BUF_LEN, "i=%d j=%d ", i, j);
			TEST_ASSERT(
				dense_matrix_compare(&dense_c, &dense_c_unrolled_parallel) == 1,
				str_buf);

			csr_matrix_free(&csr_a);
			csr_matrix_free(&csr_b);
			dense_matrix_free(&dense_c);
			dense_matrix_free(&dense_c_unrolled);
			dense_matrix_free(&dense_c_parallel);
			dense_matrix_free(&dense_c_unrolled_parallel);
		}
	}
}

void csr_mvm_test(void) {
	csr_t csr_a = { { 0 } };
	vec_t vector_b = { 0 };

	vec_t vector_c = { 0 };
	vec_t vector_c_unrolled = { 0 };
	vec_t vector_c_parallel = { 0 };
	vec_t vector_c_unrolled_parallel = { 0 };

	int i, j;

	for (i = 1; i <= 10; i++) {
		for (j = 1; j <= i * i; j++) {

			srand(0);

			csr_matrix_generate(&csr_a, i, i, j);
			vector_generate(&vector_b, i);

			csr_matrix_vector_mul(&csr_a, &vector_b, &vector_c);

			csr_matrix_vector_mul_unrolled(&csr_a, &vector_b,
				&vector_c_unrolled);

			csr_matrix_vector_mul_parallel(&csr_a, &vector_b,
				&vector_c_parallel);

			csr_matrix_vector_mul_unrolled_parallel(&csr_a, &vector_b,
				&vector_c_unrolled_parallel);

			snprintf(str_buf, STR_BUF_LEN, "i=%d j=%d ", i, j);
			TEST_ASSERT(vector_compare(&vector_c, &vector_c_unrolled) == 1,
				str_buf);

			snprintf(str_buf, STR_BUF_LEN, "i=%d j=%d ", i, j);
			TEST_ASSERT(vector_compare(&vector_c, &vector_c_parallel) == 1,
				str_buf);

			snprintf(str_buf, STR_BUF_LEN, "i=%d j=%d ", i, j);
			TEST_ASSERT(
				vector_compare(&vector_c, &vector_c_unrolled_parallel) == 1,
				str_buf);

			csr_matrix_free(&csr_a);
			vector_free(&vector_b);
			vector_free(&vector_c);
			vector_free(&vector_c_unrolled);
			vector_free(&vector_c_parallel);
			vector_free(&vector_c_unrolled_parallel);
		}
	}
}

int main(int argc, char *argv[]) {

	/* We have random in this test. We want to have a deterministic test. */
	srand(0);

	printf("\n");

	TEST_RUN(bacic2x2);
	TEST_RUN(csr_mmm_test);
	TEST_RUN(csr_mvm_test);

	TEST_RESULT(argv[0]);

	return EXIT_SUCCESS;
}
