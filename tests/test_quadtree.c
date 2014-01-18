#include <stdio.h>
#include <stdlib.h>

#include "test_framework.h"
#include "../src/qt_matrix.h"

/**
 * I had a bug in computing the quadtree height.
 * The bug has been fixed, but this test is so fast I will not delete it.
 */
static void quadtree_height_test(void) {

	qt_matrix_t qt_matrix;
	int size;
	int sm_size;

	size = 8;
	sm_size = 2;
	qt_matrix_init(&qt_matrix, size, size, 1, sm_size);
	snprintf(str_buf, STR_BUF_LEN, "size=%d, qt height=%d, sm_size=%d", size,
		qt_matrix.height, sm_size);
	TEST_ASSERT(qt_matrix.height == 2, str_buf);
	qt_matrix_free(&qt_matrix);

	size = 8;
	sm_size = 4;
	qt_matrix_init(&qt_matrix, size, size, 1, sm_size);
	snprintf(str_buf, STR_BUF_LEN, "size=%d, qt height=%d, sm_size=%d", size,
		qt_matrix.height, sm_size);
	TEST_ASSERT(qt_matrix.height == 1, str_buf);
	qt_matrix_free(&qt_matrix);

	size = 16;
	sm_size = 2;
	qt_matrix_init(&qt_matrix, size, size, 1, sm_size);
	snprintf(str_buf, STR_BUF_LEN, "size=%d, qt height=%d, sm_size=%d", size,
		qt_matrix.height, sm_size);
	TEST_ASSERT(qt_matrix.height == 3, str_buf);
	qt_matrix_free(&qt_matrix);

	size = 16;
	sm_size = 4;
	qt_matrix_init(&qt_matrix, size, size, 1, sm_size);
	snprintf(str_buf, STR_BUF_LEN, "size=%d, qt height=%d, sm_size=%d", size,
		qt_matrix.height, sm_size);
	TEST_ASSERT(qt_matrix.height == 2, str_buf);
	qt_matrix_free(&qt_matrix);

	size = 16;
	sm_size = 8;
	qt_matrix_init(&qt_matrix, size, size, 1, sm_size);
	snprintf(str_buf, STR_BUF_LEN, "size=%d, qt height=%d, sm_size=%d", size,
		qt_matrix.height, sm_size);
	TEST_ASSERT(qt_matrix.height == 1, str_buf);
	qt_matrix_free(&qt_matrix);
}

/**
 * This test loads a matrix from a file with the MatrixMarket format in both
 * quadtree and dense formats. The quadtree matrix is converted to dense format.
 * Original dense matrix and the dense matrix from the quadtree are compared.
 */
static void quadtree_load(const char *filename, int sm_size) {

	qt_matrix_t qt_matrix;
	dense_matrix_t den_from_qt_matrix;
	dense_matrix_t den_matrix;

	qt_matrix_load_mm(&qt_matrix, filename, sm_size);
	qt_matrix_to_dense(&qt_matrix, &den_from_qt_matrix);

	dense_matrix_load_mm(&den_matrix, filename);

	dense_to_html(&den_from_qt_matrix, "load_qt.html");
	dense_to_html(&den_matrix, "load_den.html");

	qt_matrix_print(&qt_matrix);

	snprintf(str_buf, STR_BUF_LEN, "file: %s", filename);
	TEST_ASSERT(dense_matrix_compare(&den_matrix, &den_from_qt_matrix) == 1,
		str_buf);

	qt_matrix_free(&qt_matrix);
	dense_matrix_free(&den_from_qt_matrix);
	dense_matrix_free(&den_matrix);
}

static void quadtree_load_test(void) {

	int sm_size;

	for (sm_size = 2; sm_size <= 32; sm_size *= 2) {
		quadtree_load("./tests/matrices/32x32_1024nz_01.mtx", sm_size);
		quadtree_load("./tests/matrices/32x32_1024nz_02.mtx", sm_size);

		if (sm_size > 16)
			continue;

		quadtree_load("./tests/matrices/16x16_256nz_01.mtx", sm_size);
		quadtree_load("./tests/matrices/16x16_256nz_02.mtx", sm_size);

		if (sm_size > 8)
			continue;

		quadtree_load("./tests/matrices/8x8_16nz_01.mtx", sm_size);
		quadtree_load("./tests/matrices/8x8_16nz_02.mtx", sm_size);

		quadtree_load("./tests/matrices/8x8_64nz.mtx", sm_size);

		if (sm_size > 4)
			continue;

		quadtree_load("./tests/matrices/4x4_4nz_01.mtx", sm_size);
		quadtree_load("./tests/matrices/4x4_4nz_02.mtx", sm_size);

		quadtree_load("./tests/matrices/4x4_10nz_01.mtx", sm_size);

		quadtree_load("./tests/matrices/4x4_16nz_01.mtx", sm_size);
		quadtree_load("./tests/matrices/4x4_16nz_02.mtx", sm_size);
	}
}

static void quadtree_mul(const char *filename_a, const char *filename_b,
	int sm_size, int html) {

	qt_matrix_t qtr_a;
	qt_matrix_t qtr_b;
	dense_matrix_t den_a;
	dense_matrix_t den_b;
	dense_matrix_t den_qtr_c;
	dense_matrix_t den_den_c;

	dense_matrix_t tmp_a;
	dense_matrix_t tmp_b;

	dense_matrix_load_mm(&den_a, filename_a);
	dense_matrix_load_mm(&den_b, filename_b);
	dense_matrix_matrix_mul(&den_a, &den_b, &den_den_c);
	dense_matrix_free(&den_a);
	dense_matrix_free(&den_b);

	qt_matrix_load_mm(&qtr_a, filename_a, sm_size);
	qt_matrix_load_mm(&qtr_b, filename_b, sm_size);

	qt_matrix_to_dense(&qtr_a, &tmp_a);
	qt_matrix_to_dense(&qtr_b, &tmp_b);

	if (html) {
		dense_to_html(&tmp_a, "qtr_a.html");
		dense_to_html(&tmp_b, "qtr_b.html");
	}

	qt_matrix_matrix_mul(&qtr_a, &qtr_b, &den_qtr_c);

	qt_matrix_print(&qtr_a);
	qt_matrix_print(&qtr_b);

	dense_matrix_free(&tmp_a);
	dense_matrix_free(&tmp_b);
	qt_matrix_free(&qtr_a);
	qt_matrix_free(&qtr_b);

	if (html) {
		dense_to_html(&den_qtr_c, "den_qtr_c.html");
		dense_to_html(&den_den_c, "den_den_c.html");
	}

	snprintf(str_buf, STR_BUF_LEN, "files: a=%s, b=%s sm_size=%d", filename_a,
		filename_b, sm_size);
	TEST_ASSERT(dense_matrix_compare(&den_qtr_c, &den_den_c) == 1, str_buf);

	dense_matrix_free(&den_qtr_c);
	dense_matrix_free(&den_den_c);
}

static void quadtree_mul_test(void) {

	int sm_size;

	/**
	 * TODO: multiply matrices with less than 100 % nnz
	 */

	for (sm_size = 2; sm_size <= 4; sm_size *= 2) {
		quadtree_mul("./tests/matrices/4x4_4nz_01.mtx",
			"./tests/matrices/4x4_4nz_02.mtx", sm_size, 1);

		quadtree_mul("./tests/matrices/4x4_4nz_01.mtx",
			"./tests/matrices/4x4_4nz_02.mtx", sm_size, 1);

		quadtree_mul("./tests/matrices/8x8_16nz_01.mtx",
			"./tests/matrices/8x8_16nz_02.mtx", sm_size, 1);

		quadtree_mul("./tests/matrices/16x16_256nz_01.mtx",
			"./tests/matrices/16x16_256nz_02.mtx", sm_size, 1);
	}

	for (sm_size = 2; sm_size <= 128; sm_size *= 2) {
		quadtree_mul("./tests/matrices/128_01.mtx",
			"./tests/matrices/128_02.mtx", sm_size, 0);
	}
}

static void csr(void) {
	csr_matrix_t csr_a;
	csr_matrix_t csr_b;
	dense_matrix_t den_c;

	csr_matrix_load_mm(&csr_a, "./tests/matrices/4096_01.mtx");
	csr_matrix_load_mm(&csr_b, "./tests/matrices/4096_02.mtx");

	printf("CSR MMM time: %lf <<<\n",
		csr_matrix_matrix_mul(&csr_a, &csr_b, &den_c));

	csr_matrix_free(&csr_a);
	csr_matrix_free(&csr_b);
	dense_matrix_free(&den_c);
}

static void quadtree_mul_notest(const char *filename_a, const char *filename_b,
	int sm_size) {
	qt_matrix_t qtr_a;
	qt_matrix_t qtr_b;

	dense_matrix_t den_qtr_c;

	qt_matrix_load_mm(&qtr_a, filename_a, sm_size);
	qt_matrix_load_mm(&qtr_b, filename_b, sm_size);

	printf("QUADTREE-CSR MMM submatrix_size=%d time: %lf <<<\n", sm_size,
		qt_matrix_matrix_mul(&qtr_a, &qtr_b, &den_qtr_c));

	qt_matrix_free(&qtr_a);
	qt_matrix_free(&qtr_b);

	dense_matrix_free(&den_qtr_c);

}

static void quadtree_vs_csr(void) {
//	quadtree_mul_notest("./tests/matrices/2048x2048_4194304nz_01.mtx",
//		"./tests/matrices/2048x2048_4194304nz_02.mtx", 128);

	quadtree_mul_notest("./tests/matrices/4096_01.mtx",
		"./tests/matrices/4096_02.mtx", 128);
	quadtree_mul_notest("./tests/matrices/4096_01.mtx",
		"./tests/matrices/4096_02.mtx", 256);
	quadtree_mul_notest("./tests/matrices/4096_01.mtx",
		"./tests/matrices/4096_02.mtx", 512);
	quadtree_mul_notest("./tests/matrices/4096_01.mtx",
		"./tests/matrices/4096_02.mtx", 1024);
}

int main(int argc, char *argv[]) {

	/* There is no random in this test. */
	/* srand(0); */

	TEST_RUN(quadtree_height_test);
	TEST_RUN(quadtree_load_test);
	TEST_RUN(quadtree_mul_test);
	TEST_RESULT(argv[0]);

//	csr();
//	quadtree_vs_csr();

	return EXIT_SUCCESS;
}
