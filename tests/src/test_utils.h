/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
 * https://github.com/nesro/sparse-matrices
 */

#ifndef TEST_UTILS_H_
#define TEST_UTILS_H_

#define MTX_DIR "./tests/matrices/"
#define MTX_GEN_DIR MTX_DIR "generated/"

typedef struct test_matrix {
	char path[50];
	int width;
	int height;
	int nnz;
} test_matrix_t;

typedef struct test_matrices_pair {
	test_matrix_t a;
	test_matrix_t b;
} test_matrices_pair_t;

/* tm = test matrices */
/* tp = test pair */
extern test_matrix_t tm_small[];
extern test_matrix_t kat_tm[];
extern test_matrices_pair_t tm_pairs[];

const test_matrix_t *foreach_matrix(const test_matrix_t *);
const test_matrices_pair_t *foreach_pair(const test_matrices_pair_t *);

#define FOREACH_MATRIX(tm, tm_a) \
	const test_matrix_t *tm; \
	while ((tm = foreach_matrix(tm_a)) != NULL)

#define FOREACH_PAIRS(tp, tm_p) \
	const test_matrices_pair_t *tp; \
	while ((tp = foreach_pair(tm_p)) != NULL)

#endif /* TEST_UTILS_H_ */
