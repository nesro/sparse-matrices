/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
 * https://github.com/nesro/sparse-matrices
 */
#include <stdio.h>
#include "test_utils.h"

test_matrix_t tm_small[] = {
//
		{ MTX_DIR "2x2_4nz_01.mtx", 2, 2, 4 }, /**/
		{ MTX_DIR "2x2_4nz_02.mtx", 2, 2, 4 }, /**/
		{ { 0 } } /**/
};

test_matrices_pair_t tm_pairs[] = {
//
		{ /**/
		{ MTX_DIR "2x2_4nz_01.mtx", 2, 2, 4 }, /**/
		{ MTX_DIR "2x2_4nz_02.mtx", 2, 2, 4 } /**/
		}, /**/

		{ /**/
		{ MTX_DIR "4x4_4nz_01.mtx", 4, 4, 4 }, /**/
		{ MTX_DIR "4x4_4nz_02.mtx", 4, 4, 4 } /**/
		}, /**/

		{ /**/
		{ MTX_DIR "4x4_4nz_03.mtx", 4, 4, 4 }, /**/
		{ MTX_DIR "4x4_4nz_04.mtx", 4, 4, 4 } /**/
		}, /**/

		{ /**/
		{ MTX_DIR "4x4_4nz_07.mtx", 4, 4, 4 }, /**/
		{ MTX_DIR "4x4_4nz_08.mtx", 4, 4, 4 } /**/
		}, /**/

		{ /**/
		{ MTX_DIR "4x4_4nz_01.mtx", 4, 4, 4 }, /**/
		{ MTX_DIR "4x4_4nz_02.mtx", 4, 4, 4 } /**/
		}, /**/

		{ /**/
		{ MTX_DIR "4x4_8nz_01.mtx", 4, 4, 8 }, /**/
		{ MTX_DIR "4x4_8nz_01.mtx", 4, 4, 8 } /**/
		}, /**/

		{ /**/
		{ MTX_DIR "4x4_16nz_01.mtx", 4, 4, 16 }, /**/
		{ MTX_DIR "4x4_16nz_02.mtx", 4, 4, 16 } /**/
		}, /**/

		{ { { 0 } } } /**/
};

const test_matrix_t *foreach_matrix(const test_matrix_t *tm_array) {
	static int i = 0;

	if (tm_array[i].path[0] != '\0') {
		i++;
		return &tm_array[i - 1];
	} else {
		i = 0;
		return NULL;
	}
}

const test_matrices_pair_t *foreach_pair(const test_matrices_pair_t *tm_p_array) {
	static int i = 0;

	if (tm_p_array[i].a.path[0] != '\0') {
		i++;
		return &tm_p_array[i - 1];
	} else {
		i = 0;
		return NULL;
	}
}
