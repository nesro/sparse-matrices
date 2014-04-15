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

		{ /**/
		{ MTX_GEN_DIR "16x16_01.mtx", 16, 16, 1024 }, /**/
		{ MTX_GEN_DIR "16x16_02.mtx", 16, 16, 1024 } /**/
		}, /**/

		{ /**/
		{ MTX_GEN_DIR "32x32_01.mtx", 32, 32, 1024 }, /**/
		{ MTX_GEN_DIR "32x32_02.mtx", 32, 32, 1024 } /**/
		}, /**/

		{ /**/
		{ MTX_GEN_DIR "64x64_01.mtx", 64, 64, 4096 }, /**/
		{ MTX_GEN_DIR "64x64_02.mtx", 64, 64, 4096 }, /**/
		}, /**/

		{ /**/
		{ MTX_GEN_DIR "128x128_01.mtx", 128, 128, 16384 }, /**/
		{ MTX_GEN_DIR "128x128_02.mtx", 128, 128, 16384 }, /**/
		}, /**/

		{ { { 0 } } } /**/
};

test_matrix_t kat_tm[] = {
//
		{ MTX_DIR "64x64_2nz_01.mtx", 64, 64, 2 }, /**/
		{ MTX_DIR "64x64_2nz_02.mtx", 64, 64, 2 }, /**/
		{ MTX_DIR "64x64_4nz_01.mtx", 64, 64, 4 }, /**/
		{ MTX_GEN_DIR "64x64_01.mtx", 64, 64, 4096 }, /**/
		{ { 0 } }, /**/
};

test_matrices_pair_t kat_tm_pairs[] = {
//

		{ /**/
				{"/tmp/matrix.mtx", 4096, 4096, 167 }, /**/
				{"/tmp/matrix2.mtx", 4096, 4096, 167 }, /**/
				}, /**/

				{ { { 0 } } }, /**/

		{ /**/
		{ MTX_DIR "64x64_1nz_01.mtx", 64, 64, 2 }, /**/
		{ MTX_DIR "64x64_1nz_01.mtx", 64, 64, 2 }, /**/
		}, /**/


		{ /**/
		{ MTX_DIR "64x64_1nz_01.mtx", 64, 64, 2 }, /**/
		{ MTX_DIR "64x64_1nz_02.mtx", 64, 64, 2 }, /**/
		}, /**/

		{ /**/
		{ MTX_DIR "64x64_2nz_01.mtx", 64, 64, 2 }, /**/
		{ MTX_DIR "64x64_2nz_02.mtx", 64, 64, 2 }, /**/
		}, /**/


		{ /**/
		{ MTX_DIR "64x64_2nz_01.mtx", 64, 64, 2 }, /**/
		{ MTX_DIR "64x64_2nz_03.mtx", 64, 64, 2 }, /**/
		}, /**/


		{ /**/
		{ MTX_GEN_DIR "64x64_01.mtx", 64, 64, 4096 }, /**/
		{ MTX_GEN_DIR "64x64_02.mtx", 64, 64, 4096 }, /**/
		}, /**/


		{ /**/
		{ MTX_GEN_DIR "128x128_01.mtx", 128, 128, 16384 }, /**/
		{ MTX_GEN_DIR "128x128_02.mtx", 128, 128, 16384 }, /**/
		}, /**/

		{ /**/
		{ MTX_GEN_DIR "256x256_01.mtx", 256, 256, 65536 }, /**/
		{ MTX_GEN_DIR "256x256_02.mtx", 256, 256, 65536 }, /**/
		}, /**/

		{ /**/
		{ MTX_GEN_DIR "512x512_01.mtx", 512, 512, 262144 }, /**/
		{ MTX_GEN_DIR "512x512_02.mtx", 512, 512, 262144 }, /**/
		}, /**/

		{ { { 0 } } }, /**/

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
