/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
 * https://github.com/nesro/sparse-matrices
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "virtual_matrix.h"
#include "utils.h"
#include "mmio.h"
#include "mm_load.h"

int cmp_item(const void * a, const void * b) {
	return (((mm_item_t*) a)->value - ((mm_item_t*) b)->value);
}

mm_file_t *mm_load(const char *filename, int round_to_power_2) {

	mm_file_t *mm_file;
	MM_typecode matcode;
	FILE *f;
	int i;

	int is_symmetric;
	int is_pattern;
	int rounded_size;

	if ((f = fopen(filename, "r")) == NULL) {
		fdie("File %s doesn't exists. Exiting.\n", filename);
	}

	if (mm_read_banner(f, &matcode) != 0) {
		die("Could not process Matrix Market banner.\n");
	}

	if (mm_is_complex(matcode) && mm_is_matrix(matcode) && mm_is_sparse(matcode)) {
		fdie(
				"Sorry, this application does not support Market Market type: " "[%s]\n",
				mm_typecode_to_str(matcode));
	}

	if (mm_is_symmetric(matcode))
		is_symmetric = 1;
	else
		is_symmetric = 0;

	if (mm_is_pattern(matcode))
		is_pattern = 1;
	else
		is_pattern = 0;

	mm_file = malloc(sizeof(mm_file_t));

	if (mm_read_mtx_crd_size(f, &mm_file->height, &mm_file->width,
			&mm_file->nnz) != 0) {
		die("wrong mm_read_mtx_crd_size");
	}

	/* TODO: for dense, coo and csr matrices, there is no reason
	 * to round them up. but this is for testing much easier */
	if (1 || round_to_power_2) {

		if (round_to_power_2 != 666) {

			rounded_size = maxi(mm_file->height, mm_file->width);
			rounded_size--;
			rounded_size |= rounded_size >> 1;
			rounded_size |= rounded_size >> 2;
			rounded_size |= rounded_size >> 4;
			rounded_size |= rounded_size >> 8;
			rounded_size |= rounded_size >> 16;
			rounded_size++;
			mm_file->height = rounded_size;
			mm_file->width = rounded_size;
		}
	}

	if (is_symmetric)
		mm_file->data_size = mm_file->nnz * 2;
	else
		mm_file->data_size = mm_file->nnz;

	mm_file->data = malloc(mm_file->data_size * sizeof(mm_item_t));

	for (i = 0; i < mm_file->nnz; i++) {
		if (is_pattern) {
			if (fscanf(f, "%d %d\n", &mm_file->data[i].row,
					&mm_file->data[i].col) != 2) {
				fdie("fscanf failed for file %s (pattern) and item no. %i\n",
						filename, i);
			}
			mm_file->data[i].value = 1;
		} else {
			if (fscanf(f, "%d %d "DPF"\n", &mm_file->data[i].row,
					&mm_file->data[i].col, &mm_file->data[i].value) != 3) {
				fdie("fscanf failed for file %s and item no. %i\n", filename,
						i);
			}
		}

		mm_file->data[i].row--;
		mm_file->data[i].col--;

		/* XXX: numerical stability? */
//		mm_file->data[i].value = floor(1.5*mm_file->data[i].value);
//		mm_file->data[i].value /= 1.5;
		/* FIXME: all symetric, and sort */
		if (is_symmetric && mm_file->data[i].row < mm_file->data[i].col) {
			fdie("Bad symmetric value at %d\n", i);
		}

		if (is_symmetric && mm_file->data[i].row != mm_file->data[i].col) {
			i++;
			mm_file->data[i].row = mm_file->data[i - 1].col;
			mm_file->data[i].col = mm_file->data[i - 1].row;
			mm_file->data[i].value = mm_file->data[i - 1].value;
			mm_file->nnz++;
		}
	}

	/* TODO: a better sort pls and maybe we should sort every time with O(n)
	 * check whether the data are sorted or not */
	/* WARN: this is O(n^2) */
	if (is_symmetric) {
#if 0
		/* bubble sort */
		int j;
		int swap_r;
		int swap_c;
		datatype_t swap_v;
		for (i = 0; i < (mm_file->nnz - 1); i++) {
			for (j = 0; j < mm_file->nnz - i - 1; j++) {
				if (mm_file->data[j].row > mm_file->data[j + 1].row) {

					swap_r = mm_file->data[j].row;
					swap_c = mm_file->data[j].col;
					swap_v = mm_file->data[j].value;

					mm_file->data[j].row = mm_file->data[j + 1].row;
					mm_file->data[j].col = mm_file->data[j + 1].col;
					mm_file->data[j].value = mm_file->data[j + 1].value;

					mm_file->data[j + 1].row = swap_r;
					mm_file->data[j + 1].col = swap_c;
					mm_file->data[j + 1].value = swap_v;
				}
			}
		}
#else
		qsort(mm_file->data, mm_file->nnz, sizeof(mm_item_t), cmp_item);
#endif
	}

#if 0

	for (i = 0; i < mm_file->nnz; i++) {
		printf("i=%d y=%d, x=%d, v=%lf\n", i, mm_file->data[i].row,
				mm_file->data[i].col, mm_file->data[i].value);
	}

#endif

	if (f != stdin)
		fclose(f);

	return mm_file;
}

void mm_free(mm_file_t *mm_file) {
	if (mm_file != NULL)
		free(mm_file->data);

	free(mm_file);
}
