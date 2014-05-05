/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
 * https://github.com/nesro/sparse-matrices
 */

#include <stdlib.h>
#include <stdio.h>

#include "virtual_matrix.h"
#include "utils.h"
#include "mmio.h"
#include "mm_load.h"

mm_file_t *mm_load(const char *filename, int round_to_power_2) {

	mm_file_t *mm_file;

	MM_typecode matcode;
	FILE *f;
	int i;
	int j;
	int is_symmetric;
	int rounded_size;

	/* for sort */
	int swap_r;
	int swap_c;
	datatype_t swap_v;

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

	mm_file = malloc(sizeof(mm_file_t));

	if (mm_read_mtx_crd_size(f, &mm_file->height, &mm_file->width,
			&mm_file->nnz) != 0) {
		die("wrong mm_read_mtx_crd_size");
	}

	if (round_to_power_2) {
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

	if (is_symmetric)
		mm_file->data_size = mm_file->nnz * 2;
	else
		mm_file->data_size = mm_file->nnz;

	mm_file->data = malloc(mm_file->data_size * sizeof(mm_item_t));

	for (i = 0; i < mm_file->nnz; i++) {

		if (fscanf(f, "%d %d "DPF"\n", &mm_file->data[i].row,
				&mm_file->data[i].col, &mm_file->data[i].value) != 3) {
			fdie("fscanf failed for file %s and item no. %i\n", filename, i);
		}

		mm_file->data[i].row--;
		mm_file->data[i].col--;

		/* FIXME: all symetric, and sort */
		if (is_symmetric && mm_file->data[i].row < mm_file->data[i].col) {
			mm_file->data[i].value = 0;
		}

		if (is_symmetric && mm_file->data[i].row != mm_file->data[i].col) {
			i++;
			mm_file->data[i].row = mm_file->data[i - 1].col;
			mm_file->data[i].col = mm_file->data[i - 1].row;
			mm_file->data[i].value = mm_file->data[i - 1].value;
			mm_file->nnz++;
		}
	}

	/* TODO: better sort pls */
	/* bubble sort */
	if (is_symmetric) { /* we should sort every time */
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
