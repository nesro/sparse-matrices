#include <stdlib.h>
#include <stdio.h>

#include "virtual_matrix.h"
#include "utils.h"
#include "mmio.h"
#include "mm_load.h"

mm_file_t *mm_load(const char *filename) {

	mm_file_t *mm_file;

	MM_typecode matcode;
	FILE *f;
	int i;

	if ((f = fopen(filename, "r")) == NULL) {
		fdie("File %s doesn't exists. Exiting.\n", filename);
	}

	if (mm_read_banner(f, &matcode) != 0) {
		die("Could not process Matrix Market banner.\n");
	}

	if (mm_is_complex(matcode) && mm_is_matrix(matcode) && mm_is_sparse(matcode)) {
		fdie("Sorry, this application does not support Market Market type: "
			"[%s]\n", mm_typecode_to_str(matcode));
	}

	mm_file = malloc(sizeof(mm_file_t));

	if (mm_read_mtx_crd_size(f, &mm_file->height, &mm_file->width,
		&mm_file->items) != 0) {
		die("wrong mm_read_mtx_crd_size");
	}

	mm_file->data = malloc(mm_file->items * sizeof(mm_item_t));

	for (i = 0; i < mm_file->items; i++) {
		if (fscanf(f, "%d %d %lg\n", &mm_file->data[i].row,
			&mm_file->data[i].col, &mm_file->data[i].value) != 3) {
			fdie("fscanf failed for file %s and item no. %i\n", filename, i);
		}
	}

	if (f != stdin)
		fclose(f);

	return mm_file;
}

void mm_free(mm_file_t *mm_file) {
	if (mm_file != NULL)
		free(mm_file->data);

	free(mm_file);
}
