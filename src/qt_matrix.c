/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
 * https://github.com/nesro/sparse-matrices
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "omp.h"
#include "qt_matrix.h"

void qt_matrix_init(qt_matrix_t *qt_matrix, int width, int height, int nnz,
	int sm_size) {

	qt_matrix->info.w = width;
	qt_matrix->info.h = height;
	qt_matrix->info.nnz = nnz;

	qt_matrix->root = calloc(1, sizeof(qt_node_t));
	assert(qt_matrix->root != NULL);

	qt_matrix->height = ceil(
		(log((width * height) / (sm_size * sm_size)) / (double) log(4)));

	qt_matrix->blocks = 0;

	qt_matrix->sm_size = sm_size;

	qt_matrix->v = NULL;
	qt_matrix->ci = NULL;
	qt_matrix->rp = NULL;

	/* for loading purposes only */
	qt_matrix->last_index_v = 0;
	qt_matrix->last_index_rp = 0;

	/* for information purposes only */
	qt_matrix->filename = NULL;
}

static void qt_matrix_free_node(qt_node_t *node) {
	if (node == NULL)
		return;

	free(node->sm);

	qt_matrix_free_node(node->tl);
	qt_matrix_free_node(node->tr);
	qt_matrix_free_node(node->bl);
	qt_matrix_free_node(node->br);

	free(node);
}

void qt_matrix_free(qt_matrix_t *qt_matrix) {
	free(qt_matrix->v);
	free(qt_matrix->ci);
	free(qt_matrix->rp);
	qt_matrix_free_node(qt_matrix->root);
}

static void inner_qt_matrix_to_dense(qt_matrix_t *qt_matrix, qt_node_t *node,
	dense_matrix_t *dense_matrix) {

	int i;
	int j;
	int rp;

	if (node->sm != NULL) {
		/* CSR -> DENSE */
		for (i = 0; i < qt_matrix->sm_size; i++) {
			rp = node->sm->irp + i;

			for (j = qt_matrix->rp[rp]; j < qt_matrix->rp[rp + 1]; j++) {
				dense_matrix->v[ /**/
				node->sm->y + i /**/
				][ /**/
				node->sm->x + qt_matrix->ci[node->sm->iv + j] /**/
				] += qt_matrix->v[node->sm->iv + j];
			}
		}

		return;
	}

	if (node->tl != NULL)
		inner_qt_matrix_to_dense(qt_matrix, node->tl, dense_matrix);
	if (node->tr != NULL)
		inner_qt_matrix_to_dense(qt_matrix, node->tr, dense_matrix);
	if (node->bl != NULL)
		inner_qt_matrix_to_dense(qt_matrix, node->bl, dense_matrix);
	if (node->br != NULL)
		inner_qt_matrix_to_dense(qt_matrix, node->br, dense_matrix);
}

void qt_matrix_to_dense(qt_matrix_t *qt_matrix, dense_matrix_t *dense_matrix) {
	dense_matrix_init(dense_matrix, qt_matrix->info.w, qt_matrix->info.h);
	inner_qt_matrix_to_dense(qt_matrix, qt_matrix->root, dense_matrix);
}

qt_submatrix_t *qt_submatrix_get(qt_matrix_t *qt_matrix, int item_row,
	int item_col) {

	int i;
	qt_node_t **tmp_node = &qt_matrix->root;

	int half_width = qt_matrix->info.w / 2;
	int half_height = qt_matrix->info.h / 2;
	int delkoef = qt_matrix->info.w / 2;

	for (i = 0;; i++) {
		if (*tmp_node == NULL) {
			*tmp_node = calloc(1, sizeof(qt_node_t));

			if (*tmp_node == NULL) {
				fprintf(stderr, "tmp_node == NULL\n");
				exit(1);
			}
		}

		if (i == qt_matrix->height) {
			break;
		}

		delkoef /= 2;

		if (item_row < half_height) {
			if (item_col < half_width) {
				tmp_node = &((*tmp_node)->tl);
				half_width -= delkoef;
				half_height -= delkoef;
			} else {
				tmp_node = &((*tmp_node)->tr);
				half_width += delkoef;
				half_height -= delkoef;
			}
		} else {
			if (item_col < half_width) {
				tmp_node = &((*tmp_node)->bl);
				half_width -= delkoef;
				half_height += delkoef;
			} else {
				tmp_node = &((*tmp_node)->br);
				half_width += delkoef;
				half_height += delkoef;
			}
		}
	}

	if ((*tmp_node)->sm == NULL) {
		(*tmp_node)->sm = calloc(1, sizeof(qt_submatrix_t));
		assert((*tmp_node)->sm != NULL);

		qt_matrix->blocks++;
		(*tmp_node)->sm->x = (item_col / qt_matrix->sm_size)
			* qt_matrix->sm_size;
		(*tmp_node)->sm->y = (item_row / qt_matrix->sm_size)
			* qt_matrix->sm_size;

	}

	assert((*tmp_node)->sm->x <= item_col);
	assert(item_col < (*tmp_node)->sm->x + qt_matrix->sm_size);
	assert((*tmp_node)->sm->y <= item_row);
	assert(item_row < (*tmp_node)->sm->y + qt_matrix->sm_size);

	return (*tmp_node)->sm;
}

/**
 * @param ma
 * @param mb
 * @param a
 * @param b
 * @param c
 */
static void qt_node_mul(qt_matrix_t *ma, qt_matrix_t *mb, qt_node_t *a,
	qt_node_t *b, dense_matrix_t *c) {

	if (a == NULL || b == NULL)
		return;

	if (a->sm != NULL) {
		assert(b->sm != NULL);

		int r;
		int rpa;
		int rpb;

		int j;
		int k;

		/*
		 * I want to optimize this part a little bit.
		 * The good readable part is in the if 0 block.
		 */
#if 0
		for (r = 0; r < ma->sm_size; r++) {
			rpa = a->sm->irp + r;
			for (j = ma->rp[rpa]; j < ma->rp[rpa + 1]; j++) {
				rpb = b->sm->irp + ma->ci[b->sm->iv + j];
				for (k = mb->rp[rpb]; k < mb->rp[rpb + 1]; k++) {

					c->v[
					a->sm->y + r /**/
					][ /**/
					b->sm->x + mb->ci[b->sm->iv + k] /**/
					] += ma->v[a->sm->iv + j] * mb->v[b->sm->iv + k];
				}
			}
		}
#else
		datatype_t *cptr;

		for (r = 0; r < ma->sm_size; r++) {

			rpa = a->sm->irp + r;
			cptr = &c->v[a->sm->y + r][b->sm->x];

			for (j = a->sm->iv + ma->rp[rpa]; j < a->sm->iv + ma->rp[rpa + 1];
				j++) {

				rpb = b->sm->irp + ma->ci[j];

				for (k = b->sm->iv + mb->rp[rpb];
					k < b->sm->iv + mb->rp[rpb + 1]; k++) {

					cptr[mb->ci[k]] += ma->v[j] * mb->v[k];
				}
			}
		}

//		datatype_t **cptr = NULL;
//
//		for (r = 0; r < ma->sm_size; r++) {
//			**cptr = c->v[a->sm->y + r][b->sm->x];
//			rpa = a->sm->irp + r;
//			for (j = ma->rp[rpa]; j < ma->rp[rpa + 1]; j++) {
//				rpb = b->sm->irp + ma->ci[b->sm->iv + j];
//				for (k = b->sm->iv + mb->rp[rpb];
//					k < b->sm->iv + mb->rp[rpb + 1]; k++) {
//
//					**(cptr + mb->ci[k]) += ma->v[a->sm->iv + j] * mb->v[+k];
//				}
//			}
//		}
#endif

		return;
	}

	/*
	 |A|B| * |E|F| = |A*E + B*G|A*F + B*H|
	 |C|D|   |G|H|   |C*E + D*G|C*F + D*H|
	 */

	qt_node_mul(ma, mb, a->tl, b->tl, c); /* AE */
	qt_node_mul(ma, mb, a->tr, b->bl, c); /* BG */

	qt_node_mul(ma, mb, a->tl, b->tr, c); /* AF */
	qt_node_mul(ma, mb, a->tr, b->br, c); /* BH */

	qt_node_mul(ma, mb, a->bl, b->tl, c); /* CE */
	qt_node_mul(ma, mb, a->br, b->bl, c); /* DG */

	qt_node_mul(ma, mb, a->bl, b->tr, c); /* CF */
	qt_node_mul(ma, mb, a->br, b->br, c); /* DH */
}

double qt_matrix_matrix_mul(qt_matrix_t *a, qt_matrix_t *b, dense_matrix_t *c) {

	double start_time;
	double end_time;

	assert(a->info.w == b->info.w);
	assert(a->info.h == b->info.h);

	dense_matrix_init(c, a->info.w, a->info.h);

	start_time = omp_get_wtime();
	qt_node_mul(a, b, a->root, b->root, c);
	end_time = omp_get_wtime();

	return end_time - start_time;
}

/******************************************************************************/

/******************************************************************************/

/******************************************************************************/

static void compute_blocks_beginnings(qt_matrix_t *qt_matrix, qt_node_t *node) {

	if (node->sm != NULL) {
		node->sm->iv = qt_matrix->last_index_v;
		node->sm->irp = qt_matrix->last_index_rp;

		qt_matrix->last_index_v += node->sm->nnz;
		qt_matrix->last_index_rp += qt_matrix->sm_size + 1;

		return;
	}

	if (node->tl != NULL)
		compute_blocks_beginnings(qt_matrix, node->tl);
	if (node->tr != NULL)
		compute_blocks_beginnings(qt_matrix, node->tr);
	if (node->bl != NULL)
		compute_blocks_beginnings(qt_matrix, node->bl);
	if (node->br != NULL)
		compute_blocks_beginnings(qt_matrix, node->br);
}

static void compute_blocks_rp(qt_matrix_t *qt_matrix, qt_node_t *node) {

	int i;
	int sum;
	int tmp;

	if (node->sm != NULL) {
		sum = 0;
		for (i = 0; i < qt_matrix->sm_size; i++) {
			tmp = qt_matrix->rp[node->sm->irp + i];
			qt_matrix->rp[node->sm->irp + i] = sum;
			sum += tmp;
		}

		qt_matrix->rp[node->sm->irp + qt_matrix->sm_size] = node->sm->nnz;

		return;
	}

	if (node->tl != NULL)
		compute_blocks_rp(qt_matrix, node->tl);
	if (node->tr != NULL)
		compute_blocks_rp(qt_matrix, node->tr);
	if (node->bl != NULL)
		compute_blocks_rp(qt_matrix, node->bl);
	if (node->br != NULL)
		compute_blocks_rp(qt_matrix, node->br);
}

double qt_matrix_load_mm(qt_matrix_t *qt_matrix, const char *filename,
	int sm_size) {

	double start_time;
	double end_time;

	int i;

	MM_typecode matcode; /* type of matrix in MM format */
	FILE *f; /* file pointer */
	int h; /* height */
	int w; /* width */
	int nnz; /* number of nonzero elements */

	int *coo_x; /* coo matrix, x coordinates */
	int *coo_y; /* coo matrix, y coordinates */
	datatype_t *coo_v; /* coo matrix, values */

	qt_submatrix_t *tmp_sm; /* temporal submatrix */

	/**************************************************************************/

	if ((f = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "File %s doesn't exists. Exiting.\n", filename);
		exit(1);
	}

	if (mm_read_banner(f, &matcode) != 0) {
		printf("Could not process Matrix Market banner.\n");
		exit(1);
	}

	if (mm_is_complex(matcode) && mm_is_matrix(matcode) && mm_is_sparse(matcode)) {
		printf("Sorry, this application does not support "
			"Market Market type: [%s]\n", mm_typecode_to_str(matcode));
		exit(1);
	}

	if (mm_read_mtx_crd_size(f, &h, &w, &nnz) != 0) {
		perror("mm_read_mtx_crd_size");
		exit(1);
	}

	/* XXX: */
	assert(h == w);
	assert(is_power_of_two(h));
	assert(nnz > 0);

	qt_matrix_init(qt_matrix, w, h, nnz, sm_size);
	qt_matrix->filename = filename;

	coo_x = malloc(nnz * sizeof(int));
	coo_y = malloc(nnz * sizeof(int));
	coo_v = malloc(nnz * sizeof(datatype_t));
	assert(coo_x != NULL && coo_y != NULL && coo_v != NULL);

	start_time = omp_get_wtime();

	for (i = 0; i < nnz; i++) {

		if (fscanf(f, "%d %d %lg\n", &coo_y[i], &coo_x[i], &coo_v[i]) != 3) {
			perror("fscanf");
			exit(1);
		}

		coo_x[i]--;
		coo_y[i]--;

		tmp_sm = qt_submatrix_get(qt_matrix, coo_y[i], coo_x[i]);
		tmp_sm->nnz++;
	}
	fclose(f);

	compute_blocks_beginnings(qt_matrix, qt_matrix->root);

	qt_matrix->v = calloc(nnz, sizeof(datatype_t));
	qt_matrix->ci = calloc(nnz, sizeof(int));

	qt_matrix->rp = calloc(qt_matrix->blocks * (qt_matrix->sm_size + 1),
		sizeof(int));

	for (i = 0; i < nnz; i++) {
		tmp_sm = qt_submatrix_get(qt_matrix, coo_y[i], coo_x[i]);

		qt_matrix->v[tmp_sm->iv + tmp_sm->i] = coo_v[i];
		qt_matrix->ci[tmp_sm->iv + tmp_sm->i] = coo_x[i] % qt_matrix->sm_size;
		qt_matrix->rp[tmp_sm->irp + (coo_y[i] % qt_matrix->sm_size)]++;
		tmp_sm->i++;
	}

	compute_blocks_rp(qt_matrix, qt_matrix->root);

	free(coo_x);
	free(coo_y);
	free(coo_v);

	end_time = omp_get_wtime();
	return end_time - start_time;

}

static void qt_matrix_print_inner(qt_matrix_t *qt_matrix, qt_node_t *node,
	int depth) {

	depth += 8;

	if (node->sm != NULL) {
		assert(node->tl == NULL);
		assert(node->tr == NULL);
		assert(node->bl == NULL);
		assert(node->br == NULL);

		printf("%*ssubmatrix:\n", depth, "");
		printf("%*s> x   = %d\n", depth, "", node->sm->x);
		printf("%*s> y   = %d\n", depth, "", node->sm->y);
		printf("%*s> nnz = %d\n", depth, "", node->sm->nnz);
		printf("%*s> iv  = %d\n", depth, "", node->sm->iv);
		printf("%*s> irp = %d\n", depth, "", node->sm->irp);

		return;
	}

	if (node->tl != NULL) {
		assert(node->sm == NULL);

		printf("%*s\n", depth, "tl");
		qt_matrix_print_inner(qt_matrix, node->tl, depth);
	}

	if (node->tr != NULL) {
		assert(node->sm == NULL);

		printf("%*s\n", depth, "tr");
		qt_matrix_print_inner(qt_matrix, node->tr, depth);
	}

	if (node->bl != NULL) {
		assert(node->sm == NULL);

		printf("%*s\n", depth, "bl");
		qt_matrix_print_inner(qt_matrix, node->bl, depth);
	}

	if (node->br != NULL) {
		assert(node->sm == NULL);

		printf("%*s\n", depth, "br");
		qt_matrix_print_inner(qt_matrix, node->br, depth);
	}
}

void qt_matrix_print(qt_matrix_t *qt_matrix) {
	int i;

	if (!PRINT)
		return;

	if (qt_matrix == NULL) {
		fprintf(stderr, "qt_matrix == NULL\n");
		return;
	}

	printf("matrix from file: %s\n",
		qt_matrix->filename != NULL ? qt_matrix->filename : "NULL");

	printf("sm_size: %d\n", qt_matrix->sm_size);
	printf("quadtree height = %d\n", qt_matrix->height);

	printf("v:  ");
	for (i = 0; i < qt_matrix->info.nnz; i++) {
		printf("%02.lf,", qt_matrix->v[i]);
		if ((i + 1) % qt_matrix->sm_size == 0)
			printf(" _ ");
	}
	printf("\nci: ");
	for (i = 0; i < qt_matrix->info.nnz; i++) {
		printf("%02d,", qt_matrix->ci[i]);
		if ((i + 1) % qt_matrix->sm_size == 0)
			printf(" _ ");
	}
	printf("\nrp: ");
	for (i = 0; i < qt_matrix->blocks * (qt_matrix->sm_size + 1); i++) {
		printf("%d,", qt_matrix->rp[i]);
		if ((i + 1) % (qt_matrix->sm_size + 1) == 0)
			printf(" _ ");
	}

	printf("nodes:\n");
	qt_matrix_print_inner(qt_matrix, qt_matrix->root, 0);

	printf("\n");
}

time_record_t qt_matrix_mm_mul(const char *matrix_a, const char *matrix_b, int sm_size) {

	qt_matrix_t a;
	qt_matrix_t b;
	dense_matrix_t c;
	time_record_t tr;

	tr.load_a = qt_matrix_load_mm(&a, matrix_a, sm_size);
	tr.load_b = qt_matrix_load_mm(&b, matrix_b, sm_size);

	tr.multiplication = qt_matrix_matrix_mul(&a, &b, &c);

	qt_matrix_free(&a);
	qt_matrix_free(&b);
	dense_matrix_free(&c);

	return tr;
}

