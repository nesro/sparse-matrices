/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013, 2014
 * https://github.com/nesro/sparse-matrices
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "omp.h"
#include "virtual_matrix.h"
#include "qdt_matrix.h"
#include "den_matrix.h"

static vm_vmt_t qdt_vmt = { /**/
(reset_t) NULL, /**/
(free_t) qdt_free, /**/
(mm_load_t) NULL,/**/
(mm_save_t) NULL, /**/
(print_t) NULL, /**/
(compare_t) NULL, /**/
(distance_t) NULL, /**/
(convert_t) NULL, /**/
(mul_t) qdt_mul, /**/
};

void qdt_vm_init(qdt_matrix_t **qdt, va_list va) {

	int va_flag = 0;
	int width;
	int height;
	int nnz;
	int sm_size;

	sm_size = va_get_int(va, 1, &va_flag);
	nnz = va_get_int(va, 1, &va_flag);
	height = va_get_int(va, -1, &va_flag);
	width = va_get_int(va, -1, &va_flag);

	qdt_init(qdt, width, height, nnz, sm_size);
}

void qdt_init(qdt_matrix_t **qdt, int width, int height, int nnz, int sm_size) {

	*qdt = malloc(sizeof(qdt_matrix_t));
	assert(*qdt != NULL);

	(*qdt)->_.type = QDT;
	(*qdt)->_.f = qdt_vmt;
	(*qdt)->_.w = width;
	(*qdt)->_.h = height;

	(*qdt)->_.nnz = nnz;

	(*qdt)->root = calloc(1, sizeof(qdt_node_t));
	assert((*qdt)->root != NULL);

	(*qdt)->height = ceil(
			(log((width * height) / (sm_size * sm_size)) / (double) log(4)));

	(*qdt)->blocks = 0;
	(*qdt)->sm_size = sm_size;

#if QDT_CSR
	(*qdt)->v = NULL;
	(*qdt)->ci = NULL;
	(*qdt)->rp = NULL;

	/* for loading purposes only */
	(*qdt)->last_index_v = 0;
	(*qdt)->last_index_rp = 0;
#endif
}

static void qdt_free_node(qdt_node_t *node) {

	if (node == NULL)
		return;

	free(node->sm);

	qdt_free_node(node->tl);
	qdt_free_node(node->tr);
	qdt_free_node(node->bl);
	qdt_free_node(node->br);

	free(node);
}

void qdt_free(qdt_matrix_t *qdt_matrix) {
#if QDT_CSR
	free(qdt_matrix->v);
	free(qdt_matrix->ci);
	free(qdt_matrix->rp);
#endif /* QDT_CSR */

	qdt_free_node(qdt_matrix->root);
}

/******************************************************************************/

static void inner_qt_matrix_to_dense(qdt_matrix_t *qt_matrix, qdt_node_t *node,
		den_matrix_t *dense_matrix) {

	int i;
	int j;
	int rp;

	if (node->sm != NULL) {
#if QDT_CSR
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
#else /* QDT_CSR */
		den_offset_addto(
				(dense_matrix_t *) node->sm->m.f.convert(node->sm->m, DEN),
				dense_matrix, 0, 0, node->sm->y, node->sm->x);
#endif /* QDT_CSR */

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

void qdt_to_dense(qdt_matrix_t *qt_matrix, den_matrix_t *dense_matrix) {
//	FIXME
//	 den_matrix_init(&dense_matrix, (den_matrix_init_t ) { qt_matrix->info.w,
//			qt_matrix->info.h, 1 });
	inner_qt_matrix_to_dense(qt_matrix, qt_matrix->root, dense_matrix);
}

qdt_submatrix_t *qdt_submatrix_get(qdt_matrix_t *qt_matrix, int item_row,
		int item_col) {

	int i;
	qdt_node_t **tmp_node = &qt_matrix->root;

	int half_width = qt_matrix->_.w / 2;
	int half_height = qt_matrix->_.h / 2;
	int delkoef = qt_matrix->_.w / 2;

	for (i = 0;; i++) {
		if (*tmp_node == NULL) {
			*tmp_node = calloc(1, sizeof(qdt_node_t));

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
		(*tmp_node)->sm = calloc(1, sizeof(qdt_submatrix_t));
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
static void qdt_node_mul(qdt_matrix_t *ma, qdt_matrix_t *mb, qdt_node_t *a,
		qdt_node_t *b, den_matrix_t **c) {

	if (a == NULL || b == NULL)
		return;

	if (a->sm != NULL) {
		assert(b->sm != NULL);

		int r;
		int rpa;
		int rpb;

		int j;
		int k;

#if QDT_CSR
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
			cptr = &(*c)->v[a->sm->y + r][b->sm->x];

			for (j = a->sm->iv + ma->rp[rpa]; j < a->sm->iv + ma->rp[rpa + 1];
					j++) {

				rpb = b->sm->irp + ma->ci[j];

				for (k = b->sm->iv + mb->rp[rpb];
						k < b->sm->iv + mb->rp[rpb + 1]; k++) {

					cptr[mb->ci[k]] += ma->v[j] * mb->v[k];
				}
			}
		}

#endif /* QDT_CSR */

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

	qdt_node_mul(ma, mb, a->tl, b->tl, c); /* AE */
	qdt_node_mul(ma, mb, a->tr, b->bl, c); /* BG */

	qdt_node_mul(ma, mb, a->tl, b->tr, c); /* AF */
	qdt_node_mul(ma, mb, a->tr, b->br, c); /* BH */

	qdt_node_mul(ma, mb, a->bl, b->tl, c); /* CE */
	qdt_node_mul(ma, mb, a->br, b->bl, c); /* DG */

	qdt_node_mul(ma, mb, a->bl, b->tr, c); /* CF */
	qdt_node_mul(ma, mb, a->br, b->br, c); /* DH */
}

double qdt_mul(qdt_matrix_t *a, qdt_matrix_t *b, den_matrix_t **c,
		char flag /* unused */) {

	double start_time;
	double end_time;

	assert(a->_.w == b->_.w);
	assert(a->_.h == b->_.h);

	if (*c == NULL)
		vm_create((vm_t **) c, DEN, 1, a->_.w, a->_.w);

	start_time = omp_get_wtime();
	qdt_node_mul(a, b, a->root, b->root, c);
	end_time = omp_get_wtime();

	return end_time - start_time;
}

/******************************************************************************/

/******************************************************************************/

/******************************************************************************/

static void compute_blocks_beginnings(qdt_matrix_t *qdt_matrix,
		qdt_node_t *node) {

	if (node->sm != NULL) {
		node->sm->iv = qdt_matrix->last_index_v;
		node->sm->irp = qdt_matrix->last_index_rp;

		qdt_matrix->last_index_v += node->sm->nnz;
		qdt_matrix->last_index_rp += qdt_matrix->sm_size + 1;

		return;
	}

	if (node->tl != NULL)
		compute_blocks_beginnings(qdt_matrix, node->tl);
	if (node->tr != NULL)
		compute_blocks_beginnings(qdt_matrix, node->tr);
	if (node->bl != NULL)
		compute_blocks_beginnings(qdt_matrix, node->bl);
	if (node->br != NULL)
		compute_blocks_beginnings(qdt_matrix, node->br);
}

static void compute_blocks_rp(qdt_matrix_t *qdt_matrix, qdt_node_t *node) {

	int i;
	int sum;
	int tmp;

	if (node->sm != NULL) {
		sum = 0;
		for (i = 0; i < qdt_matrix->sm_size; i++) {
			tmp = qdt_matrix->rp[node->sm->irp + i];
			qdt_matrix->rp[node->sm->irp + i] = sum;
			sum += tmp;
		}

		qdt_matrix->rp[node->sm->irp + qdt_matrix->sm_size] = node->sm->nnz;

		return;
	}

	if (node->tl != NULL)
		compute_blocks_rp(qdt_matrix, node->tl);
	if (node->tr != NULL)
		compute_blocks_rp(qdt_matrix, node->tr);
	if (node->bl != NULL)
		compute_blocks_rp(qdt_matrix, node->bl);
	if (node->br != NULL)
		compute_blocks_rp(qdt_matrix, node->br);
}

/******************************************************************************/

void qdt_from_mm(qdt_matrix_t **qdt, const char *file, va_list va) {

	int va_flag = 0;
	int sm_size;

	sm_size = va_get_int(va, -1, &va_flag);

	assert(sm_size != -1);

	qdt_load_mm(qdt, file, sm_size);
}

double qdt_load_mm(qdt_matrix_t **qdt_matrix, const char *filename, int sm_size) {

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

	qdt_submatrix_t *tmp_sm; /* temporal submatrix */

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

	qdt_init(qdt_matrix, w, h, nnz, sm_size);
//	qdt_matrix->filename = filename;

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

		tmp_sm = qdt_submatrix_get(*qdt_matrix, coo_y[i], coo_x[i]);
		tmp_sm->nnz++;
	}
	fclose(f);

	compute_blocks_beginnings(*qdt_matrix, (*qdt_matrix)->root);

	(*qdt_matrix)->v = calloc(nnz, sizeof(datatype_t));
	(*qdt_matrix)->ci = calloc(nnz, sizeof(int));

	(*qdt_matrix)->rp = calloc(
			(*qdt_matrix)->blocks * ((*qdt_matrix)->sm_size + 1), sizeof(int));

	for (i = 0; i < nnz; i++) {
		tmp_sm = qdt_submatrix_get(*qdt_matrix, coo_y[i], coo_x[i]);

		(*qdt_matrix)->v[tmp_sm->iv + tmp_sm->i] = coo_v[i];
		(*qdt_matrix)->ci[tmp_sm->iv + tmp_sm->i] = coo_x[i]
				% (*qdt_matrix)->sm_size;
		(*qdt_matrix)->rp[tmp_sm->irp + (coo_y[i] % (*qdt_matrix)->sm_size)]++;
		tmp_sm->i++;
	}

	compute_blocks_rp(*qdt_matrix, (*qdt_matrix)->root);

	free(coo_x);
	free(coo_y);
	free(coo_v);

	end_time = omp_get_wtime();
	return end_time - start_time;

}

static void qt_matrix_print_inner(qdt_matrix_t *qdt_matrix, qdt_node_t *node,
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
		qt_matrix_print_inner(qdt_matrix, node->tl, depth);
	}

	if (node->tr != NULL) {
		assert(node->sm == NULL);

		printf("%*s\n", depth, "tr");
		qt_matrix_print_inner(qdt_matrix, node->tr, depth);
	}

	if (node->bl != NULL) {
		assert(node->sm == NULL);

		printf("%*s\n", depth, "bl");
		qt_matrix_print_inner(qdt_matrix, node->bl, depth);
	}

	if (node->br != NULL) {
		assert(node->sm == NULL);

		printf("%*s\n", depth, "br");
		qt_matrix_print_inner(qdt_matrix, node->br, depth);
	}
}

void qt_matrix_print(qdt_matrix_t *qdt_matrix) {
	int i;

	if (!PRINT)
		return;

	if (qdt_matrix == NULL) {
		fprintf(stderr, "qt_matrix == NULL\n");
		return;
	}

	printf("matrix from file: %s\n",
			qdt_matrix->filename != NULL ? qdt_matrix->filename : "NULL");

	printf("sm_size: %d\n", qdt_matrix->sm_size);
	printf("quadtree height = %d\n", qdt_matrix->height);

	printf("v:  ");
	for (i = 0; i < qdt_matrix->_.nnz; i++) {
		printf("%02.lf,", qdt_matrix->v[i]);
		if ((i + 1) % qdt_matrix->sm_size == 0)
			printf(" _ ");
	}
	printf("\nci: ");
	for (i = 0; i < qdt_matrix->_.nnz; i++) {
		printf("%02d,", qdt_matrix->ci[i]);
		if ((i + 1) % qdt_matrix->sm_size == 0)
			printf(" _ ");
	}
	printf("\nrp: ");
	for (i = 0; i < qdt_matrix->blocks * (qdt_matrix->sm_size + 1); i++) {
		printf("%d,", qdt_matrix->rp[i]);
		if ((i + 1) % (qdt_matrix->sm_size + 1) == 0)
			printf(" _ ");
	}

	printf("nodes:\n");
	qt_matrix_print_inner(qdt_matrix, qdt_matrix->root, 0);

	printf("\n");
}

//time_record_t qt_matrix_mm_mul(const char *matrix_a, const char *matrix_b,
//		int sm_size) {
//
//	qdt_matrix_t a;
//	qdt_matrix_t b;
//	den_matrix_t c;
//	time_record_t tr;
//
//	tr.load_a = qdt_load_mm(&a, matrix_a, sm_size);
//	tr.load_b = qdt_load_mm(&b, matrix_b, sm_size);
//
//	tr.multiplication = qdt_mul(&a, &b, &c);
//
//	qdt_free(&a);
//	qdt_free(&b);
//	den_matrix_free(&c);
//
//	return tr;
//}

