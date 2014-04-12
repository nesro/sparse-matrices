/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
 * https://github.com/nesro/sparse-matrices
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

#include "mm_load.h"
#include "virtual_matrix.h"
#include "kat_matrix.h"
#include "den_matrix.h"

static vm_vmt_t kat_vmt = { /**/
(reset_t) NULL, /**/
(free_t) kat_free, /**/
(mm_load_t) NULL,/**/
(mm_save_t) NULL, /**/
(print_t) NULL, /**/
(compare_t) NULL, /**/
(distance_t) NULL, /**/
(convert_t) NULL, /**/
(mul_t) NULL, /**/
};

void kat_print_node(kat_node_t *kat_node, int depth) {

	int i;
	int j;

	if (kat_node->node_type != INNER) {
		printf("%*s END \n", depth, "");
		return;
	}

	for (i = 0; i < KAT_N; i++) {
		for (j = 0; j < KAT_N; j++) {
			printf("%*s[%d][%d] = ", depth, "", i, j);

			if (kat_node->node.knp[i][j] != NULL) {
				printf(" !\n");
				kat_print_node(kat_node->node.knp[i][j], depth + 1);
			} else {
				printf("NULL\n");
			}
		}
	}

}

/**
 * kat_vm_init(vm_t **,
 */
void kat_vm_init(kat_matrix_t **kat, va_list va) {

	int va_flag = 0;
	int width;
	int height;
	int nnz;
	int sm_size;

	sm_size = va_get_int(va, 1, &va_flag);
	nnz = va_get_int(va, 1, &va_flag);
	height = va_get_int(va, -1, &va_flag);
	width = va_get_int(va, -1, &va_flag);

	kat_init(kat, width, height, nnz, sm_size);
}

void kat_init(kat_matrix_t **kat, int width, int height, int nnz, int sm_size) {

	/* TODO: add zeros to the matrix achieve this */
	assert(width == height);
	assert(is_power_of_two(width));

	assert(is_power_of_two(sm_size));
	assert(width*height > sm_size*KAT_K);

	*kat = calloc(1, sizeof(kat_matrix_t));
	assert(*kat != NULL);

	(*kat)->_.type = KAT;
	(*kat)->_.f = kat_vmt;
	(*kat)->_.w = width;
	(*kat)->_.h = height;

	(*kat)->sm_size = sm_size;
	(*kat)->_.nnz = nnz;

	(*kat)->root = calloc(1, sizeof(kat_node_t));
	assert((*kat)->root != NULL);
	(*kat)->root->node_type = INNER;

	(*kat)->height =
			ceil(
					(log((width * height) / (sm_size * sm_size))
							/ (double) log(KAT_K)));

	_s_debugf(KAT_DEBUG,
			"initializing kat_matrix_t with: n=%d, nnz=%d, sm_size=%d, height=%d\n",
			width, nnz, sm_size, (*kat)->height);

}

/******************************************************************************/

static void kat_node_free(kat_node_t *kat_node) {

	int i;
	int j;

	if (kat_node->node_type == INNER)
		for (i = 0; i < KAT_N; i++)
			for (j = 0; j < KAT_N; j++)
				if (kat_node->node.knp[i][j])
					kat_node_free(kat_node->node.knp[i][j]);

	free(kat_node);
}

void kat_free(kat_matrix_t *kat) {

	kat_node_free(kat->root);
	free(kat->v);

#if KAT_CSR
	free(kat->ci);
	free(kat->rp);
#endif /* KAT_CSR */

	free(kat);
}

/******************************************************************************/

static kat_node_t *kat_get_node(kat_matrix_t *kat, int y, int x) {

	int i;
	kat_node_t **tmp_node = &kat->root;
	int block_size = kat->_.w / KAT_N;
	int by = 0; /* block */
	int bx = 0;
	int oy; /* offset */
	int ox;

	_s_debugf(KAT_DEBUG, "Searching for node for [y,x]=[%d,%d]\n", y, x);

	for (i = 0; i < kat->height; i++) {

		oy = (y / block_size);
		ox = (x / block_size);
		y -= oy * block_size;
		x -= ox * block_size;
		by += oy;
		bx += ox;

		if (*tmp_node == NULL) {
			_s_debug(KAT_DEBUG, "creating a new inner node!\n");

			*tmp_node = calloc(1, sizeof(kat_node_t));
			assert(tmp_node != NULL);

			(*tmp_node)->node_type = INNER;
		}

		tmp_node = &((*tmp_node)->node.knp[oy][ox]);
		block_size /= KAT_N;

		_s_debugf(KAT_DEBUG, "next node is at oy=%d ox=%d bs=%d\n", oy, ox,
				block_size);
	}

	if (*tmp_node == NULL) {

		_s_debug(KAT_DEBUG, "creating a leaf node\n");

		*tmp_node = calloc(1, sizeof(kat_node_t));
		assert(tmp_node != NULL);

		(*tmp_node)->node_type = LEAF;
		(*tmp_node)->node.sm.y = by * kat->sm_size;
		(*tmp_node)->node.sm.x = bx * kat->sm_size;
	}

	_s_debugf(KAT_DEBUG, "Node found. by=%d bx=%d y=%d x=%d\n",
			(*tmp_node)->node.sm.y, (*tmp_node)->node.sm.x,
			(*tmp_node)->node.sm.y, (*tmp_node)->node.sm.x);

	return *tmp_node;
}

/*
 * At this point, all leaf nodes has their nnz variable final.
 */
void kat_determine_node(kat_matrix_t *kat, kat_node_t *kat_node) {

	int i;
	int j;

	_s_debugf(KAT_DEBUG, "search %d\n", kat_node->node_type);

	if (kat_node->node_type == LEAF) {

		if (1||((kat->sm_size * kat->sm_size) / (double) kat_node->node.sm.nnz)
				> KAT_DENSE_TRESHOLD) {
			_s_debugf(KAT_DEBUG, "node at y=%d x=%d is dense\n",
					kat_node->node.sm.y, kat_node->node.sm.x);
			kat_node->node_type = KAT_N_DEN;
			kat->den_blocks++;
			kat->den_blocks_nnz += kat_node->node.sm.n_nnz;
		} else {
			_s_debug(KAT_DEBUG, "node is csr\n");
			kat_node->node_type = KAT_N_CSR;
			kat->csr_blocks++;
		}

		return;
	}

	for (i = 0; i < KAT_N; i++)
		for (j = 0; j < KAT_N; j++)
			if (kat_node->node.knp[i][j])
				kat_determine_node(kat, kat_node->node.knp[i][j]);
			else
				_s_debugf(KAT_DEBUG, "y=%d x=%d is null\n", i, j);

}

/*
 * At this point, all leaf nodes has their nnz variable final and determined
 * type.
 */
__attribute__((optimize("unroll-loops")))
static void kat_prepare(kat_matrix_t *kat, kat_node_t *kat_node) {

	int i;
	int j;

	if (kat_node->node_type == KAT_N_DEN || kat_node->node_type == KAT_N_CSR) {
		switch (kat_node->node_type) {
		case KAT_N_DEN:
			_s_debugf(KAT_DEBUG, "setting up dense node %x\n", 0);
			_s_debugf(KAT_DEBUG, "(*kat)->last_v = %p\n", kat->last_v);

			kat_node->node.sm.v = kat->last_v;
			kat->last_v = &(kat->last_v[kat->sm_size * kat->sm_size]);

			_s_debugf(KAT_DEBUG, "NEW: (*kat)->last_v = %p\n", kat->last_v);
			break;
		case KAT_N_CSR:
//			/*FIXME*/
//			(*kat)->csr_blocks++;
//			(*kat)->last_v += kat_node->node.sm.nnz;
//			(*kat)->last_index_ci += kat_node->node.sm.nnz;
//			(*kat)->last_index_rp += (*kat)->sm_size + 1;
			break;
		case UNDEF:
		case INNER:
		case LEAF:
		default:
			assert(0);
			break;
		}

		return;
	}

	for (i = 0; i < KAT_N; i++)
		for (j = 0; j < KAT_N; j++)
			if (kat_node->node.knp[i][j])
				kat_prepare(kat, kat_node->node.knp[i][j]);
}

/*
 * Compute rp of csr nodes.
 */
__attribute__((optimize("unroll-loops")))
static void kat_csr_rp(kat_matrix_t *kat, kat_node_t *kat_node) {

	int i;
	int j;
	int tmp;
	int sum;

	if (kat_node->node_type == KAT_N_CSR) {
		sum = 0;

		for (i = 0; i < kat->sm_size; i++) {
			tmp = kat_node->node.sm.s.csr.rp[i];
			kat_node->node.sm.s.csr.rp[i] = sum;
			sum += tmp;
		}

		kat_node->node.sm.s.csr.rp[kat->sm_size] = kat_node->node.sm.nnz;

		return;
	}

	for (i = 0; i < KAT_N; i++)
		for (j = 0; j < KAT_N; j++)
			if (kat_node->node.knp[i][j]
					&& kat_node->node.knp[i][j]->node_type == INNER)
				kat_csr_rp(kat, kat_node->node.knp[i][j]);
}

/******************************************************************************/

void kat_from_mm(kat_matrix_t **kat, const char *file, va_list va) {

	int va_flag = 0;
	int sm_size;

	sm_size = va_get_int(va, -1, &va_flag);

	assert(sm_size != -1);

	kat_load_mm(kat, file, sm_size);
}

double kat_load_mm(kat_matrix_t **kat, const char *filename, int sm_size) {

	double start_time;
	double end_time;
	mm_file_t *mm_file;
	kat_node_t *tmp_node;
	int i;

	mm_file = mm_load(filename);

	/*
	 * Start timer after mm_file is loaded. We don't want to IO operations
	 * in our timer.
	 */
	start_time = omp_get_wtime();

	kat_init(kat, mm_file->width, mm_file->height, mm_file->nnz, sm_size);

	for (i = 0; i < mm_file->nnz; i++) {

		_s_debugf(KAT_DEBUG, "i=%d v=%lf\n", i, mm_file->data[i].value);

		tmp_node = kat_get_node(*kat, mm_file->data[i].row,
				mm_file->data[i].col);
		tmp_node->node.sm.nnz++;
	}

	/*
	 * Allocating memory for matrix data.
	 */
	/*  - (*kat)->den_blocks_nnz */
	(*kat)->v = calloc(
			100 + ((*kat)->den_blocks * sm_size * sm_size) + (*kat)->_.nnz,
			sizeof(datatype_t));
	(*kat)->last_v = (*kat)->v;
	assert((*kat)->last_v == (*kat)->v);
	assert(*((*kat)->last_v) == *((*kat)->v));

	(*kat)->ci = calloc((*kat)->_.nnz, sizeof(int));
	(*kat)->rp = calloc((*kat)->blocks * ((*kat)->sm_size + 1), sizeof(int));

	_s_debug(KAT_DEBUG, "a\n");
	kat_determine_node(*kat, (*kat)->root);
	_s_debug(KAT_DEBUG, "b\n");
	kat_prepare(*kat, (*kat)->root);
	_s_debug(KAT_DEBUG, "c\n");

	kat_print_node((*kat)->root, 0);

	for (i = 0; i < mm_file->nnz; i++) {
		tmp_node = kat_get_node(*kat, mm_file->data[i].row,
				mm_file->data[i].col);

		switch (tmp_node->node_type) {
		case KAT_N_DEN:

			printf("Will write to dense sumbatrix. %d, %d, %p\n",
					((mm_file->data[i].row % (*kat)->sm_size) * (*kat)->sm_size)
							+ (mm_file->data[i].col % (*kat)->sm_size),
					((*kat)->den_blocks * sm_size + sm_size) + (*kat)->_.nnz
							- (*kat)->den_blocks_nnz,
					(void*) (&(tmp_node->node.sm.v)));

			_s_debugf(KAT_DEBUG, "tmp_node->node.sm.v = %p\n",
					tmp_node->node.sm.v);
			_s_debugf(KAT_DEBUG, "(*kat)->v = %p\n", (*kat)->v);
			_s_debugf(KAT_DEBUG, "(*kat)->last_v = %p\n", (*kat)->last_v);

//			assert(tmp_node->node.sm.v == (*kat)->v);

			/*
			 * tmp_node->node.sm.v is a pointer to the (*kat)->v array.
			 * It's a one dimensional array, but dense matrix is a two
			 * dimensional array. Thus we need to manually compute the
			 * position. (i.e. we are simulating v[a][b] via v[a*size + b])
			 *
			 * If dense matrices would be 100% dense, we could use v[n_nnz].
			 * But we allow some zero elements in dense matrix.
			 */
			tmp_node->node.sm.v[((mm_file->data[i].row % (*kat)->sm_size)
					* (*kat)->sm_size)
					+ (mm_file->data[i].col % (*kat)->sm_size)] =
					mm_file->data[i].value;

			_s_debug(KAT_DEBUG, "kk");

			break;
		case KAT_N_CSR:
			tmp_node->node.sm.v[tmp_node->node.sm.n_nnz] =
					mm_file->data[i].value;
			tmp_node->node.sm.s.csr.ci[tmp_node->node.sm.n_nnz] =
					(mm_file->data[i].col % (*kat)->sm_size);
			tmp_node->node.sm.s.csr.rp[(mm_file->data[i].row % (*kat)->sm_size)]++;
			break;
		default:
			_s_debugf(KAT_DEBUG, "Node type should be specific. node_type=%d\n",
					tmp_node->node_type);
			assert(0);
			break;
		}

		tmp_node->node.sm.n_nnz++;
	}

#if KAT_CSR
	kat_csr_rp(*kat, (*kat)->root);
#endif /* KAT_CSR */

	end_time = omp_get_wtime();
	mm_free(mm_file);
	return end_time - start_time;
}
