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
(convert_t) kat_convert, /**/
(mul_t) kat_mul, /**/
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

__attribute__((optimize("unroll-loops")))
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

//	int fy = 0;
//	int fx = 0;

	_s_debugf(KAT_DEBUG, "Searching for node for [y,x]=[%d,%d]\n", y, x);

	for (i = 0; i < kat->height; i++) {

	_s_debugf(KAT_DEBUG, "block_size=%d]\n",block_size);

	if (block_size < 1) block_size = 1;


		oy = (y / block_size);
		_s_debugf(KAT_DEBUG, "--- computing oy=(y/block_size) as %d=(%d/%d)\n",
				oy, y, block_size);

		ox = (x / block_size);
		y -= oy * block_size;
		x -= ox * block_size;
		by += oy * block_size;
		bx += ox * block_size;

		if (*tmp_node == NULL) {
			_s_debug(KAT_DEBUG, "creating a new inner node!\n");

			*tmp_node = calloc(1, sizeof(kat_node_t));
			assert(tmp_node != NULL);

			(*tmp_node)->node_type = INNER;
		}

		tmp_node = &((*tmp_node)->node.knp[oy][ox]);
		block_size /= KAT_N;

		_s_debugf(KAT_DEBUG,
				"next node is at oy=%d ox=%d bs=%d y=%d x=%d by=%d bx=%d\n", oy,
				ox, block_size, y, x, by, bx);
	}

	if (*tmp_node == NULL) {

		*tmp_node = calloc(1, sizeof(kat_node_t));
		assert(tmp_node != NULL);

		/*
		 * The submatrix will be sparse when created. Because elements are
		 * assigned after after all blocks are created, when nnz of submatrix
		 * is high, we switch it to dense format.
		 */
		(*tmp_node)->node_type = KAT_N_CSR;
		(*tmp_node)->node.sm.y = ((int) (by / kat->sm_size)) * kat->sm_size;
		(*tmp_node)->node.sm.x = ((int) (bx / kat->sm_size)) * kat->sm_size;
		kat->blocks++;

		_s_debugf(KAT_DEBUG, "Creating a leaf node at y=%d x=%d\n",
				(*tmp_node)->node.sm.y, (*tmp_node)->node.sm.x);
	}

	_s_debugf(KAT_DEBUG, "Node found. by=%d bx=%d y=%d x=%d\n",
			(*tmp_node)->node.sm.y, (*tmp_node)->node.sm.x,
			(*tmp_node)->node.sm.y, (*tmp_node)->node.sm.x);

	return *tmp_node;
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

			kat_node->node.sm.v = kat->last_v;
			kat->last_v = &(kat->last_v[kat->sm_size * kat->sm_size]);

			break;
		case KAT_N_CSR:

			kat_node->node.sm.v = kat->last_v;
			kat->last_v = &(kat->last_v[kat_node->node.sm.nnz]);

			kat_node->node.sm.s.csr.ci = kat->last_ci;
			kat->last_ci = &(kat->last_ci[kat_node->node.sm.nnz]);

			kat_node->node.sm.s.csr.rp = kat->last_rp;
			kat->last_rp = &(kat->last_rp[kat->sm_size + 1]);
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

	if (kat_node->node_type == INNER)
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

		if (tmp_node->node_type != KAT_N_DEN) {
			(*kat)->v_length++;
			tmp_node->node.sm.nnz++;
		}

		if (tmp_node->node_type
				!= KAT_N_DEN&& tmp_node->node.sm.nnz > KAT_DENSE_TRESHOLD) {

			_s_debugf(KAT_DEBUG, "tmp_node at y=%d x=%d is now DENSE\n",
					tmp_node->node.sm.y, tmp_node->node.sm.x);

			tmp_node->node_type = KAT_N_DEN;
			(*kat)->den_blocks++;
			(*kat)->v_length += sm_size * sm_size - tmp_node->node.sm.nnz;
		}

	}

	/*
	 * Allocating memory for matrix data.
	 */
	(*kat)->v = calloc((*kat)->v_length, sizeof(datatype_t));
	(*kat)->last_v = (*kat)->v;

	/*
	 * XXX: here we are computing size of csr arrays. if there will be
	 * another format, this part must be redesigned
	 */
	(*kat)->ci = calloc((*kat)->_.nnz, sizeof(int));
	(*kat)->rp = calloc(
			((*kat)->blocks - (*kat)->den_blocks) * ((*kat)->sm_size + 1),
			sizeof(int));

	kat_prepare(*kat, (*kat)->root);

//	kat_print_node((*kat)->root, 0);

	/*
	 * We'll put each element to the right place in its submatrix.
	 */
	for (i = 0; i < mm_file->nnz; i++) {
		tmp_node = kat_get_node(*kat, mm_file->data[i].row,
				mm_file->data[i].col);

		switch (tmp_node->node_type) {
		case KAT_N_DEN:

//			printf("Will write to dense sumbatrix. %d, %d, %p\n",
//					((mm_file->data[i].row % (*kat)->sm_size) * (*kat)->sm_size)
//							+ (mm_file->data[i].col % (*kat)->sm_size),
//					((*kat)->den_blocks * sm_size + sm_size) + (*kat)->_.nnz
//							- (*kat)->den_blocks_nnz,
//					(void*) (&(tmp_node->node.sm.v)));

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
			_s_debugf(KAT_DEBUG,
					"Adding %lf at y=%d x=%d i.e. index=%d node.y=%d mode.x=%d to DEN <<<<<<<<\n",
					mm_file->data[i].value, mm_file->data[i].row,
					mm_file->data[i].col,
					((mm_file->data[i].row % (*kat)->sm_size) * (*kat)->sm_size)
							+ (mm_file->data[i].col % (*kat)->sm_size),
					tmp_node->node.sm.y, tmp_node->node.sm.x);

			tmp_node->node.sm.v[((mm_file->data[i].row % (*kat)->sm_size)
					* (*kat)->sm_size)
					+ (mm_file->data[i].col % (*kat)->sm_size)] =
					mm_file->data[i].value;

			break;
		case KAT_N_CSR:

			_s_debugf(KAT_DEBUG, "adding %lf at y=%d x=%d to CSR\n",
					mm_file->data[i].value, mm_file->data[i].row,
					mm_file->data[i].col);

			tmp_node->node.sm.v[tmp_node->node.sm.n_nnz] =
					mm_file->data[i].value;
			tmp_node->node.sm.s.csr.ci[tmp_node->node.sm.n_nnz] =
					(mm_file->data[i].col % (*kat)->sm_size);
			tmp_node->node.sm.s.csr.rp[(mm_file->data[i].row % (*kat)->sm_size)]++;
			tmp_node->node.sm.n_nnz++;
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

/******************************************************************************/

static void kat_node_to_dense(kat_matrix_t *kat, kat_node_t *kat_node,
		den_matrix_t *den) {

	int i;
	int j;

	if (kat_node->node_type == KAT_N_CSR) {
		return;
	}

	if (kat_node->node_type == KAT_N_DEN) {
		for (i = 0; i < kat->sm_size; i++)
			for (j = 0; j < kat->sm_size; j++)
				den->v[i][j] += kat_node->node.sm.v[i * kat->sm_size + j];
		return;
	}

	if (kat_node->node_type == INNER)
		for (i = 0; i < KAT_N; i++)
			for (j = 0; j < KAT_N; j++)
				if (kat_node->node.knp[i][j])
					kat_node_to_dense(kat, kat_node->node.knp[i][j], den);
}

vm_t *kat_convert(kat_matrix_t *kat, vm_type_t type) {

	vm_t *vm = NULL;

	switch (type) {
	case DEN:
		vm_create(&vm, DEN, 1, kat->_.w, kat->_.h);
		kat_node_to_dense(kat, kat->root, (den_matrix_t *) vm);
		break;
	default:
		fdie("unknown format to convert %d\n", type);
		break;
	}

	/* XXX: */
	/* qdt->_.f.free(qdt); */

	return vm;
}

/******************************************************************************/

__attribute__((optimize("unroll-loops")))
static void kat_mul_den_den(const kat_matrix_t *ma, const kat_matrix_t *mb,
		const kat_node_t *a, const kat_node_t *b, den_matrix_t *c) {

	int i;
	int j;
	int k;
	int sms = ma->sm_size;

	_s_debugf(KAT_DEBUG,
			"mul_den_den a=%p a->node_type=%d a->node.sm.v=%p b->node.sm.v=%p\n",
			(void* )a, a->node_type, (void* )a->node.sm.v,
			(void* )b->node.sm.v);

	for (i = 0; i < sms; i++) {
		for (j = 0; j < sms; j++) {
			for (k = 0; k < sms; k++) {

				_s_debugf(KAT_DEBUG,
						"mul %lf at y=%d x=%d with %lf at y=%d x=%d to y=%d x=%d\n",
						a->node.sm.v[i * sms + k], i, k,
						b->node.sm.v[k * sms + j], k, j, i + a->node.sm.y,
						j + b->node.sm.x);

				c->v[i + a->node.sm.y][j + b->node.sm.x] += a->node.sm.v[i * sms
						+ k] * b->node.sm.v[k * sms + j];
			}
		}
	}
}

static void kat_mul_den_csr(const kat_matrix_t *ma, const kat_matrix_t *mb,
		const kat_node_t *a, const kat_node_t *b, den_matrix_t *c) {
	printf("cccc\n\n\n\n");

}
static void kat_mul_csr_den(const kat_matrix_t *ma, const kat_matrix_t *mb,
		const kat_node_t *a, const kat_node_t *b, den_matrix_t *c) {
	printf("bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb\n\n\n\n");

}
static void kat_mul_csr_csr(const kat_matrix_t *ma, const kat_matrix_t *mb,
		const kat_node_t *a, const kat_node_t *b, den_matrix_t *c) {
	printf("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n\n\n\n");
}

__attribute__((optimize("unroll-loops")))
static void kat_node_mul(const kat_matrix_t *ma, const kat_matrix_t *mb,
		const kat_node_t *a, const kat_node_t *b, den_matrix_t *c, int depth) {

	int i;
	int j;
	int k;

	for (i = 0; i < KAT_N; i++) {
		for (j = 0; j < KAT_N; j++) {
			for (k = 0; k < KAT_N; k++) {
				if (a->node.knp[i][k] && b->node.knp[k][j]) {

					if (a->node.knp[i][k]->node_type == KAT_N_DEN) {
						if (b->node.knp[k][j]->node_type == KAT_N_DEN) {
							kat_mul_den_den(ma, mb, a->node.knp[i][k],
									b->node.knp[k][j], c);
						} else if (b->node.knp[k][j]->node_type == KAT_N_CSR) {
							kat_mul_den_csr(ma, mb, a->node.knp[i][k],
									b->node.knp[k][j], c);
						}
						continue;
					}

					if (a->node.knp[i][k]->node_type == KAT_N_CSR) {
						if (b->node.knp[k][j]->node_type == KAT_N_DEN) {
							kat_mul_csr_den(ma, mb, a->node.knp[i][k],
									b->node.knp[k][j], c);
						} else if (b->node.knp[k][j]->node_type == KAT_N_CSR) {
							kat_mul_csr_csr(ma, mb, a->node.knp[i][k],
									b->node.knp[k][j], c);
						}
						continue;
					}

					kat_node_mul(ma, mb, a->node.knp[i][k], b->node.knp[k][j],
							c, depth);
				}
			}
		}
	}
}

double kat_mul(const kat_matrix_t *a, const kat_matrix_t *b, den_matrix_t **c,
		char flag /* unused */) {

	double start_time;
	double end_time;

	assert(a->_.w == b->_.w);
	assert(a->_.h == b->_.h);

	if (*c == NULL)
		vm_create((vm_t **) c, DEN, 1, a->_.w, a->_.w);

	start_time = omp_get_wtime();
	kat_node_mul(a, b, a->root, b->root, *c, 0);
	end_time = omp_get_wtime();

	return end_time - start_time;
}

