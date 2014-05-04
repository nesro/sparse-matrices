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

	/*
	 * XXX: be strict here. It's for testing purposes and we want to avoid
	 * edge cases.
	 */
	assert(width == height);
	assert(is_power_of_two(width));
	assert(is_power_of_two(sm_size));
	assert(is_power_of_two(KAT_N));

	*kat = calloc(1, sizeof(kat_matrix_t));
	assert(*kat != NULL);
	(*kat)->_.object_size += sizeof(kat_matrix_t);

	(*kat)->_.type = KAT;
	(*kat)->_.f = kat_vmt;
	(*kat)->_.w = width;
	(*kat)->_.h = height;

	(*kat)->sm_size = sm_size;
	(*kat)->_.nnz = nnz;

	(*kat)->root = NULL;

	/*
	 * This is not needed.
	 */
	(*kat)->height =
			ceil(
					(log((width * height) / (sm_size * sm_size))
							/ (double) log(KAT_K)));

	printf("_kat init h=%d w=%d height=%d\n", width, height, (*kat)->height);

	// not really
//	assert(width < sm_size * KAT_N);

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

#if KAT_DEBUG

static int is_uniq_pos(kat_matrix_t *kat, kat_node_t *kat_node, int y, int x) {

	int i;
	int j;
	int result = 0;

	if (kat_node->node_type == KAT_N_CSR || kat_node->node_type == KAT_N_DEN) {
		if (kat_node->node.sm.y == y && kat_node->node.sm.x == x)
		return 1;
		else
		return 0;
	}

	if (kat_node->node_type == INNER) {
		for (i = 0; i < KAT_N; i++)
		for (j = 0; j < KAT_N; j++)
		if (kat_node->node.knp[i][j])
		result += is_uniq_pos(kat, kat_node->node.knp[i][j], y, x);

		return result;
	}

	assert(0);
	return 1;
}
#endif
/*
 * This function return a pointer to the right block. Automatically build
 * the path to the block.
 */
static kat_node_t *kat_get_node(kat_matrix_t *kat, int y, int x) {

#if KAT_DEBUG
	int tmp_height = 0;
#endif

	/*
	 * Pointer to a node. We'll use it to traverse through the tree.
	 */
	kat_node_t **tmp_node = &kat->root;

	/*
	 * This is a temporary block to determine in which part we have to traverse.
	 * We're starting at the full width of the matrix and dividing it by the
	 * KAT_N number. When submatrix_size * KAT_N exceeds this block_size, we
	 * know that we're in the last node and may create a submatrix.
	 */
	int block_start_y = 0;
	int block_start_x = 0;
	int block_size = kat->_.w;

	/*
	 * We can compute the sumbatrix's y and x coords even before traversing
	 * the tree.
	 */
	int sm_y = ((int) (y / kat->sm_size)) * kat->sm_size;
	int sm_x = ((int) (x / kat->sm_size)) * kat->sm_size;

	/*
	 * Position of the next node when traversing the tree.
	 */
	int node_y;
	int node_x;

	_s_debugf(KAT_DEBUG, "item: y=%d, x=%d in block y=%d, x=%d\n", y, x, sm_y,
			sm_x);

	/*
	 * Now traverse to the last inner block.
	 */
	while (kat->sm_size * KAT_N < block_size) {

		assert(y >= block_start_y);
		assert(block_start_y + block_size >= y);
		assert(x >= block_start_x);
		assert(block_start_x + block_size >= x);

		node_y = (y - block_start_y) / (block_size / KAT_N);
		node_x = (x - block_start_x) / (block_size / KAT_N);

		assert(0 <= node_y);
		assert(node_y <= KAT_N);
		assert(0 <= node_x);
		assert(node_x <= KAT_N);

		if (*tmp_node == NULL) {
			_s_debugf(KAT_DEBUG,
					"creating an inner node knp[%d][%d] when block_start_y=%d, block_start_x=%d, block_size=%d\n",
					node_y, node_x, block_start_y, block_start_x, block_size);

			*tmp_node = calloc(1, sizeof(kat_node_t));
			kat->_.object_size += sizeof(kat_node_t);
			assert(*tmp_node != NULL);

			(*tmp_node)->node_type = INNER;
			kat->nodes_inner++;
		}

		tmp_node = &((*tmp_node)->node.knp[node_y][node_x]);

		block_start_y += node_y * (block_size / KAT_N);
		block_start_x += node_x * (block_size / KAT_N);
		block_size /= KAT_N;

#if KAT_DEBUG
		tmp_height++;

		if (tmp_height >= kat->height) {
			printf("tree to hight=%d max=%d", tmp_height, kat->height);
			exit(1);
		}
#endif
	}

	/*
	 * Now we have to compute real offset with submatrix size.
	 * The last node before submatrix may be not full.
	 */
	assert(y >= block_start_y);
	assert(block_start_y + block_size >= y);
	assert(x >= block_start_x);
	assert(block_start_x + block_size >= x);

	node_y = (y - block_start_y) / kat->sm_size;
	node_x = (x - block_start_x) / kat->sm_size;
	block_start_y += node_y * kat->sm_size;
	block_start_x += node_x * kat->sm_size;

	if (*tmp_node == NULL) {
		_s_debugf(KAT_DEBUG,
				"creating a LAST inner node knp[%d][%d] when block_start_y=%d, block_start_x=%d, block_size=%d\n",
				node_y, node_x, block_start_y, block_start_x, block_size);

		*tmp_node = calloc(1, sizeof(kat_node_t));
		kat->_.object_size += sizeof(kat_node_t);
		assert(*tmp_node != NULL);

		(*tmp_node)->node_type = INNER;
		kat->nodes_inner++;
	}
	tmp_node = &((*tmp_node)->node.knp[node_y][node_x]);

	/*
	 * We should traverse to the same location as we computed before.
	 */
	assert(block_start_y == sm_y);
	assert(block_start_x == sm_x);

	if (*tmp_node == NULL) {

#if KAT_DEBUG
		if (is_uniq_pos(kat, kat->root, sm_y, sm_x) == 1) {
			printf("there is block y=%d, x=%d", sm_y, sm_x);
			exit(1);
		}
#endif

		*tmp_node = calloc(1, sizeof(kat_node_t));
		kat->_.object_size += sizeof(kat_node_t);
		assert(*tmp_node != NULL);

		/*
		 * The submatrix will be sparse when created. Because elements are
		 * assigned after after all blocks are created, when nnz of submatrix
		 * is high, we switch it to dense format.
		 */
		(*tmp_node)->node_type = KAT_N_CSR;
		(*tmp_node)->node.sm.y = sm_y;
		(*tmp_node)->node.sm.x = sm_x;
		kat->blocks++;
		kat->nodes_csr++;

		_s_debugf(KAT_DEBUG, "Creating a leaf node at y=%d x=%d\n",
				(*tmp_node)->node.sm.y, (*tmp_node)->node.sm.x);
	}

	_s_debug(KAT_DEBUG, "returning the node\n");

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

			kat_node->node.sm.s.den.v = malloc(
					kat->sm_size * sizeof(datatype_t *));

			for (i = 0; i < kat->sm_size; i++)
				kat_node->node.sm.s.den.v[i] = &(kat->last_v[i * kat->sm_size]);

			kat->last_v = &(kat->last_v[kat->sm_size * kat->sm_size]);

			break;
		case KAT_N_CSR:

			kat_node->node.sm.s.csr.v = kat->last_v;
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
		_s_debugf(KAT_DEBUG, "computing block at y=%d, x=%d\n",
				kat_node->node.sm.y, kat_node->node.sm.x);

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
				if (kat_node->node.knp[i][j])
					kat_csr_rp(kat, kat_node->node.knp[i][j]);
}

/******************************************************************************/

void kat_print_node(kat_node_t *kat_node, int depth) {

	int i;
	int j;

	if (kat_node->node_type != INNER) {
		printf("%*s END type=%d \n", depth, "", kat_node->node_type);
		return;
	}

	for (i = 0; i < KAT_N; i++) {
		for (j = 0; j < KAT_N; j++) {
			printf("%*s[%d][%d] = \n", depth, "", i, j);
			if (kat_node->node.knp[i][j] != NULL) {
				switch (kat_node->node.knp[i][j]->node_type) {
				case KAT_N_DEN:
					printf("  den mul y=%d, x=%d\n",
							kat_node->node.knp[i][j]->node.sm.y,
							kat_node->node.knp[i][j]->node.sm.x);
					break;
				case KAT_N_CSR:
					printf("  csr mul y=%d, x=%d\n",
							kat_node->node.knp[i][j]->node.sm.y,
							kat_node->node.knp[i][j]->node.sm.x);
					break;
				case INNER:
					printf("  inner\n");
					kat_print_node(kat_node->node.knp[i][j], depth + 1);
					break;
				default:
					assert(0);
					break;
				}
			} else {
				printf("NULL\n");
			}
		}
	}
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
	kat_node_t *tn; /* temporary node */
	int i;

	mm_file = mm_load(filename, 1);

	/*
	 * Start timer after mm_file is loaded. We don't want to IO operations
	 * in our timer.
	 */
	start_time = omp_get_wtime();

	kat_init(kat, mm_file->width, mm_file->height, mm_file->nnz, sm_size);

	for (i = 0; i < mm_file->nnz; i++) {

		_s_debugf(KAT_DEBUG, "i=%d v="DPF"\n", i, mm_file->data[i].value);

		tn = kat_get_node(*kat, mm_file->data[i].row, mm_file->data[i].col);

		/*
		 * If the matrix is dense, there is no reason to count its elements.
		 */
		if (tn->node_type != KAT_N_DEN) {
			(*kat)->v_length++;
			tn->node.sm.nnz++;
		}

		/*
		 * Every block is CSR by default. If there is too many elemtents
		 * inside, we'll treat it like a dense matrix.
		 */
		if (tn->node_type != KAT_N_DEN && tn->node.sm.nnz > KAT_DENSE_TRESHOLD) {

			_s_debugf(KAT_DEBUG, "tmp_node at y=%d x=%d is now DENSE\n",
					tn->node.sm.y, tn->node.sm.x);

			tn->node_type = KAT_N_DEN;
			(*kat)->nodes_csr--;
			(*kat)->nodes_den++;
			(*kat)->v_length += sm_size * sm_size - tn->node.sm.nnz;
		}
	}

	/*
	 * Allocating memory for the matrix data.
	 */
	(*kat)->v = calloc((*kat)->v_length, sizeof(datatype_t));
	_s_debugf(1, "kat->v size = %ld\n", (*kat)->v_length);
	assert((*kat)->v !=NULL);

	(*kat)->_.object_size += (*kat)->v_length * sizeof(datatype_t);
	(*kat)->last_v = (*kat)->v;

	/*
	 * XXX: here we are computing size of CSR arrays. if there will be
	 * another format, this part must be redesigned
	 */
	(*kat)->ci = calloc((*kat)->_.nnz, sizeof(int));
	assert((*kat)->ci!=NULL);
	(*kat)->_.object_size += (*kat)->_.nnz, sizeof(int);
	(*kat)->last_ci = (*kat)->ci;
	(*kat)->rp = calloc(
			((*kat)->blocks - (*kat)->den_blocks) * ((*kat)->sm_size + 1),
			sizeof(int));
	assert((*kat)->rp!=NULL);
	(*kat)->_.object_size += ((*kat)->blocks - (*kat)->den_blocks)
			* ((*kat)->sm_size + 1) * sizeof(int);
	(*kat)->last_rp = (*kat)->rp;

	kat_prepare(*kat, (*kat)->root);

	/*
	 * We'll put each element to the right place in its submatrix.
	 */
	for (i = 0; i < mm_file->nnz; i++) {

		tn = kat_get_node(*kat, mm_file->data[i].row, mm_file->data[i].col);

		switch (tn->node_type) {
		case KAT_N_DEN:
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
					"Adding "DPF" at y=%d x=%d i.e. index=%d node.y=%d mode.x=%d to DEN <<<<<<<<\n",
					mm_file->data[i].value, mm_file->data[i].row,
					mm_file->data[i].col,
					((mm_file->data[i].row % (*kat)->sm_size) * (*kat)->sm_size)
							+ (mm_file->data[i].col % (*kat)->sm_size),
					tn->node.sm.y, tn->node.sm.x);

			/* newden */
			tn->node.sm.s.den.v/**/
			[(mm_file->data[i].row % (*kat)->sm_size)]/**/
			[(mm_file->data[i].col % (*kat)->sm_size)]/**/
			= mm_file->data[i].value;

			break;
		case KAT_N_CSR:
			_s_debugf(KAT_DEBUG, "adding "DPF" at y=%d x=%d to CSR. n_nnz=%d\n",
					mm_file->data[i].value, mm_file->data[i].row,
					mm_file->data[i].col, tn->node.sm.n_nnz);

			tn->node.sm.s.csr.v[tn->node.sm.n_nnz] = mm_file->data[i].value;
			tn->node.sm.s.csr.ci[tn->node.sm.n_nnz] = (mm_file->data[i].col
					% (*kat)->sm_size);
			tn->node.sm.s.csr.rp[(mm_file->data[i].row % (*kat)->sm_size)]++;
			break;
		default:
			_s_debugf(KAT_DEBUG, "Node type should be specific. node_type=%d\n",
					tn->node_type);
			assert(0);
			break;
		}

		tn->node.sm.n_nnz++;
	}

#if KAT_DEBUG
	kat_print_node((*kat)->root, 0);
#endif

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
		_s_debugf(KAT_DEBUG, "converting CSR node at y=%d, x=%d\n",
				kat_node->node.sm.y, kat_node->node.sm.x);

		for (i = 0; i <= kat->sm_size; i++) {
			_s_debugf(KAT_DEBUG, "csr i=%d rp=%d\n", i,
					kat_node->node.sm.s.csr.rp[i]);
		}
		for (i = 0; i < kat_node->node.sm.n_nnz; i++) {
			_s_debugf(KAT_DEBUG, "csr i=%d ci=%d v="DPF"\n", i,
					kat_node->node.sm.s.csr.ci[i],
					kat_node->node.sm.s.csr.v[i]);
		}

		for (i = 0; i < kat->sm_size; i++) {
			for (j = kat_node->node.sm.s.csr.rp[i];
					j < kat_node->node.sm.s.csr.rp[i + 1]; j++) {
				_s_debugf(KAT_DEBUG, "i=%d, j=%d\n", i, j);
				den->v /**/
				[kat_node->node.sm.y + i] /**/
				[kat_node->node.sm.x + kat_node->node.sm.s.csr.ci[j]] +=
						kat_node->node.sm.s.csr.v[j];
			}
		}
		return;
	}

	if (kat_node->node_type == KAT_N_DEN) {
		for (i = 0; i < kat->sm_size; i++)
			for (j = 0; j < kat->sm_size; j++)
				den->v[kat_node->node.sm.y + i][kat_node->node.sm.x + j] +=
						kat_node->node.sm.s.den.v[i][j];
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

	int l;
	int m;
	int k;
	int sms = ma->sm_size;

	_s_debugf(KAT_DEBUG,
			"mul_den_den a=%p a->node_type=%d a->node.sm.v=%p b->node.sm.v=%p\n",
			(void* )a, a->node_type, (void* )a->node.sm.s.den.v,
			(void* )b->node.sm.s.den.v);

	for (l = 0; l < sms; l++) {
		for (m = 0; m < sms; m++) {
			for (k = 0; k < sms; k++) {

				_s_debugf(KAT_DEBUG,
						"mul "DPF" at y=%d x=%d with "DPF" at y=%d x=%d to y=%d x=%d\n",
						a->node.sm.s.den.v[l][k], l, k,
						b->node.sm.s.den.v[k][m], k, m, l + a->node.sm.y,
						m + b->node.sm.x);

				c->v[l + a->node.sm.y][m + b->node.sm.x] +=
						a->node.sm.s.den.v[l][k] * b->node.sm.s.den.v[k][m];
			}
		}
	}
}

static void kat_mul_den_csr(const kat_matrix_t *ma, const kat_matrix_t *mb,
		const kat_node_t *a, const kat_node_t *b, den_matrix_t *c) {

	int i;
	int j;
	int k;
	int sms = ma->sm_size;

	for (i = 0; i < sms; i++) {
		for (j = 0; j < sms; j++) {
			for (k = b->node.sm.s.csr.rp[j]; k < b->node.sm.s.csr.rp[j + 1];
					k++) {

//					_s_debugf(KAT_DEBUG,
//							"mul %lf at y=%d x=%d with %lf at y=%d x=%d to y=%d x=%d\n",
//							a->node.sm.v[i * sms + k], i, k,
//							b->node.sm.v[k * sms + j], k, j, i + a->node.sm.y,
//							j + b->node.sm.x);
//				_s_debugf(KAT_DEBUG, "j=%d", j);


				c->v[i + a->node.sm.y][b->node.sm.s.csr.ci[k] + b->node.sm.x] +=
						a->node.sm.s.den.v[i][j] * b->node.sm.s.csr.v[k];
			}
		}
	}

}

static inline void kat_mul_csr_den(const kat_matrix_t *ma,
		const kat_matrix_t *mb, const kat_node_t *a, const kat_node_t *b,
		den_matrix_t *c) {

	int r;
	int ac;
	int i;

//	_s_debugf(KAT_DEBUG,
//			"mul_csr_den a=%p a->node_type=%d a->node.sm.v=%p b->node.sm.v=%p\n",
//			(void* )a, a->node_type, (void* )a->node.sm.v,
//			(void* )b->node.sm.v);

	for (r = 0; r < ma->sm_size; r++)
		for (ac = a->node.sm.s.csr.rp[r]; ac < a->node.sm.s.csr.rp[r + 1]; ac++)
			for (i = 0; i < ma->sm_size; i++) {

//				_s_debugf(KAT_DEBUG,
//						">> mul "DPF" at y=%d x=%d with "DPF" at y=%d x=%d to y=%d x=%d\n",
//						a->node.sm.v[ac], 0, 0,
//						b->node.sm.v[ma->sm_size * a->node.sm.s.csr.ci[ac] + i],
//						0, 0, a->node.sm.y + r, b->node.sm.x + i);

				c->v[a->node.sm.y + r][b->node.sm.x + i] +=
						a->node.sm.s.csr.v[ac]
								* b->node.sm.s.den.v[a->node.sm.s.csr.ci[ac]][i];
			}
}

static inline void kat_mul_csr_csr(const kat_matrix_t *ma,
		const kat_matrix_t *mb, const kat_node_t *a, const kat_node_t *b,
		den_matrix_t *c) {

	int r;
	int ac;
	int bc;

	for (r = 0; r < ma->sm_size; r++)
		for (ac = a->node.sm.s.csr.rp[r]; ac < a->node.sm.s.csr.rp[r + 1]; ac++)
			for (bc = b->node.sm.s.csr.rp[a->node.sm.s.csr.ci[ac]];
					bc < b->node.sm.s.csr.rp[a->node.sm.s.csr.ci[ac] + 1]; bc++)
				c->v[a->node.sm.y + r][b->node.sm.x + b->node.sm.s.csr.ci[bc]] +=
						a->node.sm.s.csr.v[ac] * b->node.sm.s.csr.v[bc];

}

/*
 * Because of having multiple types of leaves, we need to determine which
 * function we'll call.
 */
__attribute__((optimize("unroll-loops")))
static void kat_node_mul(const kat_matrix_t *ma, const kat_matrix_t *mb,
		const kat_node_t *a, const kat_node_t *b, den_matrix_t *c) {

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
							c);
				}
			}
		}
	}
}

/******************************************************************************/
/*
 * K-ary tree matrix multiplied by a vector.
 */

__attribute__((optimize("unroll-loops")))
static inline void kat_mul_den_vec(const kat_matrix_t *ma, const vec_t *b,
		const kat_node_t *a, vec_t *c) {

	int i;
	int j;

	for (i = 0; i < ma->sm_size; i++)
		for (j = 0; j < ma->sm_size; j++)
			c->v[i + a->node.sm.y] += a->node.sm.s.den.v[i][j]
					* b->v[a->node.sm.x + j];
}

__attribute__((optimize("unroll-loops")))
static inline void kat_mul_csr_vec(const kat_matrix_t *ma, const vec_t *b,
		const kat_node_t *a, vec_t *c) {

	int r;
	int ac;

	_s_debugf(KAT_DEBUG, "c->v = %p\n", (void * )c->v);

	for (r = 0; r < ma->sm_size; r++)
		for (ac = a->node.sm.s.csr.rp[r]; ac < a->node.sm.s.csr.rp[r + 1]; ac++)
			c->v[a->node.sm.y + r] += /**/
			a->node.sm.s.csr.v[ac] * b->v[/**/
			a->node.sm.x + /**/
			a->node.sm.s.csr.ci[ac]];
}
#if KAT_DEBUG
long int tmp_knv = 0;
long int tmp_knvden = 0;
long int tmp_knvcsr = 0;
long int tmp_knvinn = 0;
#endif

/*
 * Multiply k-ary tree matrix nodes with a vector.
 */
__attribute__((optimize("unroll-loops")))
static void mul_katnode_vec(const kat_matrix_t *ma, const vec_t *mb,
		const kat_node_t *a, vec_t *c) {

	int i;
	int j;

#if KAT_DEBUG
	if (tmp_knv++ % 1 == 0) {
		printf(
				"deep=%d, dnv=%ld, tmp_knvden=%ld, tmp_knvcsr=%ld, tmp_knvinn=%ld\n",
				deep, tmp_knv, tmp_knvden, tmp_knvcsr, tmp_knvinn);
		fflush(stdout);
	}
#endif

	for (i = 0; i < KAT_N; i++) {
		for (j = 0; j < KAT_N; j++) {
			if (a->node.knp[i][j] != NULL) {
				switch (a->node.knp[i][j]->node_type) {
				case KAT_N_DEN:
#if KAT_DEBUG
					tmp_knvden++;
					printf("  den mul y=%d, x=%d\n",
							a->node.knp[i][j]->node.sm.y,
							a->node.knp[i][j]->node.sm.x);
#endif
					kat_mul_den_vec(ma, mb, a->node.knp[i][j], c);
					break;
				case KAT_N_CSR:
#if KAT_DEBUG
					tmp_knvcsr++;
					printf("  csr mul y=%d, x=%d\n",
							a->node.knp[i][j]->node.sm.y,
							a->node.knp[i][j]->node.sm.x);
#endif
					kat_mul_csr_vec(ma, mb, a->node.knp[i][j], c);
					break;
				case INNER:
#if KAT_DEBUG
					tmp_knvinn++;
#endif
					mul_katnode_vec(ma, mb, a->node.knp[i][j], c);
					break;
				default:
					assert(0);
					break;
				}
			}
		}
	}
}
/******************************************************************************/

double kat_mul(const kat_matrix_t *a, const vm_t *b, vm_t **c,
		char flag /* unused */) {

	double start_time;
	double end_time;

	switch (b->type) {
	case VEC:
		assert(a->_.w == b->h);
		vec_init((vec_t **) c, a->_.w);
		start_time = omp_get_wtime();
		mul_katnode_vec(a, (vec_t *) b, a->root, (vec_t *) *c);
		end_time = omp_get_wtime();
		break;
	case KAT:
		assert(a->_.w == b->w);
		assert(a->_.h == b->h);
		den_matrix_init((den_matrix_t **) c, a->_.w, b->w, 1);
		start_time = omp_get_wtime();
		kat_node_mul(a, (kat_matrix_t *) b, a->root, ((kat_matrix_t*) b)->root,
				(den_matrix_t *) *c);
		end_time = omp_get_wtime();
		break;
	default:
		fprintf(stderr, "Unknown matrix type: %d\n", b->type);
		exit(1);
		return 0.;
	}

	return end_time - start_time;
}

/******************************************************************************/

