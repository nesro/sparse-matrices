/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
 * https://github.com/nesro/sparse-matrices
 */

#include "den_matrix.h"

#ifndef KAT_MATRIX_H_
#define KAT_MATRIX_H_

#define KAT_DEBUG 0

/*
 * The actual k of k-ary tree is N*N
 * If it's a constant, we may be sure that compiler will unroll it,
 * but we cannot modify it when testing.
 *
 * Is k is not set and it's not a constant, it's set by KAT_N_DEFAULT.
 *
 * TODO: kat node pointers will have to be calloc'd. so for now, it just
 * will be a constant
 */
#define KAT_N_IS_CONSTANT 1
#define KAT_N_DEFAULT 2

#if KAT_N_IS_CONSTANT
#define KAT_N KAT_N_DEFAULT
#define KAT_K (KAT_N*KAT_N)
#define GET_KAT_N(kat) \
	((void)0)
#else /* KAT_N_IS_CONSTANT */
#define GET_KAT_N(kat) \
	int KAT_N = kat->kat_n; \
	int KAT_K = kat->kat_k /* no ; */
#endif /* KAT_N_IS_CONSTANT */

/*
 * If a submatrix will have at least KAT_DENSE_TRESHOLD of nnz, it will be
 * treated like a dense matrix.
 */
#define KAT_DENSE_TRESHOLD (2)

/*
 * If true, leaves of k-ary tree will be (also) in CSR format.
 */
#define KAT_CSR 1

/*
 * kat_node is a node of a k-ary tree matrix
 * It contains a tagged union with either:
 * - an 2D array of pointers to next nodes
 * - a csr submatrix
 */
typedef enum kat_node_type {
	UNDEF, /**/
	INNER, /**/
	LEAF, /**/
	KAT_N_DEN, /**/
	KAT_N_CSR, /**/
} kat_node_type_t;
typedef struct kat_node {
	kat_node_type_t node_type;
	union {
		struct kat_node *knp[KAT_K][KAT_K]; /* kat node pointer */
		struct {
			int x;
			int y;
			datatype_t *v;

			/*
			 * When creating a k-ary tree, this variable is used to determine
			 * whether the node will be dense or sparse. If the node is dense
			 * this variable is unused.
			 */
			int nnz;

			/*
			 * When creating a k-ary tree, this holds number of nnz already
			 * loaded.
			 */
			int n_nnz;

			/*
			 * Sparse components of submatrix. If the matrix is dense, this
			 * union is unused.
			 *
			 * XXX: For now, only the CSR format is supported.
			 */
			union {
				struct {
					int *ci;
					int *rp;
				} csr;
			} s;
		} sm; /* submatrix */
	} node;
} kat_node_t;

typedef struct kat_matrix {
	vm_t _;

	kat_node_t *root;
	datatype_t *v; /* values */
	datatype_t *last_v;
	int v_length;

	int sm_size;
	int height;
	int blocks;

	int den_blocks;
	int den_blocks_nnz;

#if KAT_N_IS_CONSTANT == 0
	int kat_k;
	int kat_n;
#endif

#if KAT_CSR
	int *ci; /* column indices */
	int *last_ci;

	int *rp; /* row pointers */
	int *last_rp;

	int csr_blocks;
#endif /* KAT_CSR */

	/* private */
	datatype_t *rows_block;
} kat_matrix_t;

void kat_vm_init(kat_matrix_t **kat, va_list va);
void kat_init(kat_matrix_t **kat, int width, int height, int nnz, int sm_size);
void kat_free(kat_matrix_t *kat);
void kat_determine_node(kat_matrix_t *kat, kat_node_t *kat_node);
void kat_from_mm(kat_matrix_t **kat, const char *file, va_list va);
double kat_load_mm(kat_matrix_t **kat, const char *filename, int sm_size);
double kat_mul(const kat_matrix_t *a, const kat_matrix_t *b, den_matrix_t **c,
		char flag);
vm_t *kat_convert(kat_matrix_t *kat, vm_type_t type);

#endif /* KAT_MATRIX_H_ */
