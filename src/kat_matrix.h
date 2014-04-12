/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
 * https://github.com/nesro/sparse-matrices
 */

#ifndef KAT_MATRIX_H_
#define KAT_MATRIX_H_

#define KAT_DEBUG 1

/*
 * The actual k of k-ary tree is N*N;
 */
#define KAT_N 4
#define KAT_K (KAT_N*KAT_N)

/*
 * If a submatrix will have KAT_DENSE_TRESHOLD % of nnz, it will be treated
 * like a dense matrix.
 */
#define KAT_DENSE_TRESHOLD ((double)0.8)

/*
 * If true, leaves of k-ary tree will be (also) in csr format.
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

	int sm_size;
	int height;
	int blocks;

	int den_blocks;
	int den_blocks_nnz;

#if KAT_CSR
	int *ci; /* column indices */
	int *rp; /* row pointers */

	int csr_blocks;
	int last_index_ci;
	int last_index_rp;
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


#endif /* KAT_MATRIX_H_ */
