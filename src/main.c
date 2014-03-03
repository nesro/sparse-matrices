/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
 * https://github.com/nesro/sparse-matrices
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <getopt.h>

#include "utils.h"
#include "csr_matrix.h"
#include "coo_matrix.h"
#include "den_matrix.h"
#include "qdt_matrix.h"

static void usage() {

	printf("Usage: ./sparse-matrix-multiplication -f <format> -a <matrix_a> -b"
		" <matrix_b> -o [none|file] -l leaf_size\n"
		"Format could be: quadtree|csr|dense\n"
		"Matrices a and b has to be in MatrixMarket file format.\n");
}

//static matrix_type_t parse_format(const char *format_string) {
//
//	if (strcmp(format_string, "quadtree") == 0) {
//		return QUADTREE;
//	} else if (strcmp(format_string, "dense") == 0) {
//		return DENSE;
//	} else if (strcmp(format_string, "csr") == 0) {
//		return CSR;
//	} else {
//		fprintf(stderr, "Unknown matrix format.\n");
//		exit(EXIT_FAILURE);
//	}
//}
//
//int load_leaf_size(const char *string) {
//
//	char *tmp = NULL;
//	int leaf_size;
//
//	leaf_size = strtol(string, &tmp, 10);
//
//	if (tmp == NULL) {
//		fprintf(stderr, "ERROR: Cannot load leaf size. Invalid input.\n");
//		exit(EXIT_FAILURE);
//	}
//
//	if (leaf_size < 2) {
//		fprintf(stderr, "ERROR: Minimum of leaf size is 2.\n");
//		exit(EXIT_FAILURE);
//	}
//
//	if (!is_power_of_two(leaf_size)) {
//		fprintf(stderr, "ERROR: Leaf size is not power of two.\n");
//		exit(EXIT_FAILURE);
//	}
//
//	return leaf_size;
//}

int main(int argc, char *argv[]) {

	vm_t *a = NULL;
	vm_t *b = NULL;
	vm_t *c = NULL;
	double time;


//	vm_create(&vm, DEN, 10, 5, 3, VA_END);
	vm_load_mm(&a, DEN, "4x4_4nz_01.mtx");
	vm_load_mm(&b, DEN, "4x4_4nz_01.mtx");

	time = a->f.mul(a, b, &c, NAIVE | UNROLLED);
	printf("time=%lf\n", time);

	a->f.print(a);
	b->f.print(b);
	b->f.print(c);

	a->f.free(a);
	b->f.free(b);
	c->f.free(c);

//	/*usage();*/
//
//	vm_t *dm = NULL;
//
////	den_matrix_init((den_matrix_t **) &dm, (den_matrix_init_t ) { 10, 10, 1 });
//
//	vm_init(&dm, DENSE, (vm_init_t));
//
//	printf("yay %x\n", dm);
//
//	dm->vmt.free(dm);
//	printf("yey\n");

//
//	int c;
//	int index;
//	matrix_type_t format;
//	char *matrix_a = NULL;
//	char *matrix_b = NULL;
//	int leaf_size = -1;
//	time_record_t tr;
//
//	if (argc < 2) {
//		usage();
//		return EXIT_SUCCESS;
//	}
//
//	opterr = 0;
//	while ((c = getopt(argc, argv, "a:b:f:l:o:")) != -1) {
//		switch (c) {
//		case 'a':
//			matrix_a = optarg;
//			break;
//		case 'b':
//			matrix_b = optarg;
//			break;
//		case 'f':
//			format = parse_format(optarg);
//			break;
//		case 'l':
//			leaf_size = load_leaf_size(optarg);
//			break;
//		case 'o':
//			/* TODO */
//			break;
//		case '?':
//			if (optopt == 'a')
//				fprintf(stderr,
//					"ERROR: Filename is missing for the option -a.\n");
//			else if (optopt == 'b')
//				fprintf(stderr,
//					"ERROR:Filename is missing for the option -b.\n");
//			else if (optopt == 'f')
//				fprintf(stderr, "ERROR: Format is missing for the option -f. "
//					"A format could be: quadtree or csr.\n");
//			else if (isprint(optopt))
//				fprintf(stderr, "ERROR: Unknown option \"-%c\".\n", optopt);
//			else
//				fprintf(stderr, "ERROR: Unknown option character \"\\x%x\".\n",
//					optopt);
//			return EXIT_FAILURE;
//		default:
//			abort();
//		}
//	}
//
//	for (index = optind; index < argc; index++)
//		printf("WARNING: Needless argument: \"%s\"\n", argv[index]);
//
//	if (matrix_a == NULL) {
//		fprintf(stderr, "ERROR: Matrix a is not set.\n");
//		usage();
//		return EXIT_FAILURE;
//	}
//
//	if (matrix_b == NULL)
//		matrix_b = matrix_a;
//
//	if (format == QUADTREE) {
//		if (leaf_size == -1) {
//			fprintf(stderr,
//				"ERROR: Leaf size is not set. Set it with an -l option.\n");
//			return EXIT_FAILURE;
//		}
//	}
//
//	switch (format) {
//	case QUADTREE:
//		tr = qt_matrix_mm_mul(matrix_a, matrix_b, leaf_size);
//		break;
//	case CSR:
//		tr = csr_matrix_mm_mul(matrix_a, matrix_b);
//		break;
//	default:
//		abort();
//	}
//
//	printf("load_a %lf\nload_b %lf\nmultiplication %lf\n", tr.load_a, tr.load_b,
//		tr.multiplication);

	return EXIT_SUCCESS;
}
