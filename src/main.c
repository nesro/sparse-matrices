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

	printf("Usage: ./main "
			"-f <format> "
			"-a <matrix_a> "
			"-b <matrix_b> "
			"-o [none|file] "
			"-l leaf_size\n"
			"Formats: quadtree|csr|dense\n"
			"Matrices a and b has to be in the MatrixMarket file format. "
			"See: math.nist.gov/MatrixMarket/‎\n");
}

static vm_type_t parse_format(const char *format_string) {

	if (strcmp(format_string, "quadtree") == 0) {
		return QDT;
	} else if (strcmp(format_string, "dense") == 0) {
		return DEN;
	} else if (strcmp(format_string, "csr") == 0) {
		return CSR;
	} else {
		fprintf(stderr, "Unknown matrix format.\n");
		exit(EXIT_FAILURE);
	}
}

int load_leaf_size(const char *string) {

	char *tmp = NULL;
	int leaf_size;

	leaf_size = strtol(string, &tmp, 10);

	if (tmp == NULL) {
		fprintf(stderr, "ERROR: Cannot load leaf size. Invalid input.\n");
		exit(EXIT_FAILURE);
	}

	if (leaf_size < 2) {
		fprintf(stderr, "ERROR: Minimum of leaf size is 2.\n");
		exit(EXIT_FAILURE);
	}

	if (!is_power_of_two(leaf_size)) {
		fprintf(stderr, "ERROR: Leaf size is not power of two.\n");
		exit(EXIT_FAILURE);
	}

	return leaf_size;
}

void quick_test() {

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
}

int main(int argc, char *argv[]) {

	int c;
	int index;
	vm_type_t format = UNKNOWN;
	char *matrix_a = NULL;
	char *matrix_b = NULL;
	int leaf_size = -1;
	time_record_t tr;

	char *output = NULL; /* -o flag */

	vm_t *vm_a = NULL;
	vm_t *vm_b = NULL;
	vm_t *vm_c = NULL;
	int print_time = 0;
	double time;

	if (argc < 2) {
		fprintf(stderr,
				"Error: Not enough parameters. Run \"%s -h\" for more informations.\n",
				argv[0]);
		return EXIT_SUCCESS;
	}

	opterr = 0;
	while ((c = getopt(argc, argv, "a:b:f:hl:o:t")) != -1) {
		switch (c) {
		case 'a':
			matrix_a = optarg;
			break;
		case 'b':
			matrix_b = optarg;
			break;
		case 'f':
			format = parse_format(optarg);
			break;
		case 'h':
			usage();
			return EXIT_SUCCESS;
		case 'l':
			leaf_size = load_leaf_size(optarg);
			break;
		case 'o':
			output = optarg;
			break;
		case 'q':
			quick_test();
			return EXIT_SUCCESS;
		case 't':
			print_time = 1;
			break;
		case '?':
			if (optopt == 'a')
				fprintf(stderr,
						"ERROR: Filename is missing for the option -a.\n");
			else if (optopt == 'b')
				fprintf(stderr,
						"ERROR:Filename is missing for the option -b.\n");
			else if (optopt == 'f')
				fprintf(stderr, "ERROR: Format is missing for the option -f. "
						"A format could be: quadtree or csr.\n");
			else if (isprint(optopt))
				fprintf(stderr, "ERROR: Unknown option \"-%c\".\n", optopt);
			else
				fprintf(stderr, "ERROR: Unknown option character \"\\x%x\".\n",
						optopt);
			return EXIT_FAILURE;
		default:
			abort();
		}
	}

	for (index = optind; index < argc; index++)
		printf("WARNING: Needless argument: \"%s\"\n", argv[index]);

	if (matrix_a == NULL && matrix_b == NULL) {
		fprintf(stderr, "ERROR: Matrix a nor b are not set.\n");
		usage();
		return EXIT_FAILURE;
	}

	if (matrix_a == NULL && matrix_b != NULL)
		matrix_a = matrix_b;

	if (matrix_b == NULL && matrix_a != NULL)
		matrix_b = matrix_a;

	if (format == QDT) {
		if (leaf_size == -1) {
			fprintf(stderr,
					"ERROR: Leaf size is not set or has an invalid value."
							"Set it with an -l option to a power of two.\n");
			return EXIT_FAILURE;
		}
	}

	switch (format) {
	case DEN:
		break;
	case QDT:
		//tr = qt_matrix_mm_mul(matrix_a, matrix_b, leaf_size);
		break;
	case CSR:
		//tr = csr_matrix_mm_mul(matrix_a, matrix_b);
		break;
	case UNKNOWN:
		fprintf(stderr, "Format MUST be set.\n");
		return EXIT_FAILURE;
	default:
		abort();
	}

	vm_load_mm(&vm_a, format, matrix_a);
	vm_load_mm(&vm_b, format, matrix_b);

	time = vm_a->f.mul(vm_a, vm_b, &vm_c, NAIVE | UNROLLED);

	if (print_time)
		printf("%lf\n", time);

	if (output != NULL)
		vm_c->f.print(vm_c);

	printf("load_a %lf\nload_b %lf\nmultiplication %lf\n", tr.load_a, tr.load_b,
			tr.multiplication);

//	vm_c->f.tofile(vm_c, output);

	vm_a->f.free(vm_a);
	vm_b->f.free(vm_b);
	vm_c->f.free(vm_c);

	return EXIT_SUCCESS;
}
