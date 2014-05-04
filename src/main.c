/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013,2014
 * https://github.com/nesro/sparse-matrices
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <getopt.h>

#include "utils.h"
#include "virtual_matrix.h"
#include "csr_matrix.h"
#include "coo_matrix.h"
#include "bsr_matrix.h"
#include "den_matrix.h"
#include "kat_matrix.h"
#include "qdt_matrix.h"

static void usage() {
	printf("Usage: ./main\n"
			"-f [coo|csr|bsr|den|kat]\n"
			"-a <matrix_a>\n"
			"-b <matrix_b>\n"
			"-o [stdout|file]\n"
			"-s block_size\n"
			"-v print informations\n"
			"-V second matrix is a vector\n"
			"Formats:\n"
			"Matrices a and b has to be in the .mtx MatrixMarket file format.\n"
			".mtx format: math.nist.gov/MatrixMarket/‎\n"
			"See: github.com/nesro/sparse-matrices\n");
}

static vm_type_t parse_format(const char *format_string) {

	if (0) {
		(void) 0;
	} else if (strcmp(format_string, "qdt") == 0) {
		return QDT;
	} else if (strcmp(format_string, "den") == 0) {
		return DEN;
	} else if (strcmp(format_string, "csr") == 0) {
		return CSR;
	} else if (strcmp(format_string, "kat") == 0) {
		return KAT;
	} else if (strcmp(format_string, "bsr") == 0) {
		return BSR;
	} else if (strcmp(format_string, "coo") == 0) {
		return COO;
	} else {
		fprintf(stderr, "Unknown matrix format.\n");
		exit(EXIT_FAILURE);
	}

	/*NOTREACHABLE*/
	return UNKNOWN;
}

int load_block_size(const char *string) {

	char *tmp = NULL;
	int block_size;

	block_size = strtol(string, &tmp, 10);

	if (tmp == NULL) {
		fprintf(stderr, "ERROR: Cannot load block size. Invalid input.\n");
		exit(EXIT_FAILURE);
	}

	if (block_size < 1) {
		fprintf(stderr, "ERROR: Minimum of block size is 1.\n");
		exit(EXIT_FAILURE);
	}

	if (!is_power_of_two(block_size)) {
		fprintf(stderr, "ERROR: Leaf size is not power of two.\n");
		exit(EXIT_FAILURE);
	}

	return block_size;
}

int main(int argc, char *argv[]) {

	int c;
	int index;
	vm_type_t format = UNKNOWN;
	char *matrix_a = NULL;
	char *matrix_b = NULL;
	int block_size = -1;
	char *output = NULL; /* -o flag */
	vm_t *vm_a = NULL;
	vm_t *vm_b = NULL;
	vm_t *vm_c = NULL;
	double time_mul;
	int verbose = 0;
	int mul_vector = 0;

	if (argc < 2) {
		fprintf(stderr,
				"Error: Not enough parameters. Run \"%s -h\" for more informations.\n",
				argv[0]);
		return EXIT_SUCCESS;
	}

	opterr = 0;
	while ((c = getopt(argc, argv, "a:b:f:hs:o:tvV")) != -1) {
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
		case 's':
			block_size = load_block_size(optarg);
			break;
		case 'o':
			output = optarg;
			break;
		case 'v':
			verbose = 1;
			break;
		case 'V':
			mul_vector = 1;
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
		fprintf(stderr, "ERROR: Either matrix a and b are not set.\n");
		usage();
		return EXIT_FAILURE;
	}

	if (format == UNKNOWN) {
		fprintf(stderr, "ERROR: Format is not set.\n");
		usage();
		return EXIT_FAILURE;
	}

	if (matrix_b != NULL && matrix_a == NULL) {
		matrix_a = matrix_b;
		matrix_b = NULL;
	}
	if (matrix_b != NULL && matrix_a == NULL)
		matrix_b = matrix_a;

	if (vm_has_blocks(format) && block_size == -1) {
		fprintf(stderr, "ERROR: Block size is not set or has an invalid value."
				"Set it with an -s option to a power of two.\n");
		return EXIT_FAILURE;
	}

	if (vm_has_blocks(format))
		vm_load_mm(&vm_a, format, matrix_a, block_size);
	else
		vm_load_mm(&vm_a, format, matrix_a);

	if (matrix_b != NULL) {
		if (mul_vector) {
			vm_load_mm(&vm_b, VEC, matrix_b, vm_a->w);
		} else {
			if (vm_has_blocks(format))
				vm_load_mm(&vm_b, format, matrix_b, block_size);
			else
				vm_load_mm(&vm_b, format, matrix_b);
		}
	}

	if (matrix_b != NULL)
		time_mul = vm_a->f.mul(vm_a, vm_b, &vm_c, NAIVE);
	else
		time_mul = vm_a->f.mul(vm_a, vm_a, &vm_c, NAIVE);

	if (output != NULL)
		vm_c->f.mm_save(vm_c, output);

	if (verbose) {
		printf("time_mul %lf\n", time_mul);
		printf("a_size %zu\n", vm_a->object_size);

		if (format == KAT) {
			printf("kat_n %d\n", KAT_N);
			printf("kat_sm_size %d\n", ((kat_matrix_t*) vm_a)->sm_size);
			printf("kat_a_inner %d\n", ((kat_matrix_t*) vm_a)->nodes_inner);
			printf("kat_a_dense %d\n", ((kat_matrix_t*) vm_a)->nodes_den);
			printf("kat_a_csr %d\n", ((kat_matrix_t*) vm_a)->nodes_csr);
		}
	}

	vm_a->f.free(vm_a);

	if (matrix_b != NULL)
		vm_b->f.free(vm_b);

	vm_c->f.free(vm_c);

	return EXIT_SUCCESS;
}
