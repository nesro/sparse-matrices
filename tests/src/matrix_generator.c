/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
 * https://github.com/nesro/sparse-matrices
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>

double sparsity = 1.;
int diagonals = 0;
int center_diagonal = 0;
const char *output = NULL;
int verbose = 0;
int n = 0;
int w = 0;
int h = 0;
int random = 0;
int float_values = 0;
double start = 2.;

/* TODO: */
double random_double() {

	double d = 0;

	return d;
}

static void usage(const char *argv0) {
	printf("%s [OPTIONS]\n"
			" -c              add a center diagonal\n"
			" -d diagonals    number of diagonals\n"
			" -f              values will be float numbers\n"
			" -h              print this help\n"
			" -H height       height of a matrix\n"
			" -n size         width = height\n"
			" -o file         output file, default is stdout\n"
			" -r              values will be random\n"
			" -s sparsity     how many %% of matrix will be filled\n"
			" -S start        starting value"
			" -v              verbose mode\n"
			" -w width        width of a matrix\n", argv0);
}

static void generate_dense(FILE* f) {

	int i;
	int j;
	double v = start;

	if (verbose)
		printf("Generating dense matrix.\n");

	fprintf(f, "%%%%MatrixMarket matrix coordinate real general\n");
	fprintf(f, "%d %d %d\n", n, n, n * n);

	for (i = 1; i <= n; i++) {
		//printf("%04d/%04d\n", i, n);

		for (j = 1; j <= n; j++)
			fprintf(f, "%d %d %lf\n", i, j, v++);
	}

}

int main(int argc, char *argv[]) {

	int c;
	int index;

	srand(time(NULL));

	opterr = 0;
	while ((c = getopt(argc, argv, "cf:hn:o:s:S:v")) != -1) {
		switch (c) {
		case 'c':
			center_diagonal = 1;
			break;
		case 'd':
			diagonals = atoi(optarg);
			break;
		case 'h':
			usage(argv[0]);
			return EXIT_SUCCESS;
		case 'n':
			n = atoi(optarg);
			w = n;
			h = n;
			break;
		case 'o':
			output = optarg;
			break;
		case 's':
			sparsity = (atoi(optarg) / ((double) 100));
			break;
		case 'S':
			start = atof(optarg);
			break;
		case 'v':
			verbose = 1;
			break;
		case '?':
			if (optopt == 'o')
				fprintf(stderr,
						"ERROR: Filename is missing for the option -o.\n");
			else if (optopt == 'd')
				fprintf(stderr, "ERROR: Specify the count of diagonals.\n");
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

	/**********************************************************************/

	FILE *f;

	if (output != NULL)
		f = fopen(output, "w");
	else
		f = stdout;

	if (f == NULL) {
		printf("Error opening file: %s\n",
				(output == NULL) ? "stdout" : output);
		exit(1);
	}

	if (verbose) {
		printf("TODO: some useful informations\n");
	}

	if (sparsity == 1.) {
		generate_dense(f);
		goto done;
	}

	done: /**/
	if (output != NULL)
		fclose(f);
	return EXIT_SUCCESS;
}
