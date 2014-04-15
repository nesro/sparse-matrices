/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
 * https://github.com/nesro/sparse-matrices
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <errno.h>

/*
 * strdup is not a C standard function
 */
char *strdup(const char *s) {
	int n = strlen(s) + 1;
	char *d = malloc(n);
	if (d)
		(void) strcpy(d, s);
	return d;
}

/*
 * Because this is a stand-alone file, we can use a lot of global variables.
 */
double sparsity = 1.;
int diagonals = 0;
int center_diagonal = 0;
const char *output = NULL;
int verbose = 0;
int n = 0;
int w = -1;
int h = -1;
int random = 0;
int float_values = 0;
double g_start = 2.;

/************************************************************ position vector */

typedef struct position {
	int y;
	int x;
} position_t;

typedef struct mtx {
	position_t *data;
	int len;
	int items;
	int width;
	int height;
} mtx_t;

static mtx_t *mtx_init(int height, int width, double sparsity) {

	mtx_t *v;
	v = malloc(sizeof(mtx_t));
	v->items = 0;
	v->height = height;
	v->width = width;
	v->len = height * width * sparsity + 1;
	v->data = malloc(v->len * sizeof(position_t));
	return v;
}

static void mtx_free(mtx_t *mtx) {

	if (mtx == NULL)
		return;

	free(mtx->data);
	free(mtx);
}

static void mtx_add(mtx_t *mtx, int y, int x) {

	if (mtx->items >= mtx->len - 1) {

		printf("Reallocing...\n");

		position_t *new = NULL;
		mtx->len *= 2;

		new = realloc(mtx->data, mtx->len * sizeof(position_t));

		if (new != NULL) {
			mtx->data = new;
		} else {
			free(mtx->data);
			free(mtx);
			fprintf(stderr, "position_v_add: realloc %s", strerror(errno));
			exit(1);
		}
	}

	mtx->data[mtx->items].y = y;
	mtx->data[mtx->items].x = x;
	mtx->items++;
}

static int pos_cmp(const void *a, const void *b) {
	position_t *pa = ((position_t *) a);
	position_t *pb = ((position_t *) b);

	if (pa->y > pb->y)
		return 1;
	else if (pa->y < pb->y)
		return 0;
	else
		return (pa->x > pb->x);
}

static void mtx_uniqsort(mtx_t *mtx) {

	int i;
	int j;

	qsort(mtx->data, mtx->items, sizeof(int), pos_cmp);

	for (i = j = 0; i < mtx->items; i++)
		if (mtx->data[i].x != mtx->data[j].x
				&& mtx->data[i].y != mtx->data[j].y)
			mtx->data[++j] = mtx->data[i];

	mtx->items = j;

}

static void mtx_write(mtx_t *mtx, FILE *f) {

	int i;
	double v = g_start;

	//mtx_uniqsort(mtx);

	fprintf(f, "%%%%MatrixMarket matrix coordinate real general\n");
	fprintf(f, "%d %d %d\n", mtx->width, mtx->height, mtx->items);

	for (i = 0; i < mtx->items; i++) {
		fprintf(f, "%d %d %lf\n", mtx->data[i].y, mtx->data[i].x, v++);
	}
}

static void mtx_sparse_fill(mtx_t *mtx, int ay, int ax, int by, int bx,
		double sparsity) {

	int i;
	int j;
	int w = bx - ax;
	int h = by - ay;
	int space = w * h;
	int items = space * sparsity;
	double prob = (double) items / (double) space;

	printf("w=%d,h=%d,space=%d,items=%d,prob=%lf\n", w, h, space, items, prob);
	fflush(stdout);

	for (i = 0; i < h; i++) {
		for (j = 0; j < w; j++) {

//			printf("w=%d,h=%d,space=%d,items=%d,prob=%lf\n", w, h, space, items,
//					prob);
//			fflush(stdout);

			if (items == space
					|| ((double) rand() / (double) RAND_MAX) < prob) {
				items--;
				prob = (double) items / (double) space;

				mtx_add(mtx, i + 1, j + 1);

			}

			space--;

			if (items == 0) {
				//			return;
			}
		}
	}

	printf("items=%d\n", items);
	assert(items == 0);
}

static void mtx_sparse_blocks(mtx_t *mtx, char *blocks) {

	int ay;
	int ax;
	int by;
	int bx;
	int sparsity;
	char *tmp;
	char* end;

	for (tmp = strtok(blocks, ",");; tmp = strtok(NULL, ",")) {

		if (tmp == NULL)
			goto bad_format;

		ay = strtol(tmp, &end, 10);
		if (*end)
			goto bad_format;

		tmp = strtok(NULL, ",");
		ax = strtol(tmp, &end, 10);
		if (*end)
			goto bad_format;

		tmp = strtok(NULL, ",");
		by = strtol(tmp, &end, 10);
		if (*end)
			goto bad_format;

		tmp = strtok(NULL, ",");
		bx = strtol(tmp, &end, 10);
		if (*end)
			goto bad_format;

		tmp = strtok(NULL, ",");
		sparsity = strtol(tmp, &end, 10);
		if (*end)
			goto bad_format;

		mtx_sparse_fill(mtx, ay, ax, by, bx, sparsity);
	}

	return;

	bad_format: /**/
	fprintf(stderr, "Bad block format. Use -h for help.\n");
	exit(1);
}

/******************************************************************************/

static void usage(const char *argv0) {
	printf("%s [OPTIONS]\n"
			" -b              blocks in format: ay,ax,by,bx,sparsity,..."
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
			" -W width        width of a matrix\n", argv0);
}

/*
 * Fast function for generating 100% dense matrices.
 */
static void generate_dense(FILE *f) {

	int i;
	int j;
	double v = g_start;

	if (verbose)
		printf("Generating dense matrix.\n");

	fprintf(f, "%%%%MatrixMarket matrix coordinate real general\n");
	fprintf(f, "%d %d %d\n", n, n, n * n);

	for (i = 1; i <= n; i++)
		for (j = 1; j <= n; j++)
			fprintf(f, "%d %d %lf\n", i, j, v++);
}

int main(int argc, char *argv[]) {

	int c;
	int index;
	FILE *f;
	mtx_t *mtx = NULL;
	char *blocks = NULL;

	srand((unsigned) time(NULL));

	opterr = 0;
	while ((c = getopt(argc, argv, "cf:hH:n:o:s:S:vW:")) != -1) {
		switch (c) {
		case 'b':
			blocks = strdup(optarg);
			break;
		case 'c':
			center_diagonal = 1;
			break;
		case 'd':
			diagonals = atoi(optarg);
			break;
		case 'h':
			usage(argv[0]);
			return EXIT_SUCCESS;
		case 'H':
			h = atoi(optarg);
			break;
		case 'n':
			n = atoi(optarg);
			w = n;
			h = n;
			break;
		case 'o':
			output = optarg;
			break;
		case 's':
			sscanf(optarg, "%lf", &sparsity);
			break;
		case 'S':
			g_start = atof(optarg);
			break;
		case 'v':
			verbose = 1;
			break;
		case 'W':
			w = atoi(optarg);
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

	if (w < 0 || h < 0) {
		fprintf(stderr, "Width and height must be set and positive.\n");
		return EXIT_FAILURE;
	}

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

	/*
	 *  our matrix will be sparse. we need **pos now.
	 */
	mtx = mtx_init(h, w, sparsity);

	/* general sparsity */
	mtx_sparse_fill(mtx, 0, 0, h, w, sparsity);

	if (blocks != NULL)
		mtx_sparse_blocks(mtx, blocks);

	mtx_write(mtx, f);

	done: /**/
	if (output != NULL)
		fclose(f);
	mtx_free(mtx);
	free(blocks);

	return EXIT_SUCCESS;
}
