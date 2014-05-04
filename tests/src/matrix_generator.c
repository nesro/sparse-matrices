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
#include <math.h>

#define MATRIX_GENERATOR_DEBUG 1

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

static inline int max(int a, int b) {
	if (a > b)
		return a;
	else
		return b;
}

/*
 * Because this is a stand-alone file, we can use a lot of global variables.
 */
int is_symmetric = 0; /* symmetric or general */
double sparsity = 1.;
int has_diagonal = 1;
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

	/*
	 * Wolfram Mathematica highlights big items. We want to highlignt diagonals
	 * for example.
	 */
	int is_big;
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
	v->len = (4 * height * width * sparsity) + 1000;
	v->data = malloc(v->len * sizeof(position_t));
	return v;
}

static void mtx_free(mtx_t *mtx) {

	if (mtx == NULL)
		return;

	free(mtx->data);
	free(mtx);
}

static void mtx_add(mtx_t *mtx, int y, int x, int is_big) {

	if (0 && MATRIX_GENERATOR_DEBUG)
		printf("mtx_add y=%d, x=%d\n", y, x);

	/*
	 * Is item out of bounds?
	 */
	if (y < 0 || y > mtx->height - 1 || x < 0 || x > mtx->width - 1)
		return;

	/*
	 * If the matrix is symmetric, check out, if the item is under the diagonal.
	 */
	if (is_symmetric && !(y <= max(y, x) && x <= max(y, x))) {
		fprintf(stderr, "item at y=%d, x=%d is above main diagonal"
				" although the matrix is symmetric\n", y, x);
	}

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
	mtx->data[mtx->items].is_big = is_big;
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

	qsort(mtx->data, mtx->items, sizeof(position_t), pos_cmp);

	for (i = j = 0; i < mtx->items; i++) {
		if (mtx->data[i].y == mtx->data[j].y
				&& mtx->data[i].x == mtx->data[j].x)
			continue;

		mtx->data[++j] = mtx->data[i];
	}

	mtx->items = j;

}

static void mtx_write(mtx_t *mtx, FILE *f) {

	int i;
	double v = g_start;

	if (f != stdout)
		printf("Sorting...\n");

	mtx_uniqsort(mtx);

	fprintf(f, "%%%%MatrixMarket matrix coordinate real %s\n",
			(is_symmetric) ? "symmetric" : "general");
	fprintf(f, "%d %d %d\n", mtx->width, mtx->height, mtx->items);

	if (f != stdout)
		printf("Writing...\n");

	for (i = 0; i < mtx->items; i++) {
		if (mtx->data[i].is_big) {
			fprintf(f, "%d %d %lf\n", mtx->data[i].y, mtx->data[i].x,
					v + 10);
		} else {
			fprintf(f, "%d %d %lf\n", mtx->data[i].y, mtx->data[i].x, v);
		}
	}
}

static void mtx_diagonal(mtx_t *mtx, int ay, int ax, int by, int bx,
		double sparsity) {

	double i;
	double steps = sqrt(
			pow((double) (ay - by), 2) + pow((double) (ax - bx), 2));

	double step_y = abs(ay - by) / steps;
	double step_x = abs(ax - bx) / steps;

	double curr_y = ay;
	double curr_x = ax;
	/* TODO: space probabilty */
	/* double prob = (double) 1 / steps * sparsity; */

	if (MATRIX_GENERATOR_DEBUG)
		printf(
				"mtx_diagonal: ay=%d, ax=%d, steps=%lf, step_y=%lf, step_x=%lf, prob=%lf\n",
				ay, ax, steps, step_y, step_x, sparsity);

	for (i = 0; i < steps; i++) {
		curr_y += step_y;
		curr_x += step_x;

		if (((double) rand() / (double) RAND_MAX) < sparsity) {
			mtx_add(mtx, 1 + ((int) curr_y), 1 + ((int) curr_x), 1);
		}
	}

	mtx_add(mtx, ay + 1, ax + 1, 1);
	mtx_add(mtx, by + 1, bx + 1, 1);
}

static void mtx_fill(mtx_t *mtx, int ay, int ax, int by, int bx,
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

	assert(ay < by);
	assert(ax < bx);

	for (i = ay; i < by; i++) {
		for (j = ax; j < bx; j++) {

//			printf("w=%d,h=%d,space=%d,items=%d,prob=%lf\n", w, h, space, items,
//					prob);
//			fflush(stdout);

			if (items == space
					|| ((double) rand() / (double) RAND_MAX) < prob) {
				items--;
				prob = (double) items / (double) space;

				mtx_add(mtx, i + 1, j + 1, 0);

			}

			space--;

			if (items == 0) {
				//			return;
			}
		}
	}

	assert(items == 0);
}

static void mtx_sparse_items(mtx_t *mtx, char *blocks) {

	enum {
		/* ay, ax, by, bx, sparsity */
		BLOCK, /**/
		MIRRORED_BLOCK, /**/

		/* y, x, h, w, sparsity */
		BLOCK_WH, /**/
		MIRRORED_BLOCK_WH, /**/

		/* ay, ax, by, bx, sparsity */
		DIAGONAL, /**/
		MIRRORED_DIAGONAL, /**/

		/* count, size_probability, h, w, sparsity */
		RANDOM_BLOCKS, /**/
		MIRRORED_RANDOM_BLOCKS, /**/

		DIAGONAL_BLOCKS, /**/
	} type;

	int ay;
	int ax;
	int by;
	int bx;
	double sparsity;
	char *tmp;
	char* end;

	int i;
	int tmp_y;
	int tmp_x;
	int tmp_size;

	for (tmp = strtok(blocks, ",");; tmp = strtok(NULL, ",")) {

		if (tmp == NULL)
			return;

		if (strcmp(tmp, "block") == 0) {
			type = BLOCK;
		} else if (strcmp(tmp, "mblock") == 0) {
			type = MIRRORED_BLOCK;
		} else if (strcmp(tmp, "blockwh") == 0) {
			type = BLOCK_WH;
		} else if (strcmp(tmp, "mblockwh") == 0) {
			type = MIRRORED_BLOCK_WH;
		} else if (strcmp(tmp, "diagonal") == 0) {
			type = DIAGONAL;
		} else if (strcmp(tmp, "mdiagonal") == 0) {
			type = MIRRORED_DIAGONAL;
		} else if (strcmp(tmp, "rblocks") == 0) {
			type = RANDOM_BLOCKS;
		} else if (strcmp(tmp, "mrblocks") == 0) {
			type = MIRRORED_RANDOM_BLOCKS;
		} else if (strcmp(tmp, "diablocks") == 0) {
			type = DIAGONAL_BLOCKS;
		} else {
			fprintf(stderr, "matrix_sparse_items: not a valid item\n");
			goto bad_format;
		}
		tmp = strtok(NULL, ",");
		if (tmp == NULL) {
			fprintf(stderr,
					"matrix_sparse_items: Missing the first parameter!\n");
			goto bad_format;
		}

		ay = strtol(tmp, &end, 10);
		if (*end)
			goto bad_format;

		tmp = strtok(NULL, ",");
		if (tmp == NULL) {
			fprintf(stderr,
					"matrix_sparse_items: Missing the second parameter!\n");
			goto bad_format;
		}
		ax = strtol(tmp, &end, 10);
		if (*end)
			goto bad_format;

		tmp = strtok(NULL, ",");
		if (tmp == NULL) {
			fprintf(stderr,
					"matrix_sparse_items: Missing the third parameter!\n");
			goto bad_format;
		}
		by = strtol(tmp, &end, 10);
		if (*end)
			goto bad_format;

		tmp = strtok(NULL, ",");
		if (tmp == NULL) {
			fprintf(stderr,
					"matrix_sparse_items: Missing the fourth parameter!\n");
			goto bad_format;
		}
		bx = strtol(tmp, &end, 10);
		if (*end)
			goto bad_format;

		tmp = strtok(NULL, ",");
		if (tmp == NULL) {
			fprintf(stderr, "matrix_sparse_items: Missing sparsity!\n");
			goto bad_format;
		}
		sparsity = strtod(tmp, &end);
		if (*end)
			goto bad_format;

		switch (type) {
		case BLOCK:
			mtx_fill(mtx, ay, ax, by, bx, sparsity);
			break;
		case MIRRORED_BLOCK:
			mtx_fill(mtx, ay, ax, by, bx, sparsity);
			mtx_fill(mtx, ax, ay, bx, by, sparsity);
			break;
		case BLOCK_WH:
			mtx_fill(mtx, ay, ax, ay + by, ax + bx, sparsity);
			break;
		case MIRRORED_BLOCK_WH:
			mtx_fill(mtx, ay, ax, ay + by, ax + bx, sparsity);
			mtx_fill(mtx, ax, ay, ax + bx, ay + by, sparsity);
			break;
		case DIAGONAL:
			mtx_diagonal(mtx, ay, ax, by, bx, sparsity);
			break;
		case MIRRORED_DIAGONAL:
			mtx_diagonal(mtx, ay, ax, by, bx, sparsity);
			mtx_diagonal(mtx, ax, ay, bx, by, sparsity);
			break;
		case RANDOM_BLOCKS:
		case MIRRORED_RANDOM_BLOCKS:
			for (i = 0; i < ay; i++) {

				/*
				 * If we're constructing a symmetric matrix, we want the
				 * block under the main diagonal.
				 */
				if (is_symmetric) {
					do {
						tmp_y = rand() % (mtx->height - by - 1);
						tmp_x = rand() % (mtx->width - bx - 1);
					} while (!(tmp_y <= max(tmp_y, tmp_x)
							&& tmp_x <= max(tmp_y, tmp_x)));
				} else { /* !is_symmetric */
					tmp_y = rand() % (mtx->height - by - 1);
					tmp_x = rand() % (mtx->width - bx - 1);
				}

				/* make the block on more friendly position */
				if (1) {
					tmp_y -= tmp_y % by;
					tmp_x -= tmp_x % bx;
				}

				if (ax != 0)
					tmp_size = rand() % ax;
				else
					tmp_size = 0;

				mtx_fill(mtx, tmp_y, tmp_x, tmp_y + bx + tmp_size,
						tmp_x + by + tmp_size, sparsity);

				if (type == MIRRORED_RANDOM_BLOCKS) {
					mtx_fill(mtx, tmp_x, tmp_y, tmp_x + bx + tmp_size,
							tmp_y + by + tmp_size, sparsity);
				}
			}
			break;
		case DIAGONAL_BLOCKS:
			for (i = 0; i < mtx->height; i += ay) {
				mtx_fill(mtx, i, i, i + ay, i + ay, 1);
			}
			break;
		default:
			break;
		}
	}

	return;

	bad_format: /**/
	fprintf(stderr,
			"matrix_sparse_items: Bad block format. Use -h for help.\n");
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
			" -m              matrix is symmetric (mirrored)\n"
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
	fprintf(f, "%d %d %d\n", h, w, w * h);

	for (i = 1; i <= h; i++)
		for (j = 1; j <= w; j++)
			fprintf(f, "%d %d %lf\n", i, j, v++);
}

int main(int argc, char *argv[]) {

	int c;
	int index;
	FILE *f;
	mtx_t *mtx = NULL;
	char *items = NULL;

	srand((unsigned) time(NULL));

	opterr = 0;
	while ((c = getopt(argc, argv, "cf:hH:i:n:mo:s:S:vW:")) != -1) {
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
		case 'H':
			h = atoi(optarg);
			break;
		case 'i':
			items = strdup(optarg);
			break;
		case 'n':
			n = atoi(optarg);
			w = n;
			h = n;
			break;
		case 'm':
			is_symmetric = 1;
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

	if (has_diagonal) {
		mtx_diagonal(mtx, 0, 0, h, w, 1.);
	}

	/* general sparsity */
	if (sparsity > 0)
		mtx_fill(mtx, 0, 0, h, w, sparsity);

	if (items != NULL)
		mtx_sparse_items(mtx, items);

	mtx_write(mtx, f);

	done: /**/
	if (output != NULL)
		fclose(f);
	mtx_free(mtx);
	free(items);

	return EXIT_SUCCESS;
}
