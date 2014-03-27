#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {

	int n;
	int i;
	int j;
	double v;
	FILE *f;

	n = atoi(argv[2]); /* FIXME: strol */
	v = atoi(argv[3]);

	f = fopen(argv[1], "w");
	if (f == NULL) {
		printf("Error opening file!\n");
		exit(1);
	}

	fprintf(f, "%%%%MatrixMarket matrix coordinate real general\n");
	fprintf(f, "%d %d %d\n", n, n, n * n);

	for (i = 1; i <= n; i++) {
		printf("%04d/%04d\n", i, n);

		for (j = 1; j <= n; j++)
			fprintf(f, "%d %d %lf\n", i, j, v++);
	}

	fclose(f);
	return 0;
}
