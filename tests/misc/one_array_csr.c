
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

/**
 * n: width and height
 * nnz: number of nonzero elements
 * data: rows, [ci|value]
 */
typedef struct csr {
	int n;
	int nnz;
	void *data;
} csr_t;

void csr_mmm(csr_t *a, csr_t *b, void *c) {


}

void test_csr(const char *mm_file) {




}


int main(int argc, char *argv[]) {

}
