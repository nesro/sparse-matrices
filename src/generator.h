/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
 * https://github.com/nesro/sparse-matrices
 */

#include "csr_matrix.h"
#include "den_matrix.h"

#ifndef GENERATOR_H_
#define GENERATOR_H_

//TODO: move to dense_matrix.h
void generate_dense(den_matrix_t *dense, int n, int density);

void generate_mm(const char* filename, int width, int height, int density);

#endif /* GENERATOR_H_ */
