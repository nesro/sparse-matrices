/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
 * https://github.com/nesro/sparse-matrices
 */

#include "virtual_matrix.h"
#include "den_matrix.h"
#include "bsr_matrix.h"
#include "coo_matrix.h"
#include "csr_matrix.h"
#include "qdt_matrix.h"
#include "kat_matrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

void vm_create(vm_t **vm, vm_type_t type, ...) {

	va_list vl;

	va_start(vl, 0);

	switch (type) {
	case DEN:
		den_vm_init((den_matrix_t **) vm, vl);
		break;
	case BSR:
		bsr_vm_init((bsr_t **) vm, vl);
		break;
	case COO:
		coo_vm_init((coo_matrix_t **) vm, vl);
		break;
	case CSR:
		csr_vm_init((csr_t **) vm, vl);
		break;
	case QDT:
		qdt_vm_init((qdt_matrix_t **) vm, vl);
		break;
	case KAT:
		kat_vm_init((kat_matrix_t **) vm, vl);
		break;
	default:
		break;
	}

	va_end(vl);
}

void vm_load_mm(vm_t **vm, vm_type_t type, const char *mm, ...) {

	va_list vl;

	va_start(vl, 0);

	switch (type) {
	case DEN:
		den_from_mm((den_matrix_t **) vm, mm, vl);
		break;
	case BSR:
		bsr_from_mm((bsr_t **) vm, mm, vl);
		break;
	case COO:
		coo_from_mm((coo_matrix_t **) vm, mm, vl);
		break;
	case CSR:
		csr_from_mm((csr_t **) vm, mm, vl);
		break;
	case QDT:
		qdt_from_mm((qdt_matrix_t **) vm, mm, vl);
		break;
	case KAT:
		kat_from_mm((kat_matrix_t **) vm, mm, vl);
		break;
	default:
		break;
	}

	va_end(vl);
}

void vm_print(vm_t *vm) {

	printf("width: %d\n", vm->w);
	printf("height: %d\n", vm->h);

	printf("type: ");
	switch (vm->type) {
	case DEN:
		printf("dense\n");
		break;
	case BSR:
		printf("bsr\n");
		break;
	case CSR:
		printf("csr\n");
		break;
	case QDT:
		printf("qdt\n");
		break;
	case KAT:
		printf("kat\n");
		break;
	case COO:
		printf("kat\n");
		break;
	default:
		printf("<unknown or not defined>\n");
		break;
	}
}

void vm_exec(action_t action, vm_type_t type_a, vm_type_t type_b,
		const char *file_a, const char *file_b, vm_t **c, char flag, ...) {

//	vm_exec(MUL, CSR, DEN, "file1.mtx", "file2.mtx", c, MATRIX_A, LEAF_SIZE, 12);

	vm_t *a;
	vm_t *b;
	va_list va_list_a;
	va_list va_list_b;

	/* is va empty */

	/* is MATRIX_A */
	/* while token != MATRIX_A_END */

	/* is MATRIX_B */
	/* while token != MATRIX_B_END */

	vm_load_mm(&a, type_a, file_a, va_list_a);
	vm_load_mm(&b, type_b, file_b, va_list_b);

	a->f.mul(a, b, c, flag);
}
