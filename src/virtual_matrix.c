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
	case VEC:
		vec_from_mm((vec_t **) vm, mm, vl);
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
	case VEC:
		printf("vec\n");
		break;
	default:
		printf("<unknown or not defined>\n");
		break;
	}
}

int vm_has_blocks(vm_type_t type) {
	switch (type) {
	case BSR:
	case KAT:
		return 1;
	case COO:
	case CSR:
	case DEN:
		return 0;
	default:
		fdie("Unknown format: %d\n", type);
		/*NOTREACHED*/
		return -1;
	}
}

