/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
 * https://github.com/nesro/sparse-matrices
 */

#include <assert.h>

#ifndef UTILS_H_
#define UTILS_H_

#ifdef OMP_THREADS
#	define __OMP_NUM_THREADS__ num_threads(OMP_THREADS)
#else
#	define __OMP_NUM_THREADS__
#endif

#ifdef MMLUR
#define __CSR_MMM_LOOP_UNROLLING__ 1
#else
#define __CSR_MMM_LOOP_UNROLLING__ 0
#endif

#ifdef MVLUR
#define __CSR_MVM_LOOP_UNROLLING__ 1
#else
#define __CSR_MVM_LOOP_UNROLLING__ 0
#endif

typedef enum matrix_type {
	DENSE, /**/
	COO, /**/
	CSR, /**/
	QUADTREE, /**/
} matrix_type_t;

typedef struct time_record {
	double multiplication;
	double load_a;
	double load_b;
} time_record_t;

/* dynamic array **************************************************************/

/*
 typedef struct da {
 char *data;
 int size;
 int length;
 int default_value;
 } da_t;

 void da_init(da_t *da, size_t element_size, int size, int default_value);
 void da_resize(da_t *da);
 void da_free(da_t *da);

 #define	da_read(arr, type, index) (assert(index < (arr)->length), \
	((type*)(arr)->data)[(index)])

 #define da_append(arr, type, value) do { \
		((arr)->length + 1 >= (arr)->size) ? da_resize((arr)); \
		(type*)(arr)->data)[(arr)->length++] = value; \
	} while(0)*/

/******************************************************************************/

/**
 * Failure indicator.
 */
typedef int error_t;
#define ERROR_NONE 0 /* no error */
#define ERROR_NOMEMORY 1 /* malloc/calloc fail */

#define DATATYPE_FORMAT "%lf"
typedef double datatype_t;

#define MAKE_HTML 0

/* size of the <td> element in html output */
#define TD_SIZE 15

#define LOAD_VALUE 1

double random_double(void);
int possibility(double n);
int rand_is_nonzero(int space, int items);
int is_power_of_two(unsigned int x);

/* http://c.learncodethehardway.org/ */

#ifdef NDEBUG
#define debug(M, ...)
#else
#define debug(M, ...) fprintf(stderr, "DEBUG %s:%d: " M "\n", __FILE__, __LINE__, ##__VA_ARGS__)
#endif

#define clean_errno() \
	(errno == 0 ? "None" : strerror(errno))

#define log_err(M, ...) \
	fprintf(stderr, "[ERROR] (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, \
		clean_errno(), ##__VA_ARGS__)

#define log_warn(M, ...) \
	fprintf(stderr, "[WARN] (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, \
		clean_errno(), ##__VA_ARGS__)

#endif /* UTILS_H_ */
