/**
 * Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
 * https://github.com/nesro/sparse-matrices
 */

#include <assert.h>
#include <limits.h>
#include <stdarg.h>

#ifndef UTILS_H_
#define UTILS_H_

/*
 * Warning: this can not be changed alone. All unrolled functions must be
 * updated too.
 */
#define UNROLL 8

#ifdef PRINT_DEBUG
#define _PRINT_DEBUG 1
#else
#define _PRINT_DEBUG 0
#endif

/* variable argument list handling *****************************************/

#define VA_END -9999
#define VA_DEFAULT -9998

int va_get_int(va_list va, int default_value, int *va_flag);
double va_get_double(va_list va, double default_value, int *va_flag);

void error(const char *fmt, ...);

/* debugging ***************************************************************/

/*
 * status debug: if the status is true, the message is print
 */
#define _s_debug(status,fmt) do { \
	if (status) \
		fprintf(stdout, "d:%s:%s:%d(): " fmt, __FILE__ , \
			__func__ , __LINE__ ); \
		fflush(stdout); \
	} while (0)
#define _s_debugf(status,fmt, ...) do { \
	if (status) \
		fprintf(stdout, "d:%s:%s:%d(): " fmt "", __FILE__ , \
			__func__ , __LINE__ , __VA_ARGS__ ); \
		fflush(stdout); \
	} while (0)


/***************************************************************************/

#define _debug(fmt, ...) do { \
	if (1 || _PRINT_DEBUG) \
		fprintf(stderr, "d:%s:%s:%d(): " fmt "\n", __FILE__ , \
			__func__ , __LINE__ , __VA_ARGS__ ); \
		fflush(stderr); \
	} while (0)

#define _debug_var(var) do { \
	if (1 || _PRINT_DEBUG) \
		fprintf(stderr, "d:%s:%s:%d(): " #var "=%d\n", __FILE__ , \
			__func__ , __LINE__ , var ); \
		fflush(stderr); \
	} while (0)

//void die(const char *fmt, ...);

#define die(fmt) do { \
	if (1 || _PRINT_DEBUG) \
		fprintf(stderr, "DIE:%s:%s:%d(): " fmt "\n", __FILE__ , \
			__func__ , __LINE__); \
		fflush(stderr); \
		abort(); \
	} while (0)

#define fdie(fmt, ...) do { \
	if (1 || _PRINT_DEBUG) \
		fprintf(stderr, "DIE:%s:%s:%d(): " fmt "\n", __FILE__ , \
			__func__ , __LINE__ , __VA_ARGS__ ); \
		fflush(stderr); \
		abort(); \
	} while (0)

/***************************************************************************/

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
