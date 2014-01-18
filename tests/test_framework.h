/**
 * There are hundreds of c test like this. Based on macro functions. It's
 * faster to write my own version than get used to someone's else version.
 *
 * It's under development. This is not a final version.
 */

#ifndef TEST_UTILS_H_
#define TEST_UTILS_H_

#define BEGIN do {
#define END } while(0)

#define COLORS 1

#define EXIT_AFTER_FAIL 1
#define TEST_MATRICES_DIR "./tests/matrices"

/******************************************************************************/

#if COLORS == 1
#	define COLOR_RED "\x1B[31m"
#	define COLOR_GREEN "\x1B[32m"
#	define COLOR_YELLOW "\x1b[33m"
#	define COLOR_BLUE "\x1b[34m"
#	define COLOR_MAGNETA "\x1b[35m"
#	define COLOR_CYAN "\x1b[36m"
#	define COLOR_RESET "\x1B[0m"
#
#	define OK COLOR_GREEN "[OK]\n" COLOR_RESET
#	define FAIL COLOR_RED "[FAIL]\n" COLOR_RESET
#else /* COLORS == 0 */
#	define OK "[OK]\n"
#	define FAIL "[FAIL]\n"
#endif /* COLORS */

#define PRINT_OK() BEGIN \
	fprintf(stderr, OK); \
END

#define PRINT_FAIL() BEGIN \
	fprintf(stderr, FAIL); \
END

static int passed_tests = 0;
static int passed_assertions = 0;
static int failed_assertions = 0;
static int test_passed;

static int exit_after_fail = 1;

#define STR_BUF_LEN 256
static char str_buf[STR_BUF_LEN];

#define TEST_SET_MSG(string) BEGIN \
	snprintf(strbuf, STR_BUF_LEN, message); \
	END

#define TEST_SET_MSG_F(string, ...) BEGIN \
	END

#define TEST_INIT(argc, argv) BEGIN \
		if (argc >= 1 && strcmp(argv[1], "exit_after_fail") == 0) { \
			exit_after_fail = 1; \
		} else { \
			exit_after_fail = 0; \
		} \
	END

#define TEST_ASSERT(condition, message) BEGIN \
		if (condition) { \
			fprintf(stderr, "%s:%d:%s, %s (%s) ", __FILE__, __LINE__, \
				__FUNCTION__, message, #condition); \
			PRINT_OK(); \
			passed_assertions++; \
		} else { \
			fprintf(stderr, "%s:%d:%s, %s (%s) ", __FILE__, __LINE__, \
				__FUNCTION__, message, #condition); \
			PRINT_FAIL(); \
			test_passed = 1; \
			failed_assertions++; \
			if (exit_after_fail) { \
				exit(EXIT_FAILURE); \
			} \
		} \
	END

#define TEST_RUN(test) BEGIN \
	fprintf(stderr, "TEST_RUN(" COLOR_YELLOW "%s" COLOR_RESET ")\n", \
		#test); \
	test_passed = 0; \
	test(); \
	fprintf(stderr, "test: %s ", #test); \
	if (!test_passed) { \
		passed_tests++; \
		PRINT_OK(); \
	} else {\
		PRINT_FAIL(); \
	} \
END

#define TEST_RESULT(name) BEGIN \
	fprintf(stderr, "file=%s; passed_tests=%d; passed_assertions=%d; " \
		"failed_assertions=%d;\n", name, passed_tests, passed_assertions, \
		failed_assertions); \
END

#endif /* TEST_UTILS_H_ */
