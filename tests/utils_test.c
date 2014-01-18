/*
 * utils_test.c
 *
 *  Created on: Nov 30, 2013
 *      Author: n
 */

#include <stdio.h>
#include <stdlib.h>

#include "test_framework.h"

#include "../src/utils.h"

static void power_of_two_test(void) {
	snprintf(str_buf, STR_BUF_LEN, "is power of two");
	TEST_ASSERT(is_power_of_two(2) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(4) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(8) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(16) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(32) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(64) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(128) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(256) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(512) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(1024) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(2048) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(4096) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(8192) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(16384) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(32768) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(65536) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(131072) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(262144) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(524288) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(1048576) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(2097152) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(4194304) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(8388608) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(16777216) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(33554432) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(67108864) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(134217728) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(268435456) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(536870912) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(1073741824) == 1, str_buf);
	TEST_ASSERT(is_power_of_two(2147483648) == 1, str_buf);

	snprintf(str_buf, STR_BUF_LEN, "is NOT power of two");
	TEST_ASSERT(is_power_of_two(3) == 0, str_buf);
	TEST_ASSERT(is_power_of_two(5) == 0, str_buf);
	TEST_ASSERT(is_power_of_two(6) == 0, str_buf);
	TEST_ASSERT(is_power_of_two(7) == 0, str_buf);
	TEST_ASSERT(is_power_of_two(9) == 0, str_buf);
	TEST_ASSERT(is_power_of_two(10) == 0, str_buf);
	TEST_ASSERT(is_power_of_two(11) == 0, str_buf);
	TEST_ASSERT(is_power_of_two(12) == 0, str_buf);
	TEST_ASSERT(is_power_of_two(13) == 0, str_buf);
	TEST_ASSERT(is_power_of_two(14) == 0, str_buf);
	TEST_ASSERT(is_power_of_two(15) == 0, str_buf);
	TEST_ASSERT(is_power_of_two(1023) == 0, str_buf);
	TEST_ASSERT(is_power_of_two(1025) == 0, str_buf);
	TEST_ASSERT(is_power_of_two(536870911) == 0, str_buf);
	TEST_ASSERT(is_power_of_two(536870913) == 0, str_buf);
	TEST_ASSERT(is_power_of_two(2147483647) == 0, str_buf);
	TEST_ASSERT(is_power_of_two(2147483649) == 0, str_buf);
}

int main(int argc, char *argv[]) {
	srand(0);

	printf("\n");
	TEST_RUN(power_of_two_test);
	TEST_RESULT(argv[0]);

	return EXIT_SUCCESS;
}
