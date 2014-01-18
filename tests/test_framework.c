#include <stdio.h>
#include <stdlib.h>

#include "test_framework.h"

static void test01(void) {
	TEST_ASSERT(1 == 1, "one is one");
	TEST_ASSERT(2 == 2, "two is two");
	TEST_ASSERT((1 + 2) == 3, "one plus two is three");
}

static void test02(void) {
	TEST_ASSERT(1 == 2, "one is not two");
}

static void test03(void) {
	char foo[] = "foo";
	char bar[] = "bar";

	TEST_ASSERT(strcmp(foo, foo) == 0, "compare foo and foo");
	TEST_ASSERT(strcmp(foo, bar) == 0, "compare foo and bar");
}

int main(int argc, char *argv[]) {
	TEST_RUN(test01);
	TEST_RUN(test02);
	TEST_RUN(test03);
	TEST_RESULT(argv[0]);

	return EXIT_SUCCESS;
}
