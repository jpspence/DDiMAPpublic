#include "../include/cu/cu.h"

// Declarations of tests
TEST(test1);
TEST(test2);
TEST(test3);

// Collect tests into test suite
TEST_SUITE(T1) {
	TEST_ADD(test1),
	TEST_ADD(test2),
	TEST_ADD(test3),
	TEST_SUITE_CLOSURE
};
