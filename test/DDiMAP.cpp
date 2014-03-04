#include "../include/cu/cu.h"
#include <stdio.h>

TEST(test1) {
    assertTrue(1);
    printf("Hello from test1\n");
    assertEquals(1, 1);
}

TEST(test2) {
    assertFalse(0);
    assertNotEquals(1, 2);
    fprintf(stderr, "Hello from test2\n");
}

TEST(test3) {
    assertFalse(1);
}
