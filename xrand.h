#pragma once
#include <stdint.h>

/* The state must be seeded so that it is not all zero */
static uint64_t s[] = { 100, 200 };

/* xorshift128+ */
uint64_t xorshift(void) {
    uint64_t x = s[0];
    uint64_t const y = s[1];
    s[0] = y;
    x ^= x << 23; // a
    s[1] = x ^ y ^ (x >> 17) ^ (y >> 26); // b, c
    return s[1] + y;
}

void xsrand(uint64_t seed) {
    s[0] = seed;
    s[1] = 2*seed + 1;
}

float xdrand(void){
    double n = (double) xorshift();
    return n / UINT64_MAX;
}

uint64_t xrand(void) {
    uint64_t val = xorshift();
    return val;
}
