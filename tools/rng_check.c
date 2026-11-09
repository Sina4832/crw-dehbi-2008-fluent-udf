/*
 * rng_check.c - statistical check of gauss_rand() from udf/crw_dehbi2008_rev10.c
 *
 * The function body below is copied from the UDF (Section 1) with two
 * substitutions so that it runs outside ANSYS Fluent:
 *   - Tracked_Particle* is replaced by a fixed integer "address" and a time;
 *   - unsigned long is emulated as a 32-bit type, as on Windows/MSVC (LLP64),
 *     which is the compiler used to build the UDF for the reported results.
 *
 * Build and run:
 *     gcc -O2 -o rng_check tools/rng_check.c -lm
 *     ./rng_check            (default time step 1e-6 s, the UDF upper bound)
 *     ./rng_check 3.1e-5     (any other time step in seconds)
 *
 * For independent standard normal deviates the expected results are
 * mean ~ 0, variance ~ 1, and all correlations ~ 0 within about
 * +/- 2/sqrt(N) = +/- 0.0045 for N = 200000.
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef uint32_t ulong32; /* MSVC: unsigned long is 32-bit */

static double gauss_rand_copy(uint64_t address, double time, int salt)
{
    ulong32 base = (ulong32)address
                 ^ (ulong32)(uint64_t)(time * 1.0e9)
                 ^ (ulong32)((ulong32)salt * 2654435769UL);
    ulong32 sA = base ^ 0xA5A5A5A5UL;
    ulong32 sB = base ^ 0x5C5C5C5CUL;
    double u1, u2;

    sA = sA * 1664525UL + 1013904223UL;
    u1 = ((double)(sA & 0x7FFFFFFFUL) + 0.5) / (double)0x80000000UL;
    sB = sB * 22695477UL + 1UL;
    u2 = ((double)(sB & 0x7FFFFFFFUL) + 0.5) / (double)0x80000000UL;
    if (u1 < 1.0e-12) u1 = 1.0e-12;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

int main(int argc, char **argv)
{
    const int N = 200000;
    const uint64_t addr = 0x7ff6a1b2c3d0ULL; /* arbitrary heap-like address */
    double dt = (argc > 1) ? atof(argv[1]) : 1.0e-6;
    double sum = 0.0, sumsq = 0.0, cross01 = 0.0, pp = 0.0;
    int lag, n;

    for (n = 0; n < N; n++) {
        double t = n * dt;
        double x0 = gauss_rand_copy(addr, t, 0);
        double x1 = gauss_rand_copy(addr, t, 1);
        double y0 = gauss_rand_copy(addr + 0x1C0, t, 0); /* neighbouring particle */
        sum += x0; sumsq += x0 * x0; cross01 += x0 * x1; pp += x0 * y0;
    }
    printf("time step            : %g s, N = %d\n", dt, N);
    printf("mean, variance       : %+.4f, %.4f\n", sum / N, sumsq / N - (sum / N) * (sum / N));
    for (lag = 1; lag <= 5; lag++) {
        double a = 0.0;
        for (n = 0; n < N; n++)
            a += gauss_rand_copy(addr, n * dt, 0) * gauss_rand_copy(addr, (n + lag) * dt, 0);
        printf("autocorrelation lag %d: %+.4f\n", lag, a / N);
    }
    printf("components 0-1 corr  : %+.4f\n", cross01 / N);
    printf("two-particle corr    : %+.4f  (addresses 0x1C0 bytes apart, same time)\n", pp / N);
    return 0;
}
