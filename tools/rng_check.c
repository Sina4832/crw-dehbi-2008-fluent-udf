/*
 * rng_check.c - statistical check of gauss_rand() from udf/crw_dehbi2008_rev10.c
 *
 * The function body below is an exact copy of the UDF (Section 1), with the
 * Tracked_Particle pointer replaced by an integer "address" and a time so that
 * it runs outside ANSYS Fluent. Only the lowest 31 bits of each congruential
 * step are used, so the results are identical for 32-bit (MSVC) and 64-bit
 * unsigned long.
 *
 * Build and run:
 *     gcc -O2 -o rng_check tools/rng_check.c -lm
 *     ./rng_check            (time step 1e-6 s, the UDF upper bound)
 *     ./rng_check 3.1e-5     (any other time step in seconds)
 *
 * For independent standard normal deviates: mean ~ 0, variance ~ 1, all
 * correlations ~ 0 within about +/- 2/sqrt(N), and var(sum of M)/M ~ 1.
 * Test [4] repeats test [3] for 2000 pseudo-random 16-byte-aligned addresses
 * (fixed seed, so the output is reproducible) and prints the distribution.
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef uint32_t ulong32;

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

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
    const uint64_t addrs[] = {0x7ff6a1b2c3d0ULL, 0x1000ULL, 0x2a7f3e10ULL,
                              0x55d0c0ffee40ULL, 0x7f00deadbe00ULL, 0x12345670ULL};
    const uint64_t gaps[] = {0x8, 0x20, 0x1C0, 0x2000};
    const int n_addr = 6, n_gap = 4, M = 1000, B = 200;
    double dt = (argc > 1) ? atof(argv[1]) : 1.0e-6;
    double sum = 0.0, sumsq = 0.0;
    int a, g, lag, n, b, m;

    printf("time step: %g s\n\n", dt);
    for (n = 0; n < N; n++) {
        double x = gauss_rand_copy(addrs[0], n * dt, 0);
        sum += x; sumsq += x * x;
    }
    printf("[1] address %#llx, N = %d\n", (unsigned long long)addrs[0], N);
    printf("    mean, variance       : %+.4f, %.4f\n", sum / N, sumsq / N - (sum / N) * (sum / N));
    for (lag = 1; lag <= 5; lag++) {
        double c = 0.0;
        for (n = 0; n < N; n++)
            c += gauss_rand_copy(addrs[0], n * dt, 0) * gauss_rand_copy(addrs[0], (n + lag) * dt, 0);
        printf("    autocorrelation lag %d: %+.4f\n", lag, c / N);
    }

    printf("\n[2] correlation between two particles at the same time (N = %d)\n", N);
    for (g = 0; g < n_gap; g++) {
        double c = 0.0;
        for (n = 0; n < N; n++)
            c += gauss_rand_copy(addrs[0], n * dt, 0) * gauss_rand_copy(addrs[0] + gaps[g], n * dt, 0);
        printf("    address gap %#6llx: %+.4f\n", (unsigned long long)gaps[g], c / N);
    }

    printf("\n[3] var(sum of %d successive deviates)/%d, %d blocks (ideal 1)\n", M, M, B);
    for (a = 0; a < n_addr; a++) {
        double s2 = 0.0;
        for (b = 0; b < B; b++) {
            double s = 0.0;
            for (m = 0; m < M; m++) s += gauss_rand_copy(addrs[a], (double)(b * M + m) * dt, 0);
            s2 += s * s;
        }
        printf("    address %#16llx: %.2f\n", (unsigned long long)addrs[a], s2 / B / M);
    }

    {
        enum { NA = 2000, B4 = 50 };
        static double v[NA];
        uint64_t x = 0x9E3779B97F4A7C15ULL; /* xorshift64 state, fixed seed */
        for (a = 0; a < NA; a++) {
            uint64_t addr;
            double s2 = 0.0;
            x ^= x << 13; x ^= x >> 7; x ^= x << 17;
            addr = 0x7f0000000000ULL + ((x >> 20) & 0xFFFFFFFFF0ULL);
            for (b = 0; b < B4; b++) {
                double s = 0.0;
                for (m = 0; m < M; m++) s += gauss_rand_copy(addr, (double)(b * M + m) * dt, 0);
                s2 += s * s;
            }
            v[a] = s2 / B4 / M;
        }
        qsort(v, NA, sizeof v[0], cmp_double);
        printf("\n[4] as [3] for %d random 16-byte-aligned addresses, %d blocks each\n", NA, B4);
        printf("    min %.2f, 5%% %.2f, median %.2f, 95%% %.2f, max %.2f\n",
               v[0], v[NA / 20], v[NA / 2], v[NA - 1 - NA / 20], v[NA - 1]);
    }
    return 0;
}
