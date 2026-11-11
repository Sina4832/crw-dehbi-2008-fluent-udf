/*
 * rng_check.c - statistical test of gauss_rand() in udf/crw_dehbi2008_rev10.c
 *
 * gauss_rand_copy() reproduces gauss_rand() (UDF Section 1) statement by
 * statement; only the types are adapted so that it runs outside ANSYS Fluent:
 * the Tracked_Particle pointer becomes an integer "address", the particle time
 * is passed explicitly, and unsigned long becomes uint32_t. Only the lowest 31
 * bits of each congruential step are used, so the results are identical for
 * 32-bit (MSVC) and 64-bit unsigned long.
 *
 * Build and run:
 *     gcc -std=c99 -O2 -o rng_check tools/rng_check.c -lm
 *     ./rng_check            (time step 1e-6 s, the UDF upper bound)
 *     ./rng_check 3.1e-5     (any other time step in seconds)
 *
 * Tests (for independent standard normal deviates: mean ~ 0, variance ~ 1,
 * correlations ~ 0 within about +/- 2/sqrt(N), and var(sum of M)/M ~ 1):
 *   [1] mean, variance and autocorrelation (lags 1-5) for one address;
 *   [2] correlation between two addresses at the same particle time;
 *   [3] variance of the sum of M successive deviates, divided by M;
 *   [4] test [3] for 2000 pseudo-random 16-byte-aligned addresses (fixed
 *       seed, reproducible), printed as a distribution;
 *   [5] the Langevin update of the UDF (lines 397-399, 427-429) driven by the
 *       generator: time average of eta^2 per address (target 1), with an
 *       ideal reference generator run alongside (see the comment at [5]).
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

/* Reference generator for test [5] only: splitmix64 with Box-Muller. */
static uint64_t splitmix64(uint64_t *state)
{
    uint64_t z = (*state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static double ideal_gauss(uint64_t *state)
{
    double u1 = ((double)(splitmix64(state) >> 11) + 0.5) / 9007199254740992.0;
    double u2 = ((double)(splitmix64(state) >> 11) + 0.5) / 9007199254740992.0;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
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
    {
        /* [5] Langevin update of the UDF, eta <- (eta + sqrt(2 dt/tau_L) xi)
         * / (1 + dt/tau_L), started from the injection deviate (salt 11, UDF
         * line 678) and driven with salt 1, both for the wall-normal
         * component. tau_L =
         * 3.32e-3 s corresponds to tau_L+ = 10 in this case; the path covers
         * 0.088 s, the residence time at the bulk velocity, and eta^2 is
         * averaged over its second half. The ideal generator shows the spread
         * caused by the finite averaging window alone. */
        enum { NA = 2000 };
        static double v[NA], vi[NA];
        const double tau_l = 3.32e-3, t_res = 0.088;
        const double h = dt / tau_l, amp = sqrt(2.0 * h), den = 1.0 + h;
        const long n_step = (long)(t_res / dt + 0.5), n0 = n_step / 2;
        uint64_t x = 0x9E3779B97F4A7C15ULL, ref = 0x1234567ULL;
        long k;
        for (a = 0; a < NA; a++) {
            uint64_t addr;
            double e, ei, s2 = 0.0, s2i = 0.0;
            x ^= x << 13; x ^= x >> 7; x ^= x << 17;
            addr = 0x7f0000000000ULL + ((x >> 20) & 0xFFFFFFFFF0ULL);
            e = gauss_rand_copy(addr, 0.0, 11);
            ei = ideal_gauss(&ref);
            for (k = 1; k <= n_step; k++) {
                e = (e + amp * gauss_rand_copy(addr, (double)k * dt, 1)) / den;
                ei = (ei + amp * ideal_gauss(&ref)) / den;
                if (k > n0) { s2 += e * e; s2i += ei * ei; }
            }
            v[a] = s2 / (double)(n_step - n0);
            vi[a] = s2i / (double)(n_step - n0);
        }
        qsort(v, NA, sizeof v[0], cmp_double);
        qsort(vi, NA, sizeof vi[0], cmp_double);
        printf("\n[5] time-averaged eta^2 of the UDF Langevin update, "
               "%d addresses (ideal 1)\n", NA);
        printf("    tau_L = %.3g s, %ld steps, average over the second half\n",
               tau_l, n_step);
        printf("    UDF generator  : min %.2f, 5%% %.2f, median %.2f, "
               "95%% %.2f, max %.1f\n", v[0], v[NA / 20], v[NA / 2],
               v[NA - 1 - NA / 20], v[NA - 1]);
        printf("    ideal generator: min %.2f, 5%% %.2f, median %.2f, "
               "95%% %.2f, max %.1f\n", vi[0], vi[NA / 20], vi[NA / 2],
               vi[NA - 1 - NA / 20], vi[NA - 1]);
    }
    return 0;
}
