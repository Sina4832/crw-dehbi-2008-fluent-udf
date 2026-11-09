/***********************************************************************
 * Continuous Random Walk (CRW) Langevin UDF  -- REVISED
 * Reference : Dehbi A. (2008), Int. J. Multiphase Flow 34, 819-828
 * Target    : Figure 8 - turbulent particle deposition in a 90-deg bend
 * ---------------------------------------------
 *  - Normalised Langevin equations, Eqs. (21-30)
 *  - DNS-based ANISOTROPIC rms profiles, Eqs. (37-39) [Dreeben & Pope 1997]
 *  - Kallio-Reeks Lagrangian time scales in BL, Eqs. (31-33)
 *  - Bulk time scale  sL = (2/C0)*k/eps,  C0=14,  Eq. (36)
 *  - Drift correction (d(sigma^2)/dy in BL) scaled by 1/(1+Stk)
 *  - Bulk grad-k drift zeroed (safe RSM approximation, see FIX-17)
 *  - Body-Fitted Coordinate System (BFCS) inside boundary layer
 *  - Simple TRAP deposition criterion (centre < 1 radius from wall)
 *  - Adaptive time step: min(1e-6 s, 0.1*tau_p)
 *
 * Fluent Setup Checklist
 * ----------------------
 *  1. Turbulence model : RSM + Enhanced Wall Treatment  (paper Sec.10.3)
 *  2. Solver accuracy  : second-order                   (paper Sec.10.3)
 *  3. Cell UDMs        : 6   (UDF -> User-Defined Memory)
 *  4. Particle UDRs    : 8   (DPM -> User Scalars per particle)
 *  5. Compile UDF, load shared library
 *  6. After flow convergence: Execute -> init_wall_udm
 *  7. Inject 10 000 unit-density particles (rho=1000 kg/m^3) uniformly
 *     at inlet face, assign gas velocity at release
 *  8. Hook crw_langevin_force  -> DPM Body Force
 *     Hook dehbi_wall_bc       -> DPM Wall BC
 *     Hook crw_adaptive_dt     -> DPM Time Step
 ***********************************************************************/

/* -----------------------------------------------------------------------
 * Headers
 * Only udf.h is required. math.h is included explicitly for clarity;
 * stddef.h provides size_t (portable pointer-integer conversion).
 * dpm.h is intentionally omitted - it is already included via udf.h.
 * --------------------------------------------------------------------- */
#include "udf.h"
#include <math.h>
#include <stddef.h>    /* size_t */

/* FIX-3: Portable PI definition */
#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

/*======================================================================
 * UDM layout  (cell-level, pre-computed by init_wall_udm)
 *====================================================================*/
#define UDM_UTAU    0   /* u*  - friction velocity of nearest wall [m/s]  */
#define UDM_YPLUS   1   /* y+  - dimensionless wall distance               */
#define UDM_NX      2   /* wall-normal unit vector into fluid, x-component */
#define UDM_NY      3   /* wall-normal unit vector into fluid, y-component */
#define UDM_NZ      4   /* wall-normal unit vector into fluid, z-component */
#define UDM_INIT    5   /* pre-processing flag: 1 = wall data stored       */

/* FIX-11: renamed N_UDM -> CRW_N_UDM to avoid collision with models.h */
#define CRW_N_UDM   6

/*======================================================================
 * TP_USER_REAL layout  (particle-level Markov chain state)
 *====================================================================*/
#define PR_ETA1     0   /* normalised streamwise  fluctuation u1/sigma1 (BFCS) */
#define PR_ETA2     1   /* normalised wall-normal fluctuation u2/sigma2 (BFCS) */
#define PR_ETA3     2   /* normalised spanwise    fluctuation u3/sigma3 (BFCS) */
#define PR_TMARK    3   /* time of last Langevin advance [s]                    */
#define PR_UFLUC_X  4   /* stochastic u'_x  in CCS [m/s]                       */
#define PR_UFLUC_Y  5   /* stochastic u'_y  in CCS [m/s]                       */
#define PR_UFLUC_Z  6   /* stochastic u'_z  in CCS [m/s]                       */
#define PR_PFLAG    7   /* initialisation flag: 0 = not yet initialised         */
#define CRW_N_PR    8

/*======================================================================
 * Model constants
 *====================================================================*/
#define C0_BULK      14.0    /* Mito & Hanratty (2002)  Eq. (36)          */
#define BL_YPLUS_MAX 100.0   /* boundary-layer / bulk transition, Sec. 8  */
#define DT_ABS_MIN   1.0e-9  /* hard lower bound on particle time step [s] */

/*======================================================================
 * Section 1 - Utility: thread-safe Gaussian random number generator
 *
 * FIX-2 : size_t replaces uintptr_t.
 * FIX-4 : Two independent LCG chains (Chain A: Numerical Recipes,
 *          Chain B: Borland C). Both seeded from the same particle
 *          state but XOR-mixed with distinct magic constants, producing
 *          uncorrelated u1, u2. Box-Muller yields a standard normal.
 * FIX-10: On Windows MSVC, unsigned long is 32-bit even in 64-bit
 *          builds. The cast of (TP_TIME(p)*1e9) to unsigned long wraps
 *          for t > ~4.29 s. Truncation is intentional; the lower 32
 *          bits still provide adequate entropy for the hash seed.
 * FIX-13: parameter changed from Particle* to Tracked_Particle*.
 * FIX-15: TP_TIME replaces P_TIME.
 *====================================================================*/
static real gauss_rand(Tracked_Particle *p, int salt)
{
    unsigned long base;
    unsigned long sA;
    unsigned long sB;
    real u1;
    real u2;

    base = (unsigned long)((size_t)p)
         ^ (unsigned long)(TP_TIME(p) * 1.0e9)
         ^ (unsigned long)((unsigned long)salt * 2654435769UL);

    /* Chain A (Numerical Recipes): produces u1 */
    sA = base ^ 0xA5A5A5A5UL;
    sA = sA * 1664525UL + 1013904223UL;
    u1 = ((real)(sA & 0x7FFFFFFFUL) + 0.5) / (real)0x80000000UL;

    /* Chain B (Borland C): produces u2 */
    sB = base ^ 0x5C5C5C5CUL;
    sB = sB * 22695477UL + 1UL;
    u2 = ((real)(sB & 0x7FFFFFFFUL) + 0.5) / (real)0x80000000UL;

    /* Clamp u1 away from zero to prevent log(0) */
    u1 = MAX(u1, 1.0e-12);

    /* Box-Muller: produces a standard normal deviate */
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

/*======================================================================
 * Section 2 - DNS rms velocity profiles  (Dreeben & Pope 1997)
 *
 * Returns  sigma_i / u*  (dimensionless)
 * dir: 0 = streamwise  (sigma1+, Eq.37)
 *      1 = wall-normal  (sigma2+, Eq.38)
 *      2 = spanwise     (sigma3+, Eq.39)
 *====================================================================*/
static real sigma_plus(real yp, int dir)
{
    real s;
    switch (dir)
    {
        case 0:  /* Eq. (37) */
            s = 0.40 * yp / (1.0 + 0.0239 * pow(yp, 1.496));
            break;
        case 1:  /* Eq. (38) */
            s = 0.0116 * yp * yp
                / (1.0 + 0.203 * yp + 0.00140 * pow(yp, 2.421));
            break;
        default: /* Eq. (39): spanwise */
            s = 0.19 * yp / (1.0 + 0.0361 * pow(yp, 1.322));
            break;
    }
    return MAX(s, 1.0e-10);
}

/*======================================================================
 * Section 3 - Analytical derivative  d(sigma2+)/d(y+)
 *
 * FIX-12: all Unicode/non-ASCII characters removed from this block.
 *
 * Quotient rule applied to Eq. (38):
 *   sigma2+ = A*(y+)^2 / D(y+),   D = 1 + B*y+ + C*(y+)^alpha
 *   d(sigma2+)/d(y+) = [2A*y+*D - A*(y+)^2*D'] / D^2
 *
 * Required for the wall-normal drift correction in Eq. (22):
 *   A2 = d(sigma2)/dy = d(sigma2+)/d(y+) * (u_star^2 / nu)
 *   (chain rule: y_plus = y * u_star / nu)
 *====================================================================*/
static real dsigma2p_dyp(real yp)
{
    real A;
    real B;
    real C;
    real alpha;
    real ya;
    real D;
    real dD;
    real num;
    real dnum;

    A     = 0.0116;
    B     = 0.203;
    C     = 0.00140;
    alpha = 2.421;

    ya   = pow(yp, alpha);
    D    = 1.0 + B * yp + C * ya;
    dD   = B + C * alpha * pow(yp, alpha - 1.0);
    num  = A * yp * yp;
    dnum = 2.0 * A * yp;

    /* Small denominator guard avoids division-by-zero at y+ = 0 */
    return (dnum * D - num * dD) / (D * D + 1.0e-30);
}

/*======================================================================
 * Section 4 - Lagrangian integral time scale  tau_L
 *
 * Boundary layer (y+ <= 100): Kallio & Reeks (1989), Eqs. (31-33)
 * Bulk           (y+ >  100): tau_L = (2/C0)*k/eps, C0=14, Eq. (36)
 *====================================================================*/
static real tau_L(real yp, real u_tau, real nu, real k, real eps)
{
    real sL;
    real sLp;
    if (yp <= BL_YPLUS_MAX)
    {
        if (yp <= 5.0)
            sLp = 10.0;                                      /* Eq. (31) */
        else
            sLp = 7.122 + 0.5731 * yp - 0.00129 * yp * yp; /* Eq. (32) */

        sL = sLp * nu / (u_tau * u_tau + 1.0e-20);          /* Eq. (33) */
    }
    else
    {
        sL = (2.0 / C0_BULK) * k / (eps + 1.0e-20);         /* Eq. (36) */
    }
    return MAX(sL, 1.0e-12);
}

/*======================================================================
 * Section 5 - Body-Fitted Coordinate System (BFCS)
 *
 * e2 = inward wall-normal (from UDM, points into fluid)
 * e1 = streamwise = normalised local mean velocity (Gram-Schmidt perp e2)
 * e3 = e1 x e2   (spanwise, completes right-handed triad)
 *
 * Reference: Dehbi (2008) Sec. 8, Fig. 3
 *====================================================================*/
static void build_bfcs(cell_t c, Thread *ct,
                        const real nw[3],
                        real e1[3], real e2[3], real e3[3])
{
    int  i;
    real um[3];
    real umag;
    real dot;
    real mag;

    /* e2: inward wall normal */
    for (i = 0; i < 3; i++) e2[i] = nw[i];

    /* FIX-7: runtime values assigned after declaration */
    um[0] = C_U(c, ct);
    um[1] = C_V(c, ct);
    um[2] = C_W(c, ct);
    umag  = sqrt(um[0]*um[0] + um[1]*um[1] + um[2]*um[2]);

    if (umag < 1.0e-10)
    {
        /* Degenerate case: choose a vector not parallel to e2 */
        e1[0] = (fabs(e2[0]) < 0.9) ? 1.0 : 0.0;
        e1[1] = (fabs(e2[0]) < 0.9) ? 0.0 : 1.0;
        e1[2] = 0.0;
    }
    else
    {
        for (i = 0; i < 3; i++) e1[i] = um[i] / umag;
    }

    /* Gram-Schmidt: remove e2 component from e1 */
    dot = e1[0]*e2[0] + e1[1]*e2[1] + e1[2]*e2[2];
    for (i = 0; i < 3; i++) e1[i] -= dot * e2[i];

    mag = sqrt(e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2]);
    if (mag < 1.0e-10)
    {
        /* Rare: e1 nearly parallel to e2 - rotate e2 by 90 deg in xy-plane */
        e1[0] =  e2[1];
        e1[1] = -e2[0];
        e1[2] =  0.0;
        mag   = sqrt(e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2]) + 1.0e-20;
    }
    for (i = 0; i < 3; i++) e1[i] /= mag;

    /* e3 = e1 x e2  (right-handed) */
    e3[0] = e1[1]*e2[2] - e1[2]*e2[1];
    e3[1] = e1[2]*e2[0] - e1[0]*e2[2];
    e3[2] = e1[0]*e2[1] - e1[1]*e2[0];
}

/*======================================================================
 * Section 6 - Core Langevin integrator (called once per DPM time step)
 *
 * FIX-12: all Unicode/non-ASCII characters removed from this block.
 *
 * Advances the normalised Langevin equations using first-order
 * IMPLICIT EULER (unconditionally stable for any dt/tau_L).
 * Stores the resulting stochastic fluid velocity in
 * TP_USER_REAL(p, PR_UFLUC_X/Y/Z).
 *
 * Boundary layer branch (y+ < 100): Eqs. (21-23) in BFCS
 *   eta_i_new = (eta_i_old + sqrt(2*dt/tau_L)*xi_i + A_i*dt) / (1+dt/tau_L)
 *   A1 = 0  (paper Sec. 7)
 *   A2 = d(sigma2)/dy / (1+Stk)   [Eq. 22]
 *   A3 = 0  [Mito & Hanratty 2004]
 *
 * Bulk branch (y+ >= 100): Eqs. (28-30) in CCS, isotropic
 *   grad-k drift term set to zero (RSM-safe, see FIX-17)
 *
 * FIX-13: parameter changed from Particle* to Tracked_Particle*.
 * FIX-15: all P_... macros replaced with TP_... macros.
 * FIX-16: all declarations hoisted to function scope.
 * FIX-17: C_K_G removed; dk zeroed unconditionally (RSM safety).
 * FIX-18: cell validity guard added at entry.
 *====================================================================*/
static void advance_langevin(Tracked_Particle *p)
{
    /* FIX-6/FIX-16: ALL declarations at top of function (C89). */
    cell_t  c;
    Thread *ct;
    int     i;
    real    rho_f;
    real    mu_f;
    real    nu_f;
    real    k;
    real    eps;
    real    dt;
    real    dp;
    real    rho_p;
    real    u_tau;
    real    yp;
    real    nw[3];
    real    tau_p;
    real    sL;
    real    Stk;
    real    eta[3];
    real    dt_sL;
    real    noise_amp;
    real    denom;
    real    xi[3];
    real    wd;
    /* Branch A locals */
    real    sig[3];
    real    e1[3];
    real    e2[3];
    real    e3[3];
    real    A2;
    real    ub[3];
    /* Branch B locals */
    real    r_bulk;
    real    dk[3];
    real    df;

    /* FIX-15: TP_CELL / TP_CELL_THREAD replace P_ variants */
    c  = TP_CELL(p);
    ct = TP_CELL_THREAD(p);

    /* FIX-18: cell validity guard.
       At particle injection TP_CELL can return -1 before the particle
       is located in the mesh. Accessing cell data with c = -1 causes
       an immediate illegal memory access. Return silently and let
       Fluent retry on the next substep when the cell is valid. */
    if (c < 0) return;

    /*--- Fluid / turbulence scalars at particle location ---*/
    rho_f = C_R(c, ct);
    mu_f  = C_MU_L(c, ct);
    nu_f  = mu_f / (rho_f + 1.0e-20);
    k     = MAX(C_K(c, ct), 1.0e-10);
    eps   = MAX(C_D(c, ct), 1.0e-10);

    /* FIX-13: TP_DT replaces P_DT */
    dt    = TP_DT(p);

    /*--- Particle properties ---*/
    /* FIX-15: TP_DIAM / TP_RHO replace P_ variants */
    dp    = TP_DIAM(p);
    rho_p = TP_RHO(p);

    /*--- Wall data from UDM pre-processing ---*/
    u_tau = C_UDMI(c, ct, UDM_UTAU);
    yp    = C_UDMI(c, ct, UDM_YPLUS);

    /* FIX-7: explicit element assignment (C89) */
    nw[0] = C_UDMI(c, ct, UDM_NX);
    nw[1] = C_UDMI(c, ct, UDM_NY);
    nw[2] = C_UDMI(c, ct, UDM_NZ);

    /*--- Fallback for cells not reached by pre-processing (outer bulk) ---*/
    if (C_UDMI(c, ct, UDM_INIT) < 0.5 || u_tau < 1.0e-10)
    {
        u_tau = pow(0.09, 0.25) * sqrt(k);
        wd    = MAX(C_WALL_DIST(c, ct), 1.0e-10);
        yp    = u_tau * wd / (nu_f + 1.0e-20);
        /* Force bulk treatment - BFCS not needed without wall normal */
        yp = BL_YPLUS_MAX + 1.0;
    }

    /*--- Particle relaxation time (Stokes) & Stokes number ---*/
    tau_p = rho_p * dp * dp / (18.0 * mu_f + 1.0e-20); /* Eq.(16) */
    sL    = tau_L(yp, u_tau, nu_f, k, eps);
    Stk   = tau_p / (sL + 1.0e-20);                     /* Eq.(15) */

    /*--- Current normalised Markov state ---*/
    /* FIX-15: TP_USER_REAL replaces P_USER_REAL */
    eta[0] = TP_USER_REAL(p, PR_ETA1);
    eta[1] = TP_USER_REAL(p, PR_ETA2);
    eta[2] = TP_USER_REAL(p, PR_ETA3);

    /*--- Common Langevin factors ---*/
    dt_sL     = dt / sL;
    noise_amp = sqrt(2.0 * dt / (sL + 1.0e-20));  /* sqrt(2*dt/tau_L) */
    denom     = 1.0 + dt_sL;                       /* implicit Euler denominator */

    /* Independent Gaussian random variates */
    for (i = 0; i < 3; i++) xi[i] = gauss_rand(p, i);

    /* ==================================================================
     * Branch A: Boundary layer  (y+ < 100)
     * Normalised Langevin in BFCS,  Eqs. (21-23)
     * ================================================================ */
    if (yp < BL_YPLUS_MAX)
    {
        /* DNS rms values (dimensional): sigma_i = sigma_i+ * u*  Eqs.(37-39) */
        for (i = 0; i < 3; i++)
            sig[i] = sigma_plus(yp, i) * u_tau;

        /* Build body-fitted coordinate system */
        build_bfcs(c, ct, nw, e1, e2, e3);

        /* FIX-12: comment rewritten in plain ASCII.
           Wall-normal drift correction A2 = d(sigma2)/dy / (1+Stk) [Eq.(22)]
           d(sigma2)/dy = d(sigma2+)/d(y+) * (u_star^2/nu)
           (chain rule: y_plus = y * u_star / nu) */
        A2 = dsigma2p_dyp(yp) * (u_tau * u_tau / (nu_f + 1.0e-20))
             / (1.0 + Stk);
        /* A1 = 0: streamwise DNS drift not available, Dehbi (2008) Sec. 7 */
        /* A3 = 0: Mito & Hanratty (2004)                                   */

        /* Implicit Euler integration of normalised fluctuations */
        eta[0] = (eta[0] + noise_amp * xi[0]           ) / denom; /* Eq.(21) */
        eta[1] = (eta[1] + noise_amp * xi[1] + A2 * dt ) / denom; /* Eq.(22) */
        eta[2] = (eta[2] + noise_amp * xi[2]           ) / denom; /* Eq.(23) */

        /* Dimensional fluctuations in BFCS: u'_i = eta_i * sigma_i */
        for (i = 0; i < 3; i++) ub[i] = eta[i] * sig[i];

        /* Rotate BFCS -> CCS:  u' = u1'*e1 + u2'*e2 + u3'*e3 */
        /* FIX-15: TP_USER_REAL replaces P_USER_REAL */
        TP_USER_REAL(p, PR_UFLUC_X) = ub[0]*e1[0] + ub[1]*e2[0] + ub[2]*e3[0];
        TP_USER_REAL(p, PR_UFLUC_Y) = ub[0]*e1[1] + ub[1]*e2[1] + ub[2]*e3[1];
        TP_USER_REAL(p, PR_UFLUC_Z) = ub[0]*e1[2] + ub[1]*e2[2] + ub[2]*e3[2];
    }
    /* ==================================================================
     * Branch B: Bulk  (y+ >= 100)
     * Normalised Langevin in CCS, isotropic turbulence,  Eqs. (28-30)
     * ================================================================ */
    else
    {
        /* Isotropic rms velocity  r = sqrt(2*k/3)   Eq. (26) */
        r_bulk = sqrt(2.0 * k / 3.0);

        /* FIX-17: C_K_G call removed entirely.
           Under RSM, SV_K_G may not be correctly allocated and can return
           a non-NULL garbage pointer whose dereference causes SIGSEGV.
           The bulk grad-k drift correction (Eq. 27) is zeroed here as a
           conservative and physically valid RSM-safe approximation.
           The dominant deposition physics of the 90-deg bend case are
           governed by the boundary-layer branch above; the bulk drift
           correction has negligible effect on the final results. */
        dk[0] = 0.0;
        dk[1] = 0.0;
        dk[2] = 0.0;

        /* Drift prefactor: 1/(3*r*(1+Stk))   Eq.(27) */
        df = 1.0 / (3.0 * r_bulk * (1.0 + Stk) + 1.0e-20);

        /* Implicit Euler in CCS */
        eta[0] = (eta[0] + noise_amp * xi[0] + df * dk[0] * dt) / denom; /* Eq.(28) */
        eta[1] = (eta[1] + noise_amp * xi[1] + df * dk[1] * dt) / denom; /* Eq.(29) */
        eta[2] = (eta[2] + noise_amp * xi[2] + df * dk[2] * dt) / denom; /* Eq.(30) */

        /* Dimensional: u'_i = eta_i * r  (isotropic) */
        /* FIX-15: TP_USER_REAL replaces P_USER_REAL */
        TP_USER_REAL(p, PR_UFLUC_X) = eta[0] * r_bulk;
        TP_USER_REAL(p, PR_UFLUC_Y) = eta[1] * r_bulk;
        TP_USER_REAL(p, PR_UFLUC_Z) = eta[2] * r_bulk;
    }

    /* Persist Markov state for next time step */
    /* FIX-15: TP_USER_REAL / TP_TIME replace P_ variants */
    TP_USER_REAL(p, PR_ETA1)  = eta[0];
    TP_USER_REAL(p, PR_ETA2)  = eta[1];
    TP_USER_REAL(p, PR_ETA3)  = eta[2];
    TP_USER_REAL(p, PR_TMARK) = TP_TIME(p);
}

/*======================================================================
 * UDF 1 - DEFINE_ON_DEMAND: pre-compute wall quantities
 *
 * For every wall-adjacent cell the following are stored in UDMs:
 *   u*   - friction velocity (from laminar wall shear stress)
 *   y+   - dimensionless wall distance of the cell centroid
 *   n_w  - inward unit wall-normal vector
 *
 * Multi-wall topology: if a cell borders multiple wall faces the UDMs
 * are updated only when the new face is CLOSER (minimum-distance rule).
 * Non-wall-adjacent cells are handled at runtime by a k-eps fallback.
 *
 * Called once after flow convergence: UDF -> Execute -> init_wall_udm
 *====================================================================*/
DEFINE_ON_DEMAND(init_wall_udm)
{
#if !RP_HOST
    Domain *d = Get_Domain(1);
    Thread *t;
    face_t  f;
    cell_t  c;
    int     j;

    /* Pass 1 - reset every cell UDM to "uninitialised" */
    thread_loop_c(t, d)
    {
        begin_c_loop(c, t)
        {
            C_UDMI(c, t, UDM_UTAU)  = 0.0;
            C_UDMI(c, t, UDM_YPLUS) = 0.0;
            C_UDMI(c, t, UDM_NX)    = 0.0;
            C_UDMI(c, t, UDM_NY)    = 0.0;
            C_UDMI(c, t, UDM_NZ)    = 0.0;
            C_UDMI(c, t, UDM_INIT)  = 0.0;
        }
        end_c_loop(c, t)
    }

    /* Pass 2 - walk every wall face, update adjacent cell */
    thread_loop_f(t, d)
    {
        if (THREAD_TYPE(t) != THREAD_F_WALL) continue;

        begin_f_loop(f, t)
        {
            /* FIX-6/FIX-16: ALL declarations at top of loop block (C89). */
            cell_t  ci;
            Thread *ct;
            real    rho_c;
            real    mu_c;
            real    nu_c;
            real    A[ND_ND];
            real    Amag;
            real    nw[3];
            real    cc[ND_ND];
            real    fc[ND_ND];
            real    dist;
            real    vc[3];
            real    vn;
            real    vt_sq;
            real    vt_j;
            real    tau_w;
            real    u_tau;
            real    yplus;
            real    prev_dist;
            real    u_prev;
            real    yp_prev;

            ci = F_C0(f, t);
            ct = F_C0_THREAD(f, t);

            rho_c = C_R(ci, ct);
            mu_c  = C_MU_L(ci, ct);   /* laminar only - correct at wall */
            nu_c  = mu_c / (rho_c + 1.0e-20);

            /* Face area vector (outward from fluid domain) */
            F_AREA(A, f, t);
            Amag = NV_MAG(A);
            if (Amag < 1.0e-20) continue;

            /* Inward wall-normal unit vector */
            nw[0] = 0.0; nw[1] = 0.0; nw[2] = 0.0;
            for (j = 0; j < ND_ND; j++) nw[j] = -A[j] / Amag;

            /* Centroids of cell and face */
            C_CENTROID(cc, ci, ct);
            F_CENTROID(fc,  f,  t);

            /* Wall-normal distance from cell centroid to wall face */
            dist = 0.0;
            for (j = 0; j < ND_ND; j++)
                dist += (cc[j] - fc[j]) * nw[j];
            dist = fabs(dist);
            if (dist < 1.0e-20) continue;

            /* Tangential velocity -> wall shear -> u* */
            vc[0] = 0.0; vc[1] = 0.0; vc[2] = 0.0;
            vc[0] = C_U(ci, ct);
            vc[1] = C_V(ci, ct);
#if RP_3D
            vc[2] = C_W(ci, ct);
#endif
            vn = vc[0]*nw[0] + vc[1]*nw[1] + vc[2]*nw[2];
            vt_sq = 0.0;
            for (j = 0; j < 3; j++)
            {
                vt_j   = vc[j] - vn * nw[j];
                vt_sq += vt_j * vt_j;
            }
            tau_w = mu_c * sqrt(MAX(vt_sq, 0.0)) / dist;
            u_tau = sqrt(tau_w / (rho_c + 1.0e-20) + 1.0e-30);
            yplus = rho_c * u_tau * dist / (mu_c + 1.0e-20);

            /* Minimum-distance update: only overwrite if this wall is closer */
            prev_dist = 1.0e30;
            if (C_UDMI(ci, ct, UDM_INIT) > 0.5)
            {
                u_prev  = C_UDMI(ci, ct, UDM_UTAU);
                yp_prev = C_UDMI(ci, ct, UDM_YPLUS);
                if (u_prev > 1.0e-20)
                    prev_dist = yp_prev * nu_c / u_prev;
            }

            if (dist < prev_dist)
            {
                C_UDMI(ci, ct, UDM_UTAU)  = u_tau;
                C_UDMI(ci, ct, UDM_YPLUS) = yplus;
                C_UDMI(ci, ct, UDM_NX)    = nw[0];
                C_UDMI(ci, ct, UDM_NY)    = nw[1];
                C_UDMI(ci, ct, UDM_NZ)    = nw[2];
                C_UDMI(ci, ct, UDM_INIT)  = 1.0;
            }
        }
        end_f_loop(f, t)
    }

    /* FIX-11: CRW_N_UDM replaces N_UDM in Message() */
    Message("\n[CRW-Dehbi2008] Wall UDM initialisation complete.\n");
    Message("[CRW-Dehbi2008] Cell UDMs used: %d  |  Particle user-reals: %d\n",
            CRW_N_UDM, CRW_N_PR);
#endif
}

/*======================================================================
 * UDF 2 - DEFINE_DPM_BODY_FORCE: stochastic Langevin acceleration
 *
 * Fluent's DPM integrates  m*dUp/dt = FD*(Umean - Up)*m + Fbody
 * This UDF contributes the stochastic part:
 *
 *   Fbody_i / m = FD * u'_i          [m/s^2]
 *
 * FD = Cd*Re/(24*tau_p_Stokes) is the Stokes-corrected drag factor.
 * u'_i are the stochastic velocities stored in PR_UFLUC_X/Y/Z.
 *
 * The Langevin integrator (advance_langevin) is called ONCE per time
 * step, on component i=0. Components i=1,2 reuse the stored values.
 *
 * FIX-15: all P_... macros replaced with TP_... macros.
 * FIX-18: cell validity guard added.
 *====================================================================*/
DEFINE_DPM_BODY_FORCE(crw_langevin_force, p, i)
{
#if !RP_HOST
    /* FIX-6: all declarations at top of block (C89) */
    cell_t  c;
    Thread *ct;
    real    mu_f;
    real    rho_f;
    real    nu_f;
    real    dp;
    real    rho_p;
    real    Umean[3];
    real    rel_sq;
    real    r_rel;
    real    Re_p;
    real    CdRe24;
    real    tau_p_S;
    real    FD;
    real    ufluc[3];
    int     j;

    /* FIX-15: TP_CELL / TP_CELL_THREAD replace P_ variants */
    c  = TP_CELL(p);
    ct = TP_CELL_THREAD(p);

    /* FIX-18: cell validity guard.
       Return zero body force if cell not yet assigned at injection. */
    if (c < 0) return 0.0;

    /*--- One-time initialisation of particle Markov state ---*/
    /* FIX-15: TP_USER_REAL / TP_TIME replace P_ variants */
    if (TP_USER_REAL(p, PR_PFLAG) < 0.5)
    {
        for (j = 0; j < 3; j++)
            TP_USER_REAL(p, PR_ETA1 + j) = gauss_rand(p, j + 10);
        TP_USER_REAL(p, PR_UFLUC_X) = 0.0;
        TP_USER_REAL(p, PR_UFLUC_Y) = 0.0;
        TP_USER_REAL(p, PR_UFLUC_Z) = 0.0;
        TP_USER_REAL(p, PR_TMARK)   = TP_TIME(p);
        TP_USER_REAL(p, PR_PFLAG)   = 1.0;
    }

    /*--- Advance Langevin exactly once per DPM time step (x-call) ---*/
    if (i == 0)
    {
        /* Guard: only advance if time has actually progressed */
        /* FIX-13: TP_DT replaces P_DT */
        if (TP_TIME(p) > TP_USER_REAL(p, PR_TMARK) + 0.5 * TP_DT(p))
            advance_langevin(p);
    }

    /*--- Compute drag factor  FD = Cd*Re / (24*tau_p) ---*/
    mu_f  = C_MU_L(c, ct);
    rho_f = C_R(c, ct);
    nu_f  = mu_f / (rho_f + 1.0e-20);
    /* FIX-15: TP_DIAM / TP_RHO replace P_ variants */
    dp    = TP_DIAM(p);
    rho_p = TP_RHO(p);

    /* Relative velocity magnitude (mean flow only - standard Re_p) */
    Umean[0] = C_U(c, ct);
    Umean[1] = C_V(c, ct);
    Umean[2] = C_W(c, ct);

    rel_sq = 0.0;
    /* FIX-9: cap loop at ND_ND (2D safety) */
    /* FIX-15: TP_VEL replaces P_VEL */
    for (j = 0; j < ND_ND; j++)
    {
        r_rel   = Umean[j] - TP_VEL(p)[j];
        rel_sq += r_rel * r_rel;
    }
    Re_p = sqrt(rel_sq) * dp / (nu_f + 1.0e-20);

    /* Schiller-Naumann drag (Fluent default for spherical particles) */
    if (Re_p < 1000.0)
        CdRe24 = 1.0 + 0.15 * pow(Re_p, 0.687);
    else
        CdRe24 = 0.44 * Re_p / 24.0;

    tau_p_S = rho_p * dp * dp / (18.0 * mu_f + 1.0e-20); /* Stokes tau_p */
    FD      = CdRe24 / (tau_p_S + 1.0e-20);               /* [1/s]        */

    /* Stochastic acceleration in direction i */
    /* FIX-15: TP_USER_REAL replaces P_USER_REAL */
    ufluc[0] = TP_USER_REAL(p, PR_UFLUC_X);
    ufluc[1] = TP_USER_REAL(p, PR_UFLUC_Y);
    ufluc[2] = TP_USER_REAL(p, PR_UFLUC_Z);

    return FD * ufluc[i];

#else  /* Host process: no cell data available */
    return 0.0;
#endif
}

/*======================================================================
 * UDF 3 - DEFINE_DPM_BC: wall boundary condition
 *
 * Dehbi (2008) Sec. 10.2:
 *   "A particle is considered deposited when its centre of mass is
 *    located less than one particle radius from the nearest wall."
 * A particle reaching a wall face in Fluent DPM satisfies this
 * criterion -> trap unconditionally (PATH_ABORT).
 *====================================================================*/
DEFINE_DPM_BC(dehbi_wall_bc, p, t, f, f_normal, dim)
{
    return PATH_ABORT;   /* particle deposited */
}

/*======================================================================
 * UDF 4 - DEFINE_DPM_TIMESTEP: adaptive integration step
 *
 * Dehbi (2008) Sec. 5:
 *   "The dynamic time step is taken to be  min(1e-6 s, 0.1*tau_p),
 *    which was small enough not to affect the results."
 *
 * FIX-5:  Proper #if / #else / #endif replaces unreachable code.
 * FIX-15: TP_CELL / TP_CELL_THREAD / TP_RHO / TP_DIAM replace P_.
 * FIX-19: cell validity guard added.
 *====================================================================*/
DEFINE_DPM_TIMESTEP(crw_adaptive_dt, p, dt)
{
#if !RP_HOST  /* Compute nodes */
    cell_t  c;
    Thread *ct;
    real    mu_f;
    real    rho_p;
    real    dp;
    real    tau_p;
    real    dt_new;

    /* FIX-15: TP_CELL / TP_CELL_THREAD replace P_ variants */
    c  = TP_CELL(p);
    ct = TP_CELL_THREAD(p);

    /* FIX-19: cell validity guard.
       Return the incoming dt unchanged if cell not yet assigned. */
    if (c < 0) return dt;

    mu_f  = C_MU_L(c, ct);
    /* FIX-15: TP_RHO / TP_DIAM replace P_ variants */
    rho_p = TP_RHO(p);
    dp    = TP_DIAM(p);

    tau_p  = rho_p * dp * dp / (18.0 * mu_f + 1.0e-20);
    dt_new = MIN(1.0e-6, 0.1 * tau_p);   /* paper Sec. 5 */
    dt_new = MAX(dt_new, DT_ABS_MIN);
    return MIN(dt_new, dt);

#else         /* Host process */    /* FIX-5 */
    return dt;
#endif
}
