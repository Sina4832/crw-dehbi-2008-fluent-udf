---
title: "Reproduction of the Dehbi (2008) CRW–Langevin Particle Dispersion Model in a Turbulent 90° Pipe Bend Using ANSYS Fluent and a User-Defined Function"
subtitle: "Technical report, revision 2"
author: "Sina Karbakhsh Ravari (M.Sc. candidate, Aerospace Engineering, Shiraz University). Supervisor: Dr. Ebrahim Ghahramani"
date: "November 2026"
---

## Revision note

Revision 1 of this report (March 2026) is kept unchanged in
[`original/report-dehbi-v1-2026-03.pdf`](original/report-dehbi-v1-2026-03.pdf).
Revision 2 keeps the same simulations, settings and results. It changes the
presentation in four ways:

1. The model equations are written out as they are implemented in
   [`udf/crw_dehbi2008_rev10.c`](../udf/crw_dehbi2008_rev10.c), with source
   line references.
2. The deposition fractions are tabulated from the Fluent DPM summaries that
   were already reproduced in Appendix B of revision 1.
3. Statements of revision 1 that the evidence does not support (for example
   "excellent agreement" and "true mesh independence") are replaced by what was
   actually measured.
4. A code review of the UDF is added (Section 7). It identifies two
   implementation issues that affect the physics. They are reported here, not
   corrected, because the reported results were produced with the code as it
   stands.

No simulation was re-run for this revision.

## 1. Objective

The study reproduces the 90° pipe-bend deposition case of Dehbi (2008), who
coupled a RANS flow solution with a continuous random walk (CRW) model based on
a normalised Langevin equation. Standard eddy-interaction (discrete random walk,
DRW) models in RANS solvers assume locally homogeneous, isotropic turbulence.
In wall-bounded flows they violate the well-mixed condition: particles with
low inertia drift towards regions of low turbulence intensity and accumulate
artificially near the wall, which over-predicts deposition. Dehbi's model
counters this with anisotropic near-wall velocity statistics and a drift term
proportional to the wall-normal gradient of the velocity variance.

The specific objectives were:

1. Compute the turbulent flow in the bend with the Reynolds stress model (RSM).
2. Implement the CRW–Langevin model as a compiled ANSYS Fluent UDF that
   replaces the built-in DRW dispersion model.
3. Compute the deposition fraction for particle diameters of 5 to 50 µm and
   compare it with Dehbi (2008, Fig. 8) and the experiments of Pui et al.
   (1987).

## 2. Geometry and flow conditions

| Quantity | Value | Source |
| --- | --- | --- |
| Pipe diameter $D$ | 0.02 m | revision 1, Sec. 2 |
| Bend curvature ratio $R/r$ | 4 ($R$ = 0.04 m) | revision 1, Sec. 2 |
| Straight inlet / outlet sections | 5$D$ (0.10 m) each | revision 1, Sec. 2 |
| Fluid | air, $\rho$ = 1.225 kg m⁻³, $\mu$ = 1.81 × 10⁻⁵ Pa s | revision 1, Sec. 4.1 |
| Inlet | uniform velocity, $U$ = 3.0 m s⁻¹ | revision 1, Sec. 4.1 |
| Reynolds number $Re = \rho U D/\mu$ | 4.06 × 10³ (computed from the values above; revision 1 states 4080) | computed |
| Outlet | pressure outlet, 0 Pa gauge | revision 1, Sec. 4.1 |
| Walls | stationary, no slip | revision 1, Sec. 4.1 |

![Computational domain](../figures/geometry.png)

*Figure 1. Computational domain: 90° bend with straight inlet and outlet sections.*

**Open item.** Revision 1 states that these conditions replicate Pui et al.
(1987). The geometry and Reynolds number have not been cross-checked against
the case definition in Dehbi (2008), Sec. 10, for this revision. Because the
comparison in Section 6 is made at equal particle diameter, any difference in
$D$, $U$ or $R/r$ from the benchmark case changes the Stokes number at a given
diameter and must be resolved before the comparison is considered final.

## 3. Mesh

The domain was meshed with structured hexahedral cells. The first cell height is
70 µm, which gives a wall $y^+ \approx 0.95$ (Figure 3). The near-wall region
contains 30 inflation layers with a growth ratio of 1.1. The final mesh has
about 5.3 × 10⁶ cells; a baseline mesh of about 2.5 × 10⁶ cells was used first.

| Refined mesh, 5.3 × 10⁶ cells | Baseline mesh, 2.5 × 10⁶ cells |
| --- | --- |
| ![refined](../figures/mesh-refined-5.3M.png) | ![baseline](../figures/mesh-baseline-2.5M.png) |

*Figure 2. Cross-section of the bend. The refined mesh adds resolution in the
bend and in the region $y^+ < 100$.*

![Wall y+](../figures/wall-yplus.png)

*Figure 3. Wall $y^+$ on the bend surface (colour scale 0 to 5).*

Mesh quality histograms (orthogonal quality and skewness) are in
[`figures/mesh-orthogonal-quality.png`](../figures/mesh-orthogonal-quality.png)
and [`figures/mesh-skewness.png`](../figures/mesh-skewness.png).

**Mesh sensitivity.** Revision 1 reports that the baseline mesh gave
deposition fractions for 5 and 10 µm particles that deviated from the
benchmark, while larger particles were insensitive, and that the refined mesh
removed the deviation. The deposition fractions obtained on the baseline mesh
were not recorded in revision 1. The study therefore demonstrates mesh
sensitivity for small particles but does not quantify grid convergence;
the phrase "true mesh independence" used in revision 1 is withdrawn.

## 4. Solver settings

| Item | Setting |
| --- | --- |
| Solver | ANSYS Fluent 2024 R2, pressure based, steady |
| Turbulence model | Reynolds stress model, linear pressure-strain, enhanced wall treatment (Figure S1) |
| Pressure–velocity coupling | SIMPLE |
| Spatial discretisation | least-squares gradients; second-order pressure; second-order upwind for momentum, $k$, $\varepsilon$ and Reynolds stresses (Figure S2) |
| Convergence | scaled residuals of all equations levelled off between about 10⁻³ and 10⁻⁶ after 1000 iterations (Figure 4) |
| Discrete phase | inert spheres, one-way coupling, steady tracking |
| Injection | surface injection at the inlet, 10,000 particles per diameter, randomised starting points |
| Particle density | 1000 kg m⁻³ (revision 1, Sec. 4.3) |
| DPM tracking | accuracy control (tolerance 10⁻⁵), automated high-/low-order scheme (Runge–Kutta / implicit) (Figure S4) |
| Built-in turbulent dispersion | disabled; dispersion supplied by the UDF body force |
| UDF hooks | body force `crw_langevin_force`, time step `crw_adaptive_dt`, 8 particle user scalars (Figure S3); wall boundary condition `dehbi_wall_bc`; on-demand `init_wall_udm` |
| Deposition criterion | particle trapped at first wall contact; deposition fraction = aborted / tracked |

![Residuals](../figures/residuals.png)

*Figure 4. Scaled residuals of the continuous phase.*

The Fluent dialogs reproduced in revision 1 are kept as supplementary figures
in [`figures/setup/`](../figures/setup/): S1 viscous model
(`viscous-model.png`), S2 solution methods (`solution-methods.png`), S3 DPM
UDF hooks (`dpm-udf-hooks.png`), S4 DPM numerics (`dpm-numerics.png`),
S5 injection (`injection.png`) and S6 UDF compilation environment
(`udf-compilation.png`).

**Settings to confirm.** The injection dialog in revision 1
([`figures/setup/injection.png`](../figures/setup/injection.png)) shows three
values that differ from the text and should be confirmed from the case file:

1. The injection named `particles-10micron` has a diameter of 5 × 10⁻⁶ m.
2. The particle material is `anthracite`, whose Fluent default density is
   1550 kg m⁻³, while the text states 1000 kg m⁻³. The value used depends on
   whether the material density was edited; the dialog does not show it.
3. The particle injection velocity is 0 m s⁻¹, while the text and the UDF
   set-up notes state that particles are released at the gas velocity. With
   $\tau_p \le 7.7$ ms and a 0.10 m inlet section, particles relax to the gas
   velocity before the bend, so the effect on deposition is expected to be
   small, but the text should match the setting.

## 5. CRW–Langevin model as implemented

Symbols: $u_\tau$ friction velocity, $\nu$ kinematic viscosity, $y^+ = y u_\tau/\nu$,
$k$ turbulent kinetic energy, $\varepsilon$ dissipation rate, $\tau_p$
particle relaxation time, $\tau_L$ fluid Lagrangian integral time,
$\xi$ standard normal deviate. Equation numbers in square brackets are those
cited in the source code for Dehbi (2008).

### 5.1 Particle equation of motion

Fluent integrates the particle momentum equation with Schiller–Naumann drag
computed from the mean slip velocity. The UDF adds the drag of the fluid
velocity fluctuation as a body force per unit mass (lines 644–738):

$$
\frac{F_i}{m_p} = f_D\,u'_i, \qquad
f_D = \frac{1 + 0.15\,Re_p^{0.687}}{\tau_p}, \qquad
\tau_p = \frac{\rho_p d_p^2}{18\,\mu}, \qquad
Re_p = \frac{|\mathbf{U} - \mathbf{u}_p|\,d_p}{\nu}.
$$

The fluctuation $u'_i$ is updated once per particle time step, on the call for
the first velocity component, and reused for the other two.

### 5.2 Normalised Langevin equation

The model evolves the normalised fluctuation $\eta_i = u'_i/\sigma_i$:

$$
d\eta_i = -\frac{\eta_i}{\tau_L}\,dt + A_i\,dt + \sqrt{\frac{2}{\tau_L}}\,dW_i .
$$

It is discretised semi-implicitly: the relaxation term is implicit, and the
drift and noise terms are explicit (lines 397–399, 427–429):

$$
\eta_i^{n+1} = \frac{\eta_i^{n} + \sqrt{2\,\Delta t/\tau_L}\;\xi_i + A_i\,\Delta t}{1 + \Delta t/\tau_L}.
$$

The code comments describe this scheme as "implicit Euler"; it is more
precisely a semi-implicit (linearly implicit) Euler–Maruyama scheme. For
$\Delta t \ll \tau_L$, which holds here because $\Delta t \le 10^{-6}$ s, the
stationary variance of $\eta$ is $1/(1 + \Delta t/2\tau_L) \approx 1$.

### 5.3 Near-wall branch ($y^+ < 100$)

In the near-wall branch the equations are solved in a body-fitted frame
(lines 225–275): $\mathbf{e}_2$ is the inward wall normal, $\mathbf{e}_1$ the
local mean-velocity direction made orthogonal to $\mathbf{e}_2$, and
$\mathbf{e}_3 = \mathbf{e}_1 \times \mathbf{e}_2$. The root-mean-square
fluctuations are the DNS fits attributed to Dreeben and Pope (1997)
[Eqs. 37–39] (lines 131–148):

$$
\frac{\sigma_1}{u_\tau} = \frac{0.40\,y^+}{1 + 0.0239\,(y^+)^{1.496}}, \quad
\frac{\sigma_2}{u_\tau} = \frac{0.0116\,(y^+)^2}{1 + 0.203\,y^+ + 0.00140\,(y^+)^{2.421}}, \quad
\frac{\sigma_3}{u_\tau} = \frac{0.19\,y^+}{1 + 0.0361\,(y^+)^{1.322}} .
$$

The Lagrangian time scale follows Kallio and Reeks (1989) [Eqs. 31–33]
(lines 196–214):

$$
\tau_L = \tau_L^+\,\frac{\nu}{u_\tau^2}, \qquad
\tau_L^+ =
\begin{cases}
10, & y^+ \le 5,\\
7.122 + 0.5731\,y^+ - 0.00129\,(y^+)^2, & 5 < y^+ \le 100.
\end{cases}
$$

Only the wall-normal component carries a drift term [Eq. 22]
(lines 421–422), reduced for inertial particles by the local Stokes number
$St = \tau_p/\tau_L$:

$$
A_1 = A_3 = 0, \qquad
A_2 = \frac{1}{1 + St}\,\frac{\partial \sigma_2}{\partial y}
    = \frac{1}{1 + St}\,\frac{u_\tau^2}{\nu}\,\frac{d(\sigma_2/u_\tau)}{dy^+}.
$$

The derivative $d(\sigma_2/u_\tau)/dy^+$ is evaluated analytically
(lines 163–188). The fluctuation is then $u'_i = \eta_i \sigma_i$ in the
body-fitted frame and is rotated to Cartesian components.

### 5.4 Bulk branch ($y^+ \ge 100$)

Away from the wall the turbulence is treated as isotropic
(lines 444–474), with $\sigma = (2k/3)^{1/2}$ [Eq. 26] and

$$
\tau_L = \frac{2}{C_0}\,\frac{k}{\varepsilon}, \qquad C_0 = 14 \quad \text{[Eq. 36]}.
$$

The bulk drift term of the model, proportional to $\nabla k$ [Eq. 27], is set
to zero in this implementation (modification FIX-17 in the code). The code
gives the reason that reading the cell gradient of $k$ under the RSM caused
a segmentation fault. The effect of omitting this term was not quantified.

### 5.5 Wall data and time step

The on-demand function `init_wall_udm` (lines 498–625) stores $u_\tau$,
$y^+$ and the wall normal in user-defined memory. It computes the wall shear
stress from the tangential velocity of the wall-adjacent cell divided by the
wall distance. The particle time step is
$\Delta t = \max\left[\min(10^{-6}\ \text{s},\ 0.1\,\tau_p),\ 10^{-9}\ \text{s}\right]$,
capped by Fluent's own step (lines 765–797).

## 6. Results

### 6.1 Continuous phase

The RSM solution shows the secondary flow expected in a curved pipe: the
centrifugal pressure gradient drives faster core fluid towards the outer wall
and slower near-wall fluid inwards, forming a pair of counter-rotating Dean
vortices. Revision 1 describes this flow qualitatively; no secondary-flow
contours or profiles are included in either revision.

### 6.2 Deposition fraction

Table 1 lists the deposition fractions taken from the Fluent DPM summaries
reproduced in Appendix B of revision 1 (also in
[`data/deposition_fraction.csv`](../data/deposition_fraction.csv)). The
Stokes number is $Stk = \rho_p d_p^2 U/(9\mu D)$ with the values of Section 2,
neglecting the Cunningham correction.

*Table 1. Deposition fraction (DF) per particle diameter.*

| $d_p$ (µm) | $Stk$ | Tracked | Escaped | Trapped | DF (%) | Sampling standard error (%) |
| --- | --- | --- | --- | --- | --- | --- |
| 10 | 0.092 | 10,000 | 8,329 | 1,671 | 16.7 | 0.37 |
| 20 | 0.368 | 10,000 | 4,279 | 5,721 | 57.2 | 0.49 |
| 30 | 0.829 | 10,000 | 2,142 | 7,858 | 78.6 | 0.41 |
| 40 | 1.473 | 10,000 | 546 | 9,454 | 94.5 | 0.23 |
| 50 | 2.302 | 10,000 | 199 | 9,801 | 98.0 | 0.14 |

The standard error, $\sqrt{DF(1-DF)/N}$, covers only the finite number of
particles. It does not include the effect of the random-number issue described
in Section 7.2.

Revision 1 (Fig. 10, kept as
[`figures/original-figure10-deposition.png`](../figures/original-figure10-deposition.png))
also shows points at 5 µm (about 5 %) and 15 µm (about 33 %). No DPM summary
is available for these two diameters, so they are not included in Table 1.

![Deposition fraction](../figures/deposition-fraction.png)

*Figure 5. Deposition fraction from Table 1, generated by
[`tools/plot_deposition.py`](../tools/plot_deposition.py).*

### 6.3 Comparison with the benchmark

Revision 1 plots the simulated points together with a solid curve labelled
"CRW Model", described in the text as the numerical curve of Dehbi (2008).
The curve passes through every simulated point, including changes of slope at
15 and 20 µm, and no experimental points of Pui et al. (1987) are shown.
The figure therefore does not allow the reader to judge the agreement
independently, and the statement "excellent agreement" is withdrawn until the
following comparison is made:

1. Digitise Dehbi (2008, Fig. 8) and the Pui et al. (1987) data, recording the
   digitisation uncertainty.
2. Plot them against the values in Table 1, at matched Stokes number if the
   flow conditions differ (see the open item in Section 2).
3. Report the difference at each diameter.

## 7. Code review of the UDF

The source file was reviewed line by line for this revision. The file in this
repository is unchanged from the version used for the results above.

### 7.1 The near-wall branch is applied only in wall-adjacent cells

`init_wall_udm` loops over wall faces and writes wall data only into the cell
adjacent to each face (`F_C0`, lines 522–618). Every other cell keeps the flag
`UDM_INIT = 0`. In `advance_langevin`, such cells fall back to the bulk branch:
lines 376–383 compute a $y^+$ from the wall distance and then overwrite it with
`BL_YPLUS_MAX + 1`.

As a result, the anisotropic near-wall statistics and the drift correction of
Section 5.3 act only in the first cell layer, where the cell-centre
$y^+ \approx 0.5$. With 30 inflation layers, the cells between $y^+ \approx 1$
and $y^+ = 100$ use the isotropic, drift-free bulk branch, which is the
treatment the Dehbi model was designed to replace. The reported agreement
therefore cannot be attributed to the drift correction as implemented. The
observed sensitivity to near-wall refinement (Section 3) may also be partly
caused by this behaviour, because refinement changes the thickness of the
single layer in which the near-wall branch is active.

**Recommended correction:** store $u_\tau$, $y^+$ and the wall normal for all
cells with $y^+ \le 100$, for example by assigning each cell the data of its
nearest wall face and using `C_WALL_DIST` for the distance, then re-run all
diameters.

### 7.2 The random-number generator produces correlated deviates

`gauss_rand` (lines 94–121) does not keep a random-number state. It hashes the
particle's memory address, the current particle time and a salt through one
step of two linear congruential generators and applies the Box–Muller
transform. The program [`tools/rng_check.c`](../tools/rng_check.c) runs a copy
of this function with 32-bit `unsigned long`, as on the MSVC compiler used for
the build. For a time step of 10⁻⁶ s and 2 × 10⁵ samples it gives:

| Statistic | Ideal value | Result |
| --- | --- | --- |
| Mean, variance | 0, 1 | +0.001, 1.003 |
| Autocorrelation, lags 1 to 5 | 0 ± 0.005 | +0.22, +0.06, −0.11, −0.29, −0.23 |
| Correlation between two particles at the same time (addresses 0x1C0 bytes apart) | 0 ± 0.005 | +0.13 |

The one-point statistics are correct, but successive deviates for the same
particle are strongly correlated, and particles stored close together in
memory receive correlated forcing. Serial correlation in $\xi$ changes the
effective diffusivity of the Langevin process, so the dispersion of particles
is not the one the model prescribes.

**Recommended correction:** give each particle its own generator state, for
example a 32-bit xorshift or SplitMix state stored in a particle user scalar
and seeded once at injection, and re-run `rng_check` on the new generator
before re-running the simulations.

### 7.3 Smaller points

1. The deposition comment (lines 741–752) states that a particle reaching a
   wall face is within one radius of the wall. Fluent detects wall contact when
   the particle centre reaches the face, so the trap is applied at zero
   distance, not at one radius. The difference matters most for the largest
   particles near the wall.
2. The bulk-branch comment states that omitting the $\nabla k$ drift has a
   "negligible effect on the final results". This was not tested.
3. In the fallback branch (lines 378–380) the computed $y^+$ is immediately
   overwritten; the two lines have no effect.
4. The particle Reynolds number in the drag factor uses only the first
   `ND_ND` components (3 in this 3D case) and the mean slip velocity, which is
   the usual approximation.

The code is otherwise well organised: it uses named constants for memory
indices, guards against unassigned cells at injection, follows C89 declaration
rules required by the Fluent compiler, and records each modification (FIX-2 to
FIX-19) in comments.

## 8. Conclusions

1. A CRW–Langevin dispersion model following Dehbi (2008) was implemented as a
   compiled ANSYS Fluent UDF with anisotropic near-wall statistics, a
   body-fitted frame, a Stokes-number-weighted wall-normal drift and an
   adaptive particle time step.
2. On an RSM flow field with wall $y^+ \approx 0.95$, the UDF gives deposition
   fractions that rise from 16.7 % at 10 µm to 98.0 % at 50 µm (Table 1).
3. Deposition of 5 and 10 µm particles was sensitive to the near-wall mesh;
   this sensitivity was observed but not quantified.
4. The code review shows that, in the version used, the near-wall branch acts
   only in wall-adjacent cells and the random deviates are correlated. Until
   both are corrected and the benchmark comparison of Section 6.3 is made, the
   results should be described as a working implementation of the model with
   preliminary validation, not as a validated reproduction.

## 9. Use of software tools

As stated in revision 1, the Cursor code editor was used during UDF
development to debug C syntax, memory-access errors and compilation problems
with Fluent's `udf.h` macros, and to tidy comments. The physical model,
equations and coupling logic were taken from Dehbi (2008). The analysis in
Section 7 of this revision was carried out with the help of an AI assistant and
checked against the source code.

## References

1. Dehbi, A. (2008). Turbulent particle dispersion in arbitrary wall-bounded
   geometries: A coupled CFD-Langevin-equation based approach.
   *International Journal of Multiphase Flow*, 34(9), 819–828.
   <https://doi.org/10.1016/j.ijmultiphaseflow.2008.03.001>
2. Pui, D. Y. H., Romay-Novas, F., & Liu, B. Y. H. (1987). Experimental study
   of particle deposition in bends of circular cross section.
   *Aerosol Science and Technology*, 7(3), 301–315.
3. Dreeben, T. D., & Pope, S. B. (1997). Probability density function and
   Reynolds-stress modeling of near-wall turbulent flows.
   *Physics of Fluids*, 9(1), 154–163.
4. Kallio, G. A., & Reeks, M. W. (1989). A numerical simulation of particle
   deposition in turbulent boundary layers.
   *International Journal of Multiphase Flow*, 15(3), 433–446.
