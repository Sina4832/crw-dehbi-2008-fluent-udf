---
title: "Implementation of the Dehbi (2008) CRW–Langevin Particle Dispersion Model as an ANSYS Fluent UDF"
subtitle: "Technical report: 90° pipe-bend deposition case"
author: "Sina Karbakhsh Ravari, Shiraz University. Supervisor: Dr. Ebrahim Ghahramani"
date: "November 2026"
---

## Abstract

The continuous random walk (CRW) model of Dehbi (2008) was implemented as a
compiled ANSYS Fluent user-defined function (UDF) and applied to turbulent
aerosol deposition in a 90° pipe bend at $Re \approx 4.1 \times 10^3$. The UDF
integrates a normalised Langevin equation for the fluid velocity seen by each
particle, with anisotropic near-wall statistics and a wall-normal drift
correction, and couples it to Fluent's discrete phase model through a body
force. On a Reynolds-stress-model flow field, the computed deposition fraction
rises from 16.7 % at a flow Stokes number of 0.09 to 98.0 % at 2.3. A
line-by-line review and a statistical test of the random-number generator
identify three aspects of the implementation that define the scope of the
present results.

## 1. Introduction

Discrete random walk (DRW) models sample velocity fluctuations from locally
isotropic turbulence. In inhomogeneous near-wall turbulence they produce a
spurious drift of low-inertia particles towards the wall, violating the
well-mixed condition and over-predicting small-particle deposition. Dehbi
(2008) replaced the DRW with a normalised Langevin equation that uses
anisotropic near-wall rms profiles and a drift term proportional to the
wall-normal gradient of the rms velocity. This work implements that model in
ANSYS Fluent and computes the deposition fraction in a 90° bend for particle
diameters of 10 to 50 µm.

## 2. Case set-up

*Table 1. Geometry, flow conditions and solver settings.*

| Item | Value |
| --- | --- |
| Pipe diameter $D$; bend curvature ratio $R/r$ | 0.02 m; 4 |
| Straight inlet and outlet sections | 5$D$ each |
| Fluid | air, $\rho$ = 1.225 kg m⁻³, $\mu$ = 1.81 × 10⁻⁵ Pa s |
| Inlet; outlet | uniform velocity 3.0 m s⁻¹; pressure outlet |
| Reynolds number $\rho U D/\mu$; friction Reynolds number $u_\tau r/\nu$ | 4.06 × 10³; ≈ 140 (Blasius estimate) |
| Mesh | structured hexahedral, ≈ 5.3 × 10⁶ cells, wall $y^+ \approx 1$ |
| Solver | ANSYS Fluent 2024 R2, pressure based, steady, SIMPLE |
| Turbulence model | Reynolds stress model, linear pressure–strain, enhanced wall treatment |
| Discretisation | least-squares gradients, second-order pressure, second-order upwind |
| Discrete phase | inert spheres, one-way coupling, steady tracking, 10,000 particles per diameter, built-in dispersion disabled |
| Deposition | particle stopped at first wall contact (reported by Fluent as "aborted") |

![Computational domain](../figures/geometry.png)

*Figure 1. Computational domain.*

| Refined mesh, 5.3 × 10⁶ cells | Baseline mesh, 2.5 × 10⁶ cells |
| --- | --- |
| ![refined](../figures/mesh-refined-5.3M.png) | ![baseline](../figures/mesh-baseline-2.5M.png) |

*Figure 2. Bend cross-section of the two meshes. Small-particle deposition was
sensitive to the near-wall resolution; the refined mesh was used for all
results.*

![Wall y+](../figures/wall-yplus.png)

*Figure 3. Wall $y^+$ on the bend surface.*

The scaled residuals of the continuous phase decreased to between 10⁻³ and
2 × 10⁻⁶ ([`figures/residuals.png`](../figures/residuals.png)). The Fluent
dialogs for the turbulence model, solution methods, discrete-phase numerics,
UDF hooks, injection and compilation are archived in
[`figures/setup/`](../figures/setup/).

## 3. Model implementation

Line numbers refer to [`udf/crw_dehbi2008_rev10.c`](../udf/crw_dehbi2008_rev10.c).

**Coupling.** Fluent integrates the particle momentum equation with the mean
slip velocity. The UDF adds the drag of the fluid velocity fluctuation
$u'_i$ as a body force per unit mass (lines 644–738):

$$
\frac{F_i}{m_p} = f_D\,u'_i, \qquad
f_D = \frac{1 + 0.15\,Re_p^{0.687}}{\tau_p}, \qquad
\tau_p = \frac{\rho_p d_p^2}{18\,\mu}.
$$

**Langevin equation.** The normalised fluctuation $\eta_i = u'_i/\sigma_i$
obeys $d\eta_i = -\eta_i\,dt/\tau_L + A_i\,dt + \sqrt{2/\tau_L}\,dW_i$ and is
advanced with a semi-implicit Euler–Maruyama scheme (lines 397–399, 427–429):

$$
\eta_i^{n+1} = \frac{\eta_i^{n} + \sqrt{2\,\Delta t/\tau_L}\;\xi_i + A_i\,\Delta t}{1 + \Delta t/\tau_L}.
$$

**Near-wall branch ($y^+ < 100$).** The equation is solved in a body-fitted
frame (wall normal, mean-flow direction, binormal; lines 225–275). The rms
profiles are the DNS fits given by Dehbi (2008) (lines 131–148):

$$
\frac{\sigma_1}{u_\tau} = \frac{0.40\,y^+}{1 + 0.0239\,(y^+)^{1.496}}, \quad
\frac{\sigma_2}{u_\tau} = \frac{0.0116\,(y^+)^2}{1 + 0.203\,y^+ + 0.00140\,(y^+)^{2.421}}, \quad
\frac{\sigma_3}{u_\tau} = \frac{0.19\,y^+}{1 + 0.0361\,(y^+)^{1.322}}.
$$

The Lagrangian time scale follows Kallio and Reeks (1989):
$\tau_L^+ = 10$ for $y^+ \le 5$ and
$\tau_L^+ = 7.122 + 0.5731\,y^+ - 0.00129\,(y^+)^2$ otherwise
(lines 196–214). Only the wall-normal component carries a drift, reduced by
the local particle Stokes number $St = \tau_p/\tau_L$; its derivative is
evaluated analytically (lines 163–188, 421–422):

$$
A_2 = \frac{1}{1 + St}\,\frac{\partial \sigma_2}{\partial y}, \qquad A_1 = A_3 = 0.
$$

**Bulk branch ($y^+ \ge 100$).** Isotropic turbulence with
$\sigma = (2k/3)^{1/2}$ and $\tau_L = (2/C_0)\,k/\varepsilon$, $C_0 = 14$
(lines 444–474). The $\nabla k$ drift term is omitted, because accessing the
cell gradient of $k$ under the RSM caused a solver failure.

**Wall data and time step.** The on-demand function `init_wall_udm`
(lines 498–625) stores $u_\tau$, $y^+$ and the wall normal in user-defined
memory. The particle time step is $\min(10^{-6}\ \text{s},\ 0.1\,\tau_p)$
(lines 765–797).

## 4. Results

The flow Stokes number is $Stk = \rho_p d_p^2 U/(9\mu D)$, evaluated with
$\rho_p$ = 1000 kg m⁻³. The standard error is $\sqrt{DF(1-DF)/N}$.

*Table 2. Deposition fraction (DF) from the Fluent DPM summaries
([`data/deposition_fraction.csv`](../data/deposition_fraction.csv)).*

| $d_p$ (µm) | $Stk$ | Tracked | Escaped | Deposited | DF (%) | Standard error (%) |
| --- | --- | --- | --- | --- | --- | --- |
| 10 | 0.092 | 10,000 | 8,329 | 1,671 | 16.7 | 0.37 |
| 20 | 0.368 | 10,000 | 4,279 | 5,721 | 57.2 | 0.49 |
| 30 | 0.829 | 10,000 | 2,142 | 7,858 | 78.6 | 0.41 |
| 40 | 1.473 | 10,000 | 546 | 9,454 | 94.5 | 0.23 |
| 50 | 2.302 | 10,000 | 199 | 9,801 | 98.0 | 0.14 |

![Deposition fraction](../figures/deposition-fraction.png)

*Figure 4. Deposition fraction against particle diameter and flow Stokes
number, generated by [`tools/plot_deposition.py`](../tools/plot_deposition.py).*

The deposition fraction increases monotonically with Stokes number, from
inertia-limited deposition at $Stk < 0.1$ to near-complete impaction at
$Stk > 1.5$, which is the expected qualitative behaviour for bend flows.
Deposition is counted on all wall zones.

## 5. Code verification

**Random-number generator.** The function `gauss_rand` (lines 94–121)
generates Box–Muller deviates from a hash of the particle address and the
particle time, without a stored state. The program
[`tools/rng_check.c`](../tools/rng_check.c) runs an exact copy of it.

*Table 3. Statistics of the UDF generator at $\Delta t$ = 10⁻⁶ s.*

| Statistic | Ideal | Result |
| --- | --- | --- |
| Mean, variance (2 × 10⁵ samples) | 0, 1 | +0.001, 1.003 |
| Autocorrelation, lags 1 to 5 | 0 | +0.22, +0.06, −0.11, −0.29, −0.23 |
| Variance of 1000 accumulated deviates divided by 1000, 2000 particle addresses | 1 | median 0.48; 5–95 % range 0.20 to 2.7 |

The one-point statistics are correct, while successive deviates are
correlated, so the effective turbulent diffusivity varies between particles.

**Implementation scope.** The review also establishes that:

1. the near-wall statistics are applied in the wall-adjacent cell layer and
   evaluated at the cell centre; other cells use the bulk branch;
2. the Langevin state is stored in the body-fitted frame near the wall and in
   Cartesian components in the bulk, without rotation between the two.

## 6. Limitations and further work

The present case differs from the turbulent bend experiment of Pui et al.
(1987), reported as $D$ = 8.51 mm, curvature ratio 5.6 and $Re$ = 10,000
(Sakib et al., 2023); a comparison with that benchmark is meaningful only at
matched Stokes number. The 5$D$ inlet section does not produce fully developed
flow at the bend. The archived injection dialog for the smallest size class
specifies a diameter of 5 µm and the material anthracite, whereas Table 2
assumes 10 µm and 1000 kg m⁻³, and the gravity setting is not recorded in the
archived dialogs. Further work comprises a per-particle generator
state, near-wall statistics evaluated at the particle position across the
whole $y^+ < 100$ region, a consistent frame transformation, a well-mixed
tracer test, and a comparison with the Pui et al. data at matched Stokes
number.

## 7. Conclusions

A CRW–Langevin dispersion model following Dehbi (2008) was implemented as a
compiled ANSYS Fluent UDF with anisotropic near-wall statistics, a body-fitted
frame, a Stokes-number-weighted drift correction and an adaptive particle
time step. Applied to a 90° bend, it gives deposition fractions from 16.7 % to
98.0 % over $0.09 \le Stk \le 2.3$. The accompanying verification tools
quantify the statistical properties of the random forcing and define the
improvements required for a quantitative benchmark comparison.

## Acknowledgement of software tools

The Cursor editor was used to debug C compilation errors during UDF
development.

## References

1. Dehbi, A. (2008). Turbulent particle dispersion in arbitrary wall-bounded
   geometries: A coupled CFD-Langevin-equation based approach.
   *International Journal of Multiphase Flow*, 34(9), 819–828.
   <https://doi.org/10.1016/j.ijmultiphaseflow.2008.03.001>
2. Kallio, G. A., & Reeks, M. W. (1989). A numerical simulation of particle
   deposition in turbulent boundary layers.
   *International Journal of Multiphase Flow*, 15(3), 433–446.
3. Pui, D. Y. H., Romay-Novas, F., & Liu, B. Y. H. (1987). Experimental study
   of particle deposition in bends of circular cross section.
   *Aerosol Science and Technology*, 7(3), 301–315.
4. Sakib, M. N., Shuvo, M. S., Rahman, R., & Saha, S. (2023). Particle
   deposition and fluid flow characteristics in turbulent corrugated pipe flow
   using Eulerian-Lagrangian approach. *Heliyon*, 9(3), e14603.
   <https://doi.org/10.1016/j.heliyon.2023.e14603>
