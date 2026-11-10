# Implementation of the Dehbi (2008) Continuous Random Walk Model as an ANSYS Fluent UDF

**Technical report: aerosol deposition in a 90° pipe bend**

Sina Karbakhsh Ravari, Shiraz University · Project supervisor: Dr. Ebrahim Ghahramani

Revision 3, November 2026; supersedes the project report of March 2026.

## Abstract

The continuous random walk (CRW) model of Dehbi (2008), which integrates a
normalised Langevin equation for the fluid velocity seen by a particle, was
implemented as a compiled ANSYS Fluent user-defined function (UDF). It was
applied to aerosol deposition in a 90° pipe bend (diameter 20 mm, curvature
ratio 4, Reynolds number 4.1 × 10³) on a Reynolds-stress-model flow field,
with 10⁴ particles per size class. The computed deposition fraction increases
from 16.7 % for the smallest class (nominally 10 µm; Stokes number 0.09) to
98.0 % at 50 µm (Stokes number 2.3). At equal Stokes number it differs from
the bend correlation of Pui et al. (1987) by −5.5 to +1.4 percentage points;
because the geometry and flow conditions differ from those of the
experiments, this comparison is qualitative. Source inspection and a
statistical test identify three departures from the published model: the
near-wall statistics and the drift correction act only in wall-adjacent
cells, the Langevin state is not transformed between reference frames, and
successive random deviates are correlated.

## 1. Introduction

Discrete random walk (DRW) models sample velocity fluctuations from locally
isotropic turbulence. In inhomogeneous near-wall turbulence they produce a
spurious drift of low-inertia particles towards the wall, which violates the
well-mixed condition (Thomson, 1987) and leads to over-prediction of
small-particle deposition (Dehbi, 2008). Dehbi (2008) replaced the DRW with a
normalised Langevin equation that uses anisotropic near-wall rms profiles and
a drift term proportional to the wall-normal gradient of the rms velocity.
This report documents an implementation of that model as a compiled ANSYS
Fluent UDF, its application to a 90° pipe bend for particle diameters of
10 to 50 µm, a comparison with the bend correlation of Pui et al. (1987), and
the points at which the implementation departs from the published model.

## 2. Case set-up

*Table 1. Geometry, flow conditions and solver settings.*

| Item | Value |
| --- | --- |
| Pipe diameter $D$; curvature ratio $R_0 = R/r$ (bend radius over pipe radius) | 20 mm; 4 ($R$ = 40 mm) |
| Straight inlet and outlet sections | $5D$ each |
| Fluid | air, $\rho$ = 1.225 kg m⁻³, $\mu$ = 1.81 × 10⁻⁵ Pa s |
| Inlet; outlet | uniform velocity $U$ = 3.0 m s⁻¹; pressure outlet |
| Reynolds number $\rho U D/\mu$; friction Reynolds number $u_\tau r/\nu$ | 4.06 × 10³; ≈ 140 (Blasius estimate) |
| Mesh | structured hexahedral, ≈ 5.3 × 10⁶ cells, wall $y^+ \approx 1$ |
| Solver | ANSYS Fluent 2024 R2, pressure based, steady, SIMPLE |
| Turbulence model | Reynolds stress model (RSM), linear pressure–strain, enhanced wall treatment |
| Discretisation | least-squares gradients, second-order pressure, second-order upwind |
| Discrete phase | inert spheres, one-way coupling, steady tracking, 10⁴ particles per size class, built-in dispersion model disabled |
| Particles | density 1000 kg m⁻³, released at the gas velocity (UDF header; see Section 6) |
| Deposition | particle stopped at first wall contact (reported by Fluent as "aborted") |

![Computational domain](../figures/geometry.png)

*Figure 1. Computational domain: straight inlet of length 5D, 90° bend of
curvature ratio R<sub>0</sub> = 4 and straight outlet of length 5D; D = 20 mm.*

| Refined mesh, 5.3 × 10⁶ cells | Baseline mesh, 2.5 × 10⁶ cells |
| --- | --- |
| ![refined](../figures/mesh-refined-5.3M.png) | ![baseline](../figures/mesh-baseline-2.5M.png) |

*Figure 2. Bend cross-sections of the refined and baseline meshes (screenshots
not at the same scale). The refined mesh was used for all results; results on
the baseline mesh are not archived.*

![Wall y+](../figures/wall-yplus.png)

*Figure 3. Wall y<sup>+</sup> on the bend surface; the values lie below about 1.2
(colour scale 0–5).*

After about 950 iterations the scaled residuals of the continuous phase were
between 7 × 10⁻⁴ (ε) and 2 × 10⁻⁶ and still decreasing slowly
([`figures/residuals.png`](../figures/residuals.png); the legend of the
archived screenshot is truncated). No integral-quantity monitors were
archived. The Fluent dialogs for the turbulence model, solution methods,
discrete-phase numerics, UDF hooks and injection, and the start of Fluent from
the Visual Studio command prompt used to compile the UDF, are archived in
[`figures/setup/`](../figures/setup/).

## 3. Model and implementation

Line numbers refer to
[`udf/crw_dehbi2008_rev10.c`](../udf/crw_dehbi2008_rev10.c). Its executable
code is identical to the file first uploaded to this repository (commit
`dfa9cad`, May 2026); only comments have been revised since.

**Coupling.** Fluent integrates the particle momentum equation with the mean
slip velocity. The UDF adds the drag of the fluid velocity fluctuation
$u'_i$ as a body force per unit mass (lines 644–738):

$$
\frac{F_i}{m_p} = f_D\,u'_i, \qquad
f_D = \frac{1 + 0.15\,Re_p^{0.687}}{\tau_p}, \qquad
\tau_p = \frac{\rho_p d_p^2}{18\,\mu},
$$

where $Re_p$ is the particle Reynolds number based on the mean slip velocity.
If Fluent's default spherical drag law (Morsi and Alexander, 1972) was used
for the mean drag, the two drag factors differ by up to about 6 % for
$0.1 < Re_p < 1$; the drag law selected in the case is not recorded.

**Langevin equation.** The normalised fluctuation $\eta_i = u'_i/\sigma_i$,
where $\sigma_i$ is the rms fluid velocity, obeys
$d\eta_i = -\eta_i\,dt/\tau_L + A_i\,dt + \sqrt{2/\tau_L}\,dW_i$. Here
$\tau_L$ is the Lagrangian integral time scale of the fluid seen by the
particle, $A_i$ is a drift term and $dW_i$ is a Wiener increment. The equation
is advanced with a semi-implicit Euler–Maruyama scheme (lines 397–399,
427–429), in which $\xi_i$ are standard normal deviates (Section 5):

$$
\eta_i^{n+1} = \frac{\eta_i^{n} + \sqrt{2\,\Delta t/\tau_L}\;\xi_i + A_i\,\Delta t}{1 + \Delta t/\tau_L}.
$$

**Near-wall branch (wall-adjacent cells).** Dehbi (2008) uses anisotropic
statistics for $y^+ < 100$ in a body-fitted frame. In the UDF, $e_2$ is the
wall normal pointing into the fluid, $e_1$ the local mean-flow direction made
orthogonal to $e_2$ (streamwise) and $e_3 = e_1 \times e_2$ (spanwise)
(lines 225–275). The on-demand function `init_wall_udm`
(lines 498–625) stores the friction velocity $u_\tau$, the $y^+$ of the cell
centroid and the wall normal for cells adjacent to a wall only. Every other
cell is assigned $y^+ = 101$ at run time and uses the bulk branch
(lines 374–383). The near-wall branch is therefore active only in the first
cell layer, where $y^+ \approx 0.5$–1 on the present mesh. Its rms profiles
are the DNS-based fits given by Dehbi (2008) (lines 131–148):

$$
\frac{\sigma_1}{u_\tau} = \frac{0.40\,y^+}{1 + 0.0239\,(y^+)^{1.496}},
$$

$$
\frac{\sigma_2}{u_\tau} = \frac{0.0116\,(y^+)^2}{1 + 0.203\,y^+ + 0.00140\,(y^+)^{2.421}},
$$

$$
\frac{\sigma_3}{u_\tau} = \frac{0.19\,y^+}{1 + 0.0361\,(y^+)^{1.322}}.
$$

The Lagrangian time scale follows Kallio and Reeks (1989):
$\tau_L^+ = \tau_L u_\tau^2/\nu = 10$ for $y^+ \le 5$ and
$\tau_L^+ = 7.122 + 0.5731\,y^+ - 0.00129\,(y^+)^2$ otherwise
(lines 196–214). Only the wall-normal component carries a drift, reduced by
the local Stokes number $St_L = \tau_p/\tau_L$; its derivative is evaluated
analytically (lines 163–188, 421–422):

$$
A_2 = \frac{1}{1 + St_L}\,\frac{\partial \sigma_2}{\partial y}, \qquad A_1 = A_3 = 0.
$$

**Bulk branch (all other cells).** Isotropic turbulence with
$\sigma = (2k/3)^{1/2}$ and $\tau_L = (2/C_0)\,k/\varepsilon$, $C_0 = 14$,
where $k$ is the turbulent kinetic energy and $\varepsilon$ its dissipation
rate (lines 444–474). Dehbi (2008) includes a drift proportional to $\nabla k$;
it is omitted here because accessing the cell gradient of $k$ under the RSM
caused a solver failure. Between the second cell layer ($y^+ \approx 1.5$) and
$y^+ = 100$, about 70 % of the pipe radius, the wall-normal fluctuation is
therefore modelled with the isotropic value and without drift. With $k$ taken
from the fits above, the isotropic value exceeds $\sigma_2$ by factors of
about 7 at $y^+ = 5$ and 2 at $y^+ = 30$. The well-mixed condition is not
enforced in this region, and because the extent of the near-wall branch is
set by the first-cell height, the model itself depends on the near-wall mesh.

**Wall interaction and time step.** `dehbi_wall_bc` (lines 749–752) stops a
particle when its trajectory reaches a wall face. The particle time step is
limited to $\min(10^{-6}\ \text{s},\ 0.1\,\tau_p)$ (lines 765–797), which is
10⁻⁶ s for all diameters considered.

*Table 2. Departures of the implementation from the published model.*

| Element | Dehbi (2008) | Present implementation |
| --- | --- | --- |
| Near-wall statistics and drift | applied for $y^+ < 100$ | applied only in wall-adjacent cells, at the cell centroid |
| Bulk drift | proportional to $\nabla k$ | omitted |
| Frame of the Langevin state | body-fitted near the wall, Cartesian in the bulk | the same, but the state is not transformed when a particle changes branch |
| Random forcing | independent Gaussian increments | deviates hashed from particle address and time; serially correlated (Section 5) |
| Deposition | particle centre within one particle radius of the wall | particle trajectory reaches a wall face |

## 4. Results

The Stokes number is $Stk = 2\tau_p U/D = \rho_p d_p^2 U/(9\mu D)$, based on
the bulk velocity and the pipe radius as in Pui et al. (1987), and is
evaluated with $\rho_p$ = 1000 kg m⁻³ and without Cunningham correction. The
deposition fraction is DF = $N_{dep}/N$, where $N$ is the number of tracked
particles and $N_{dep}$ the number reported as "aborted"; it equals the
deposition efficiency (one minus the penetration) of Pui et al. (1987) and
may include deposition in the straight sections (Section 6). The binomial
standard error
$\sqrt{DF(1-DF)/N}$ assumes statistically independent particle histories;
because the generator correlates particles (Section 5), it is a lower bound.

*Table 3. Deposition fraction (DF) from the Fluent DPM summaries of the size
classes whose summaries are archived
([`data/deposition_fraction.csv`](../data/deposition_fraction.csv)), and the
correlation of Pui et al. (1987) at the same Stokes number.*

| $d_p$ (µm) | $Stk$ | Tracked | Escaped | Deposited | DF (%) | Standard error (%) | Pui et al. (%) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 10ᵃ | 0.092 | 10,000 | 8,329 | 1,671 | 16.7 | 0.37 | 18.5 |
| 20 | 0.368 | 10,000 | 4,279 | 5,721 | 57.2 | 0.49 | 55.8 |
| 30 | 0.829 | 10,000 | 2,142 | 7,858 | 78.6 | 0.41 | 84.1 |
| 40 | 1.473 | 10,000 | 546 | 9,454 | 94.5 | 0.23 | 96.2ᵇ |
| 50 | 2.302 | 10,000 | 199 | 9,801 | 98.0 | 0.14 | 99.4ᵇ |

ᵃ Nominal diameter; see Section 6. ᵇ Beyond the measured range of the correlation.

![Deposition fraction](../figures/deposition-fraction.png)

*Figure 4. Deposition fraction against Stokes number (ρ<sub>p</sub> = 1000 kg m⁻³
assumed), with the correlation of Pui et al. (1987) for R<sub>0</sub> = 5.7 and
Re = 6000–10,000, drawn dashed beyond the measured range. Binomial 95 %
intervals are smaller than the markers. Generated by
[`tools/plot_deposition.py`](../tools/plot_deposition.py).*

The deposition fraction increases monotonically with $Stk$, from 16.7 % at
$Stk$ = 0.09, where particle inertia is small, to near-complete deposition at
$Stk > 1.5$. The correlation of Pui et al. (1987),
$DF = 1 - 10^{-0.963\,Stk}$, was fitted to measurements in 90° bends with
tube diameters of 5.03–8.51 mm, $R_0$ = 5.7 and $Re$ = 6000–10,000; the
largest measured Stokes number is about 1.35 (Sippola and Nazaroff, 2002,
Eq. 121, Table 8 and Fig. 6). At equal $Stk$ the computed deposition fractions
differ from the correlation by −5.5 to +1.4 percentage points. The present
case differs in curvature ratio (4), Reynolds number (4.1 × 10³) and inlet
condition (uniform velocity $5D$ upstream of the bend), its deposition
fraction may include the straight sections, and its two largest Stokes
numbers lie beyond the measured range. With $\rho_p$ = 1550 kg m⁻³, the default density of the material named
in the archived injection dialog, the correlation would give 27.1, 71.8 and
94.2 % at the three smallest diameters. Given these differences and the
departures in Table 2, the comparison is qualitative and is not a validation
of the model.

## 5. Random-number generator test

The function `gauss_rand` (lines 94–121) generates Box–Muller deviates from a
hash of the particle memory address and the particle time, without a stored
state. The program [`tools/rng_check.c`](../tools/rng_check.c) runs a
functionally identical copy (types adapted to run outside Fluent) at the time
step of the present runs.

*Table 4. Statistics of the UDF generator at Δt = 10⁻⁶ s (output of
`tools/rng_check.c`; the test number is given in brackets).*

| Statistic | Ideal | Result |
| --- | --- | --- |
| [1] Mean, variance (2 × 10⁵ samples) | 0, 1 | +0.001, 1.003 |
| [1] Autocorrelation, lags 1 to 5 | 0 | +0.22, +0.06, −0.11, −0.29, −0.23 |
| [2] Correlation of two particles at equal time, address gaps 8 B, 32 B, 448 B, 8 KiB | 0 | +0.85, −0.49, +0.13, −0.65 |
| [4] Normalised variance of the sum of 1000 successive deviates, 2000 particle addresses | 1 | median 0.48; 5–95 % range 0.20–2.7; maximum 60 |
| [5] Time average of $\eta^2$ in the Langevin update over 0.044 s ($\tau_L$ = 3.32 ms, $\tau_L^+$ = 10), 2000 particle addresses | 1 (ideal generator: median 0.94, 5–95 % range 0.51–1.70) | median 0.41; 5–95 % range 0.15–2.9; maximum 465 |

The mean and variance agree with those of a standard normal distribution.
Because the particle time advances in equal steps, successive deviates of one
particle are correlated, and deviates of two particles at the same time can be
strongly correlated. Driven by this generator, the Langevin update gives a
time-averaged variance of the modelled fluid velocity below its target for
most particles and far above it for a few. The random forcing, and with it
the effective turbulent diffusivity, therefore varies between particles, and
the standard errors in Table 3 underestimate the statistical uncertainty.

## 6. Limitations and further work

**Recorded settings.** The archived injection dialog for the smallest size
class ([`figures/setup/injection.png`](../figures/setup/injection.png))
specifies $d_p$ = 5 µm, the material "anthracite" and zero initial velocity,
whereas Table 3 assumes 10 µm and 1000 kg m⁻³ and the UDF header specifies
release at the gas velocity; at 5 µm, $Stk$ for this class would be 0.023.
Gravity, the inlet turbulence, the drag law, the DPM step limit, the solver
precision and the assignment of `dehbi_wall_bc` to the three wall zones
(inlet pipe, bend, outlet pipe) are not recorded in the archived dialogs, and
deposition was not recorded per wall zone.

**Implementation.** The departures in Table 2 are expected to matter most for
the smallest particles, whose deposition depends on turbulent dispersion
rather than on inertial impaction. Whether Fluent reuses particle memory
between particles, which would make particles share random sequences, has not
been tested.

**Further work.** A quantitative assessment of the model requires (i) a
per-particle random-number state, (ii) near-wall statistics evaluated at the
particle position throughout $y^+ < 100$, (iii) a consistent transformation of
the Langevin state between frames, (iv) a well-mixed tracer test, (v) the
confirmed injection settings, and (vi) a comparison with bend experiments at
matched curvature ratio and Reynolds number.

## 7. Conclusions

The CRW model of Dehbi (2008) was implemented as a compiled ANSYS Fluent UDF
and applied to a 90° pipe bend at $Re$ = 4.1 × 10³. The computed deposition
fraction increases from 16.7 % at $Stk$ = 0.09 to 98.0 % at $Stk$ = 2.3 and
differs from the correlation of Pui et al. (1987) by −5.5 to +1.4 percentage
points at equal Stokes number. The implementation departs from the published
model in three respects: the near-wall statistics and drift are applied only
in wall-adjacent cells, the Langevin state is not transformed between the
near-wall and bulk frames, and successive random deviates are correlated. The
results show that the model can be coupled to Fluent's discrete phase model
through body-force and boundary-condition UDFs; a quantitative test of the
model requires the changes listed in Section 6.

## Software

The Cursor code editor was used during UDF development to debug C syntax and
compilation errors related to the Fluent UDF macros. The PDF version of this
report is generated from the Markdown source with pandoc and Chromium.

## References

1. Dehbi, A. (2008). Turbulent particle dispersion in arbitrary wall-bounded
   geometries: A coupled CFD-Langevin-equation based approach.
   *International Journal of Multiphase Flow*, 34(9), 819–828.
   <https://doi.org/10.1016/j.ijmultiphaseflow.2008.03.001>
2. Kallio, G. A., & Reeks, M. W. (1989). A numerical simulation of particle
   deposition in turbulent boundary layers.
   *International Journal of Multiphase Flow*, 15(3), 433–446.
   <https://doi.org/10.1016/0301-9322(89)90012-8>
3. Morsi, S. A., & Alexander, A. J. (1972). An investigation of particle
   trajectories in two-phase flow systems. *Journal of Fluid Mechanics*,
   55(2), 193–208. <https://doi.org/10.1017/S0022112072001806>
4. Pui, D. Y. H., Romay-Novas, F., & Liu, B. Y. H. (1987). Experimental study
   of particle deposition in bends of circular cross section.
   *Aerosol Science and Technology*, 7(3), 301–315.
   <https://doi.org/10.1080/02786828708959166>
5. Sippola, M. R., & Nazaroff, W. W. (2002). *Particle deposition from
   turbulent flow: Review of published research and its applicability to
   ventilation ducts in commercial buildings* (Report LBNL-51432). Lawrence
   Berkeley National Laboratory. <https://doi.org/10.2172/815478>
6. Thomson, D. J. (1987). Criteria for the selection of stochastic models of
   particle trajectories in turbulent flow. *Journal of Fluid Mechanics*, 180,
   529–556. <https://doi.org/10.1017/S0022112087001940>
