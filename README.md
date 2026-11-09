# CRW–Langevin particle dispersion UDF for ANSYS Fluent (Dehbi, 2008)

A compiled ANSYS Fluent user-defined function (UDF) that replaces Fluent's
discrete random walk with the continuous random walk (CRW) model of
[Dehbi (2008)](https://doi.org/10.1016/j.ijmultiphaseflow.2008.03.001), applied
to turbulent aerosol deposition in a 90° pipe bend.

**Author:** Sina Karbakhsh Ravari, M.Sc. candidate in Aerospace Engineering
(Aerodynamics), Shiraz University · ORCID
[0009-0003-6017-8628](https://orcid.org/0009-0003-6017-8628)
**Project supervisor:** Dr. Ebrahim Ghahramani
**Status:** working implementation with deposition results; benchmark
comparison not yet made; three implementation issues documented and scheduled
for correction.

## Summary

- **What was built:** a ~800-line C UDF that integrates a normalised Langevin
  equation for the fluid velocity seen by each particle, with anisotropic
  near-wall statistics and a wall-normal drift term, coupled to Fluent's
  discrete phase model through a body force.
- **Result:** on a Reynolds-stress-model flow field, the deposition fraction
  rises from 16.7 % at 10 µm (Stokes number 0.09) to 98.0 % at 50 µm
  (Stokes number 2.3).
- **Limits:** the case differs from the benchmark experiment of Pui et al.
  (1987), and a line-by-line code review found three issues that affect the
  physics. Both are documented below and in the report.

![Deposition fraction](figures/deposition-fraction.png)

## Why this model

Discrete random walk (DRW) models sample velocity fluctuations from locally
isotropic turbulence. In inhomogeneous near-wall turbulence this produces a
spurious drift of low-inertia particles towards the wall, a violation of the
well-mixed condition, so the deposition of small particles is over-predicted.
Dehbi's model integrates a normalised Langevin equation with anisotropic
near-wall statistics and a wall-normal drift term that restores the well-mixed
condition.

## What the UDF implements

| Component | Implementation | Function |
| --- | --- | --- |
| Normalised Langevin equation | semi-implicit Euler–Maruyama, $\eta^{n+1} = (\eta^n + \sqrt{2\Delta t/\tau_L}\,\xi + A\Delta t)/(1+\Delta t/\tau_L)$ | `advance_langevin` |
| Near-wall statistics, $y^+<100$ (currently first cell layer only, see issue 1) | DNS-based fits for $\sigma_1,\sigma_2,\sigma_3$ from Dehbi (2008) in a body-fitted frame | `sigma_plus`, `build_bfcs` |
| Wall-normal drift | $A_2 = (\partial\sigma_2/\partial y)/(1+St)$ with $St=\tau_p/\tau_L$ the local particle Stokes number; analytical derivative | `dsigma2p_dyp` |
| Lagrangian time scale | Kallio and Reeks (1989) near the wall; $(2/C_0)\,k/\varepsilon$, $C_0=14$, in the bulk | `tau_L` |
| Coupling to Fluent | drag of the velocity fluctuation as a DPM body force, Schiller–Naumann drag | `crw_langevin_force` |
| Particle time step | $\min(10^{-6}\ \text{s},\ 0.1\tau_p)$ | `crw_adaptive_dt` |
| Deposition | particle stopped at first wall contact | `dehbi_wall_bc` |

The equations as implemented, with source-line references, are in Section 5
of the [technical report](report/crw-dehbi2008-technical-report.md)
([PDF](report/crw-dehbi2008-technical-report.pdf)).

## Case and results

Air at 3.0 m/s ($Re \approx 4.1\times10^3$) in a pipe of diameter 0.02 m with
bend curvature ratio $R/r = 4$, solved with the Reynolds stress model and
enhanced wall treatment on a hexahedral mesh of about 5.3 million cells with
wall $y^+ \approx 0.95$ (ANSYS Fluent 2024 R2). 10,000 particles were tracked
per diameter.

| $d_p$ (µm) | 10 | 20 | 30 | 40 | 50 |
| --- | --- | --- | --- | --- | --- |
| Flow Stokes number $Stk$ ¹ | 0.092 | 0.368 | 0.829 | 1.473 | 2.302 |
| Deposition fraction (%) | 16.7 | 57.2 | 78.6 | 94.5 | 98.0 |

¹ $Stk = \rho_p d_p^2 U/(9\mu D)$ with particle density 1000 kg/m³ as stated
in the original report. The density, the 10 µm injection diameter and the
gravity setting are being confirmed from the case file (report, Section 4).
The deposition fraction counts particles stopped on all walls, including the
straight sections.

Values come from the Fluent DPM summaries
([`data/deposition_fraction.csv`](data/deposition_fraction.csv)). The case
differs from the turbulent experiment of Pui et al. (1987), reported as
D = 8.51 mm, curvature ratio 5.6 and $Re$ = 10,000, so a benchmark comparison
must be made at matched Stokes number. That comparison has not yet been made,
and no agreement is claimed here (report, Sections 2 and 6.3).

## Known issues in the current UDF

A line-by-line review (report, Section 7) found three issues in
`crw_dehbi2008_rev10.c`, the version that produced the results above:

1. **Near-wall branch only in wall-adjacent cells.** Wall data are stored only
   for cells that touch a wall, so cells between $y^+\approx1$ and $y^+=100$
   use the isotropic bulk branch without drift. Where the branch is active,
   the statistics are evaluated at the cell centre rather than at the
   particle's own wall distance.
2. **Correlated random numbers.** The Gaussian generator hashes the particle
   address and time instead of keeping a state.
   [`tools/rng_check.c`](tools/rng_check.c) shows that, over 2000 simulated
   particle addresses, the variance of 1000 accumulated deviates has a median
   of 0.48 times the correct value and ranges from 0.08 to 60 times it.
3. **No frame rotation between branches.** The Langevin state is stored in the
   body-fitted frame near the wall and in Cartesian components in the bulk,
   without rotation when a particle changes branch.

In addition, the bulk drift term proportional to $\nabla k$ is set to zero,
because reading the gradient of $k$ under the RSM crashed Fluent. Revision 3 of
the UDF will correct these issues, after which all diameters will be re-run
and compared with the benchmark.

## Repository contents

```
udf/crw_dehbi2008_rev10.c                 UDF source, unchanged from the version used for the results
report/crw-dehbi2008-technical-report.md  technical report, revision 2 (source)
report/crw-dehbi2008-technical-report.pdf technical report, revision 2 (PDF of the same text)
report/original/report-dehbi-v1-2026-03.pdf  revision 1 (March 2026), superseded; some claims withdrawn in revision 2
data/deposition_fraction.csv              deposition fractions from the Fluent DPM summaries
figures/                                  geometry, mesh, y+, residuals, results; setup/ holds Fluent dialogs
tools/plot_deposition.py                  regenerates figures/deposition-fraction.png
tools/rng_check.c                         statistical test of the UDF random-number generator
```

Fluent case and data files are not included because of their size.

## How to use the UDF

1. Solve the flow with the RSM and enhanced wall treatment, second-order
   discretisation.
2. Allocate 6 user-defined memory locations, and set the discrete phase model
   to 8 user scalars per particle.
3. Compile `udf/crw_dehbi2008_rev10.c` (User-Defined → Functions → Compiled)
   and load `libudf`. The original build used the
   *x64 Native Tools Command Prompt for Visual Studio* on Windows.
4. After the flow has converged, run Execute On Demand → `init_wall_udm`.
5. Disable the built-in turbulent dispersion model. Hook `crw_langevin_force`
   (DPM body force) and `crw_adaptive_dt` (DPM time step), set the DPM wall
   condition of the pipe walls to user-defined with `dehbi_wall_bc`, and set
   the inlet and outlet to escape.
6. Inject particles uniformly over the inlet. Set the DPM maximum number of
   steps above about 10⁵, because each particle needs about 9 × 10⁴ steps of
   10⁻⁶ s to cross the domain. The deposition fraction is the number of
   "aborted" particles divided by the number tracked.

To regenerate the results figure:
`pip install -r tools/requirements.txt && python3 tools/plot_deposition.py`.
To run the generator test:
`gcc -O2 -o rng_check tools/rng_check.c -lm && ./rng_check`.

## Use of AI tools

The Cursor editor was used while developing the UDF to debug C syntax,
memory-access and compilation errors, and to tidy comments. The physical model
and equations follow Dehbi (2008). Revision 2 of the report, this README, the
two scripts in `tools/` and the code review were drafted with the help of an AI
assistant (Claude) and checked by the author against the source code and the
original Fluent output.

## References

- Dehbi, A. (2008). Turbulent particle dispersion in arbitrary wall-bounded
  geometries: A coupled CFD-Langevin-equation based approach.
  *Int. J. Multiphase Flow* 34(9), 819–828.
  [doi:10.1016/j.ijmultiphaseflow.2008.03.001](https://doi.org/10.1016/j.ijmultiphaseflow.2008.03.001)
- Pui, D. Y. H., Romay-Novas, F., Liu, B. Y. H. (1987). Experimental study of
  particle deposition in bends of circular cross section.
  *Aerosol Sci. Technol.* 7(3), 301–315.
- Kallio, G. A., Reeks, M. W. (1989). A numerical simulation of particle
  deposition in turbulent boundary layers. *Int. J. Multiphase Flow* 15(3),
  433–446.

## License and citation

The source code in `udf/` and `tools/` is released under the
[MIT License](LICENSE). The report, data and figures are the author's work;
please cite this repository if you use them (see [`CITATION.cff`](CITATION.cff)).
