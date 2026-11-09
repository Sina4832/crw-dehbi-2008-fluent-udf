# CRW–Langevin particle dispersion UDF for ANSYS Fluent (Dehbi, 2008)

A compiled ANSYS Fluent user-defined function (UDF) that replaces Fluent's
discrete random walk with the continuous random walk (CRW) model of
[Dehbi (2008)](https://doi.org/10.1016/j.ijmultiphaseflow.2008.03.001), applied
to turbulent aerosol deposition in a 90° pipe bend.

**Author:** Sina Karbakhsh Ravari, M.Sc. candidate in Aerospace Engineering
(Aerodynamics), Shiraz University · ORCID
[0009-0003-6017-8628](https://orcid.org/0009-0003-6017-8628)
**Supervisor:** Dr. Ebrahim Ghahramani
**Status:** working implementation with preliminary validation; two known
implementation issues are documented below.

![Deposition fraction](figures/deposition-fraction.png)

## Why this model

RANS eddy-interaction models assume locally isotropic turbulence. In
wall-bounded flows this violates the well-mixed condition, so small particles
accumulate artificially near walls and deposition is over-predicted. Dehbi's
model integrates a normalised Langevin equation with anisotropic near-wall
velocity statistics and a wall-normal drift term that restores the well-mixed
condition.

## What the UDF implements

| Component | Implementation | Code |
| --- | --- | --- |
| Normalised Langevin equation | semi-implicit Euler–Maruyama, $\eta^{n+1} = (\eta^n + \sqrt{2\Delta t/\tau_L}\,\xi + A\Delta t)/(1+\Delta t/\tau_L)$ | `advance_langevin` |
| Near-wall statistics ($y^+<100$) | DNS fits for $\sigma_1,\sigma_2,\sigma_3$ (Dreeben and Pope, 1997) in a body-fitted frame | `sigma_plus`, `build_bfcs` |
| Wall-normal drift | $A_2 = (\partial\sigma_2/\partial y)/(1+St)$, analytical derivative | `dsigma2p_dyp` |
| Lagrangian time scale | Kallio and Reeks (1989) near the wall; $(2/C_0)\,k/\varepsilon$, $C_0=14$, in the bulk | `tau_L` |
| Coupling to Fluent | fluctuation drag as a DPM body force with Schiller–Naumann drag | `crw_langevin_force` |
| Particle time step | $\min(10^{-6}\ \text{s},\ 0.1\tau_p)$ | `crw_adaptive_dt` |
| Deposition | trap at first wall contact | `dehbi_wall_bc` |

The full derivation, with line references, is in Section 5 of the
[technical report](report/crw-dehbi2008-technical-report.md).

## Case and results

Air at 3.0 m/s ($Re \approx 4.1\times10^3$) in a pipe of diameter 0.02 m with
$R/r = 4$, solved with the Reynolds stress model and enhanced wall treatment
on a 5.3-million-cell hexahedral mesh with wall $y^+ \approx 0.95$.
10,000 particles of density 1000 kg/m³ were tracked per diameter
(ANSYS Fluent 2024 R2).

| $d_p$ (µm) | 10 | 20 | 30 | 40 | 50 |
| --- | --- | --- | --- | --- | --- |
| Stokes number | 0.092 | 0.368 | 0.829 | 1.473 | 2.302 |
| Deposition fraction (%) | 16.7 | 57.2 | 78.6 | 94.5 | 98.0 |

Values come from the Fluent DPM summaries
([`data/deposition_fraction.csv`](data/deposition_fraction.csv)). An
independent overlay of the Dehbi (2008) and Pui et al. (1987) benchmark data
has not yet been made, so agreement with the benchmark is not claimed here;
see Section 6.3 of the report.

## Known issues in the current UDF

A line-by-line review (report, Section 7) found two issues in
`crw_dehbi2008_rev10.c`, the version that produced the results above:

1. **Near-wall branch only in wall-adjacent cells.** Wall data are stored only
   for cells touching a wall, so cells between $y^+\approx1$ and $y^+=100$ use
   the isotropic bulk branch without drift.
2. **Correlated random numbers.** The Gaussian generator hashes the particle
   address and time instead of keeping a state. [`tools/rng_check.c`](tools/rng_check.c)
   shows lag-1 to lag-5 autocorrelations between −0.29 and +0.22 and a
   correlation of +0.13 between neighbouring particles.

The bulk $\nabla k$ drift term is also set to zero (FIX-17). Corrections are
planned for a revision 3 of the UDF, followed by a re-run of all diameters.

## Repository contents

```
udf/crw_dehbi2008_rev10.c      UDF source (unchanged from the version used for the results)
report/
  crw-dehbi2008-technical-report.md   technical report, revision 2
  original/report-dehbi-v1-2026-03.pdf  revision 1 (March 2026)
data/deposition_fraction.csv   deposition fractions from the Fluent DPM summaries
figures/                       geometry, mesh, y+, residuals, results; setup/ holds Fluent dialogs
tools/plot_deposition.py       regenerates figures/deposition-fraction.png
tools/rng_check.c              statistical test of the UDF random-number generator
```

Fluent case and data files are not included because of their size.

## How to use the UDF

1. Solve the flow with the RSM and enhanced wall treatment, second-order
   discretisation.
2. Allocate 6 user-defined memory locations and set the DPM to 8 user scalars
   per particle.
3. Compile `udf/crw_dehbi2008_rev10.c` (Define → User-Defined → Functions →
   Compiled) and load `libudf`. On Windows, start Fluent from the
   *x64 Native Tools Command Prompt for Visual Studio*.
4. After the flow has converged, run Execute On Demand → `init_wall_udm`.
5. Disable the built-in turbulent dispersion model and hook
   `crw_langevin_force` (DPM body force), `crw_adaptive_dt` (DPM time step)
   and `dehbi_wall_bc` (wall boundary condition of the pipe walls).
6. Inject particles uniformly over the inlet, released at the gas velocity,
   and compute the deposition fraction as trapped / tracked.

To regenerate the results figure: `python tools/plot_deposition.py`
(requires matplotlib). To run the generator test:
`gcc -O2 -o rng_check tools/rng_check.c -lm && ./rng_check`.

## Use of software tools

The Cursor editor was used while developing the UDF to debug C syntax,
memory-access and compilation errors and to tidy comments. The physical model
and equations follow Dehbi (2008).

## References

- Dehbi, A. (2008). Turbulent particle dispersion in arbitrary wall-bounded
  geometries: A coupled CFD-Langevin-equation based approach.
  *Int. J. Multiphase Flow* 34(9), 819–828.
  [doi:10.1016/j.ijmultiphaseflow.2008.03.001](https://doi.org/10.1016/j.ijmultiphaseflow.2008.03.001)
- Pui, D. Y. H., Romay-Novas, F., Liu, B. Y. H. (1987). Experimental study of
  particle deposition in bends of circular cross section.
  *Aerosol Sci. Technol.* 7(3), 301–315.
- Dreeben, T. D., Pope, S. B. (1997). Probability density function and
  Reynolds-stress modeling of near-wall turbulent flows.
  *Phys. Fluids* 9(1), 154–163.
- Kallio, G. A., Reeks, M. W. (1989). A numerical simulation of particle
  deposition in turbulent boundary layers. *Int. J. Multiphase Flow* 15(3),
  433–446.

## License and citation

The source code (`udf/`, `tools/`) is released under the [MIT License](LICENSE).
The report, data and figures are the author's work; please cite this
repository if you use them (see [`CITATION.cff`](CITATION.cff)).
