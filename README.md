# CRW–Langevin Particle Dispersion UDF for ANSYS Fluent

[![Report](https://img.shields.io/badge/Technical_report-PDF-1F4E79?style=flat-square)](report/crw-dehbi2008-technical-report.pdf)
[![ANSYS Fluent](https://img.shields.io/badge/ANSYS_Fluent-2024_R2-FFB71B?style=flat-square)](udf/crw_dehbi2008_rev10.c)
[![Language](https://img.shields.io/badge/language-C-555555?style=flat-square)](udf/crw_dehbi2008_rev10.c)
[![License](https://img.shields.io/badge/license-MIT-2E7D32?style=flat-square)](LICENSE)

Implementation of the continuous random walk (CRW) model of
[Dehbi (2008)](https://doi.org/10.1016/j.ijmultiphaseflow.2008.03.001) as a
compiled ANSYS Fluent user-defined function, applied to turbulent aerosol
deposition in a 90° pipe bend.

| Computational domain | Bend cross-section mesh (5.3 × 10⁶ cells) |
| :---: | :---: |
| <img src="figures/geometry.png" width="380" alt="90-degree pipe bend with straight inlet and outlet sections"> | <img src="figures/mesh-refined-5.3M.png" width="380" alt="Structured hexahedral mesh of the bend cross-section"> |

**Author:** Sina Karbakhsh Ravari, Shiraz University ·
ORCID [0009-0003-6017-8628](https://orcid.org/0009-0003-6017-8628) ·
Supervisor: Dr. Ebrahim Ghahramani

**Technical report:** [Markdown](report/crw-dehbi2008-technical-report.md) ·
[PDF](report/crw-dehbi2008-technical-report.pdf)

## Overview

Discrete random walk models assume locally isotropic turbulence and
over-predict the deposition of low-inertia particles in wall-bounded flows.
The UDF replaces Fluent's built-in model with a normalised Langevin equation
for the fluid velocity seen by each particle, using anisotropic near-wall
statistics and a wall-normal drift correction, coupled to the discrete phase
model through a body force.

| Component | Implementation |
| --- | --- |
| Langevin equation | semi-implicit Euler–Maruyama scheme |
| Near-wall statistics ($y^+ < 100$) | DNS-based rms fits of Dehbi (2008) in a body-fitted frame |
| Drift correction | $A_2 = (\partial\sigma_2/\partial y)/(1+St)$, analytical derivative |
| Lagrangian time scale | Kallio and Reeks (1989) near the wall; $(2/C_0)\,k/\varepsilon$, $C_0 = 14$, in the bulk |
| Coupling | fluctuating drag as a DPM body force, Schiller–Naumann correction |
| Time step | $\min(10^{-6}\ \text{s},\ 0.1\tau_p)$ |

## Case and results

Air at 3.0 m/s ($Re \approx 4.1\times10^3$) in a pipe of diameter 0.02 m with
bend curvature ratio $R/r = 4$; Reynolds stress model with enhanced wall
treatment on a 5.3 × 10⁶-cell hexahedral mesh (ANSYS Fluent 2024 R2);
10,000 particles per diameter.

| $d_p$ (µm) | 10 | 20 | 30 | 40 | 50 |
| --- | --- | --- | --- | --- | --- |
| Flow Stokes number $Stk$ | 0.092 | 0.368 | 0.829 | 1.473 | 2.302 |
| Deposition fraction (%) | 16.7 | 57.2 | 78.6 | 94.5 | 98.0 |

<p align="center"><img src="figures/deposition-fraction.png" width="560" alt="Deposition fraction against particle diameter and flow Stokes number"></p>

The report documents the model equations with source-line references, a
statistical verification of the random-number generator, and the scope and
limitations of the present results.

## Repository structure

```
udf/       UDF source (crw_dehbi2008_rev10.c)
report/    technical report (Markdown and PDF)
data/      deposition fractions from the Fluent DPM summaries
figures/   geometry, mesh, wall y+, residuals, results; setup/ holds the Fluent dialogs
tools/     plotting script and random-number-generator test
```

## Usage

1. Solve the flow with the RSM and enhanced wall treatment.
2. Allocate 6 user-defined memory locations and 8 DPM user scalars.
3. Compile and load `udf/crw_dehbi2008_rev10.c`; run `init_wall_udm` on demand.
4. Disable the built-in turbulent dispersion. Hook `crw_langevin_force`
   (body force) and `crw_adaptive_dt` (time step), and assign `dehbi_wall_bc`
   to the walls.
5. Inject particles at the inlet with a DPM step limit above 10⁵. The
   deposition fraction is the number of "aborted" particles divided by the
   number tracked.

```
pip install -r tools/requirements.txt && python3 tools/plot_deposition.py
gcc -O2 -o rng_check tools/rng_check.c -lm && ./rng_check
```

## Citation and licence

Please cite this repository using [`CITATION.cff`](CITATION.cff). The code is
released under the [MIT License](LICENSE).
