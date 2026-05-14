Drift-Corrected CRW UDF for ANSYS Fluent

Reproduction of Dehbi (2008) Aerosol Deposition Benchmark

Overview

This repository contains a custom User-Defined Function (UDF) written in C for
ANSYS Fluent to reproduce the drift-corrected Continuous Random Walk (CRW) model proposed by Dehbi (2008) for turbulent aerosol deposition in a 90° pipe bend.

The project was developed as a CFD validation study focused on overcoming the limitations of Fluent’s standard Discrete Random Walk (DRW) particle tracking model, particularly the artificial near-wall accumulation of low-inertia particles.

The implementation integrates a generalized Langevin stochastic differential equation with Reynolds-stress-based drift corrections inside Fluent’s Discrete Phase Model (DPM).

---

Objectives

- Reproduce the benchmark deposition curve presented in Dehbi (2008)
- Validate aerosol deposition behavior against the experimental work of Pui et al. (1987)
- Implement a custom drift-corrected CRW model using Fluent UDFs
- Improve prediction accuracy for low-inertia particle deposition
- Investigate mesh sensitivity in the viscous sublayer region

---

Numerical Methodology

Geometry

- 90° pipe bend
- Pipe diameter: 2 cm
- Curvature ratio: R/r = 4
- Straight inlet and outlet extensions included

Continuous Phase

- Solver: ANSYS Fluent 2024 R2
- Pressure-based steady-state solver
- Turbulence model: Reynolds Stress Model (RSM)
- Enhanced Wall Treatment enabled
- Second-order spatial discretization

Discrete Phase

- Lagrangian particle tracking
- One-way coupling
- Particle diameters:
  - 10 μm
  - 15 μm
  - 20 μm
  - 30 μm
  - 40 μm
  - 50 μm
- 10,000 particles injected per simulation

Near-Wall Resolution

A highly refined hexahedral mesh was employed to fully resolve the viscous sublayer:

- y+ ≈ 0.95
- 30 inflation layers
- Growth rate: 1.1
- Final mesh size: ~5.3 million cells

---

UDF Implementation

The repository includes a custom C-language UDF implementing:

- Continuous Random Walk (CRW) turbulent dispersion
- Generalized Langevin equation
- Reynolds stress gradient drift correction
- Particle fluctuation velocity update
- Custom stochastic tracking for Fluent DPM

The implementation replaces Fluent’s default DRW approach to maintain physically realistic particle concentration behavior near walls.

---

Validation

The final deposition fractions showed strong agreement with:

- Experimental data from Pui et al. (1987)
- Numerical benchmark results from Dehbi (2008)

The study successfully reproduced the characteristic deposition curve reported in Figure 8 of Dehbi’s publication.

---

Repository Structure

.
├── udf/
│   └── crw_dehbi2008_rev10.c
│
├── report/
│   └── report-dehbi.pdf
│
├── figures/
│   ├── validation-curve.png
│   ├── mesh.png
│   └── yplus.png
│
└── README.md

---

Compilation

The UDF was compiled using:

- Microsoft Visual Studio x64 Native Tools
- Fluent compiled UDF workflow
- Successfully loaded as "libudf.dll"

---

Software

- ANSYS Fluent 2024 R2
- Microsoft Visual Studio
- C Language (Fluent UDF API)

---

References

Dehbi, A. (2008)

A generalized Brownian diffusion model for particle transport in turbulent flows.

Pui, D. Y. H., Romay-Novas, F., & Liu, B. Y. H. (1987)

Experimental study of particle deposition in bends of circular cross section.

---

Notes

This repository is intended for academic and research purposes related to:

- CFD
- Aerosol transport
- Turbulent particle deposition
- ANSYS Fluent UDF development
- Lagrangian stochastic modeling

---

Author

Sina Karbakhsh
