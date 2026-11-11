"""Plot the deposition fraction in data/deposition_fraction.csv against the
Stokes number, with the turbulent-bend correlation of Pui et al. (1987).

Usage (from the repository root):
    python3 tools/plot_deposition.py

Writes figures/deposition-fraction.png. Requires Python 3 and matplotlib.

Correlation: DF = 1 - 10**(-0.963 Stk), fitted to measurements in 90-degree
bends with curvature ratio 5.7 and Re = 6000-10,000 (Pui et al., 1987; as
reproduced by Sippola and Nazaroff, 2002, Eq. 121 and Table 8). It is drawn
solid up to Stk = 1.35, the largest Stokes number of the measurements, and
dashed beyond. It is shown for reference only: the present case has a
curvature ratio of 4 and Re = 4.1e3 (see the technical report, Section 4).
"""
import csv
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("Agg")

STK_MEASURED_MAX = 1.35
STK_PLOT_MAX = 2.6

root = Path(__file__).resolve().parents[1]
with open(root / "data" / "deposition_fraction.csv", newline="") as fh:
    rows = list(csv.DictReader(ln for ln in fh if not ln.startswith("#")))

d_p = [float(r["d_p_um"]) for r in rows]
stk = [float(r["stokes_number"]) for r in rows]
dep = [100.0 * float(r["deposition_fraction"]) for r in rows]


def pui(s):
    """Deposition fraction (%) from the Pui et al. (1987) correlation."""
    return 100.0 * (1.0 - 10.0 ** (-0.963 * s))


n = 200
s_meas = [STK_MEASURED_MAX * i / n for i in range(n + 1)]
s_extr = [STK_MEASURED_MAX + (STK_PLOT_MAX - STK_MEASURED_MAX) * i / n
          for i in range(n + 1)]

matplotlib.rcdefaults()  # independent of any local matplotlibrc
plt.rcParams.update({"font.size": 11, "axes.labelsize": 11,
                     "legend.fontsize": 9.5})
fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=150)
ax.plot(s_meas, [pui(s) for s in s_meas], color="0.45", lw=1.6,
        label="Pui et al. (1987) correlation")
ax.plot(s_extr, [pui(s) for s in s_extr], color="0.45", lw=1.6, ls="--",
        label="Pui et al. (1987), beyond measured range")
ax.plot(stk, dep, "o", color="#1f4e79", ms=7, zorder=3,
        label="Present CRW UDF")
for s, y, d in zip(stk, dep, d_p, strict=True):
    ax.annotate(f"{d:.0f} µm", (s, y), textcoords="offset points",
                xytext=(7, -13), fontsize=9, color="#1f4e79")
ax.set_xlabel(r"Stokes number $Stk = \rho_p d_p^2 U / (9 \mu D)$")
ax.set_ylabel("Deposition fraction (%)")
ax.set_xlim(0.0, STK_PLOT_MAX)
ax.set_ylim(0.0, 100.0)
ax.grid(True, alpha=0.3)
ax.legend(loc="lower right", frameon=False)
fig.tight_layout()
fig.savefig(root / "figures" / "deposition-fraction.png")
print("wrote figures/deposition-fraction.png")
