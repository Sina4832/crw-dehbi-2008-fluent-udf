"""Plot the deposition fraction recorded in data/deposition_fraction.csv.

Usage (from the repository root):
    python tools/plot_deposition.py

Writes figures/deposition-fraction.png. Requires Python 3 and matplotlib.
Only values backed by a Fluent DPM summary are plotted; no benchmark data are drawn,
because digitized Dehbi (2008) and Pui et al. (1987) values are not yet in this repository.
"""
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

root = Path(__file__).resolve().parents[1]
rows = []
with open(root / "data" / "deposition_fraction.csv", newline="") as fh:
    for row in csv.DictReader(line for line in fh if not line.startswith("#")):
        rows.append(row)

d = [float(r["d_p_um"]) for r in rows]
df = [100 * float(r["deposition_fraction"]) for r in rows]
err = [100 * 1.96 * float(r["std_error"]) for r in rows]

fig, ax = plt.subplots(figsize=(5.5, 4.0), dpi=150)
ax.errorbar(d, df, yerr=err, fmt="o", color="#1f4e79", ecolor="#1f4e79",
            capsize=3, markersize=6, label="Present CRW UDF (Fluent DPM summary)")
ax.set_xlabel(r"Particle diameter $d_p$ ($\mu$m)")
ax.set_ylabel("Deposition fraction (%)")
ax.set_xlim(0, 55)
ax.set_ylim(0, 100)
ax.grid(True, alpha=0.3)
ax.legend(loc="lower right", frameon=False, fontsize=8)
ax.set_title("90° bend, 10,000 particles per diameter", fontsize=9)
fig.text(0.01, 0.01, "Error bars: 95% binomial sampling interval only.", fontsize=6.5, color="#555")
fig.tight_layout()
fig.savefig(root / "figures" / "deposition-fraction.png")
print("wrote figures/deposition-fraction.png")
