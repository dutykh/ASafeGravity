#!/usr/bin/env python3
# Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
"""Figure of the README: the least damped oscillating mode of the Hayward black hole against gamma.

The frequencies are read from the reports of qnmdata/Paper3HaywardBH/ (the first line of the list
of oscillatory modes of each report, M omega), from the Schwarzschild black hole (gamma = 0) to the
extremal one (gamma = 32/27). The folder names give gamma: g0, g1, g4_27 = 4/27, g59_50 = 59/50;
the reports of the extremal black hole are those of the folder extreme/. The second report of
s = 1, l = 1, gamma = 1, computed at higher resolution, gives the same mode and is not used. From
this folder,
    python3 plot_readme_figure.py
writes hayward-fundamental-modes.png (Python 3 with NumPy and Matplotlib). The style is the stamped
copy manuscript_style.py beside it.
"""
import re
import sys
from fractions import Fraction
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import manuscript_style as hs  # noqa: E402

DATA = HERE.parent.parent / "qnmdata" / "Paper3HaywardBH"
GAMMA_EXTREMAL = Fraction(32, 27)
SERIES = [(0, 0), (0, 1), (1, 1), (2, 2)]
NAMES = {0: "scalar", 1: "electromagnetic", 2: "axial model"}


def least_damped(report: Path) -> complex:
    """The first oscillatory mode of a report, written as 0.1104549-0.1048957*I."""
    lines = report.read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("Oscillatory modes"))
    match = re.search(r"Mode #1\s*:\s*([0-9.]+)-([0-9.]+)\*I", lines[start + 1])
    return complex(float(match.group(1)), -float(match.group(2)))


def gamma_of(folder: str) -> Fraction:
    numerator, _, denominator = folder[1:].partition("_")
    return Fraction(int(numerator), int(denominator or 1))


def series(s: int, ell: int):
    points = {}
    for report in DATA.glob(f"*extreme/s{s}/L{ell}/g*/*.txt"):
        if "hires" in report.name:
            continue
        points[gamma_of(report.parent.name)] = least_damped(report)
    for report in DATA.glob(f"extreme/s{s}/L{ell}/*.txt"):
        points[GAMMA_EXTREMAL] = least_damped(report)
    gammas = sorted(points)
    return np.array([float(g) for g in gammas]), np.array([points[g] for g in gammas])


hs.configure(family="libertine")
plt.rcParams["savefig.dpi"] = 260
plt.rcParams["figure.facecolor"] = hs.PAPER
plt.rcParams["savefig.facecolor"] = hs.PAPER

fig, (left, right) = plt.subplots(1, 2, figsize=(hs.WIDTH, 2.75))
for k, (s, ell) in enumerate(SERIES):
    gamma, omega = series(s, ell)
    reference = omega[0]
    assert gamma[0] == 0.0, "every series starts at the Schwarzschild black hole"
    print(f"s = {s}, l = {ell}: gamma = {[round(g, 5) for g in gamma]}")
    print(f"    M omega at gamma = 0: {reference.real:.7f} - {-reference.imag:.7f} i, "
          f"at gamma = 32/27: {omega[-1].real:.7f} - {-omega[-1].imag:.7f} i")
    style = hs.series(k)
    label = rf"$s = {s}$, $\ell = {ell}$ ({NAMES[s]})"
    left.plot(gamma, 100 * (omega.real / reference.real - 1), lw=1.0, markersize=3.6, label=label, **style)
    right.plot(gamma, 100 * (omega.imag / reference.imag - 1), lw=1.0, markersize=3.6, **style)
for axis in (left, right):
    axis.axvline(float(GAMMA_EXTREMAL), color=hs.SLATE, lw=0.6, linestyle=(0, (1.6, 1.6)), zorder=0)
    axis.axhline(0.0, color=hs.MIST, lw=0.6, zorder=0)
    axis.set_xlim(-0.03, 1.25)
    axis.set_xlabel(r"$\gamma$ (extremal at $32/27$, dotted)")
left.set_ylabel(r"change of $\mathrm{Re}\,M\omega$ (per cent)")
right.set_ylabel(r"change of $-\mathrm{Im}\,M\omega$ (per cent)")
hs.panel(left, r"(a) oscillation frequency")
hs.panel(right, r"(b) damping rate")
hs.legend_outside(fig, where="lower", ncol=2)

out = HERE / "hayward-fundamental-modes.png"
hs.save(fig, out)
# The background is opaque: without its alpha channel, and compressed a second time without
# loss, the file is smaller.
Image.open(out).convert("RGB").save(out, optimize=True)
print("written", out)
