<!-- Author: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE) -->

# ASafeGravity: quasinormal modes of Planck stars, Bonanno–Reuter and Hayward black holes

![The least damped oscillating modes of the Hayward black hole from Schwarzschild to extremality](docs/figures/hayward-fundamental-modes.png)

The least damped oscillating mode of the Hayward black hole, read from the reports of `qnmdata/Paper3HaywardBH/`, as
the deformation parameter γ grows from the Schwarzschild black hole (γ = 0) to the extremal one (γ = 32/27, dotted
line), for the scalar field with ℓ = 0 and ℓ = 1, the electromagnetic field with ℓ = 1 and the axial model with
ℓ = 2. Panel (a) shows the relative change of the oscillation frequency Re Mω and panel (b) that of the damping rate
−Im Mω, both measured from their Schwarzschild values (0.110455 − 0.104896i, 0.292936 − 0.097660i,
0.248263 − 0.092488i and 0.373672 − 0.088962i). The deformation lowers the damping of every mode, by 16 to 22 per cent
at extremality, and raises the frequency, by up to 9.5 per cent for the axial quadrupole, except for the scalar
monopole, whose frequency rises up to γ ≈ 0.74 and then falls back almost to its Schwarzschild value.

ASafeGravity holds the codes and the computed spectra of a series of three articles on the quasinormal modes of black
holes whose metric is improved by a running Newton coupling, as in asymptotically safe gravity: Planck stars (Paper 1),
Bonanno–Reuter black holes (Paper 2) and Hayward black holes (Paper 3). For massless scalar (s = 0), electromagnetic
(s = 1) and gravitational (s = 2) perturbations, the radial equation is compactified onto the interval −1 ⩽ y ⩽ 1, its
quasinormal boundary conditions are built into a factorisation of the solution, and a Chebyshev collocation in
multiprecision arithmetic gives a quadratic eigenvalue problem whose converged eigenvalues are the quasinormal
frequencies. The repository holds the Maple worksheets of the derivations and the Maple procedures that assemble the
matrices for the three papers, and 186 reports of computed spectra. For Paper 3 it also holds the MATLAB script that
solves the eigenvalue problem with the Advanpix toolbox and writes the reports, and the Python driver that runs Maple
and MATLAB for every set of parameters; for Papers 1 and 2 the solution of the eigenvalue problem, done in MATLAB with
the Advanpix toolbox as well, is not part of it. The articles, by Davide Batic, Denys Dutykh and Fabio Scardigli
(Papers 1 and 2) and by Davide Batic and Denys Dutykh (Paper 3), are cited below.

## The three geometries

The three black holes are static and spherically symmetric, ds² = −f dt² + dr²/f + r² dΩ². Papers 1 and 2 use the
Newton coupling that runs with the radius, in Planck units (r in Planck lengths, M the dimensionless geometric mass):

```text
f(r) = 1 − 2M G(r)/r,        G(r) = r³ / ( r³ + α (r + γM) ),

Paper 1, Planck stars:              α = −41/(10π),    γ = 1/2 and 9/2,    M = 1
Paper 2, Bonanno–Reuter black holes: α = 118/(15π),    γ = 9/2
```

In Paper 1 the negative α matches the effective Newtonian potential of the one-loop effective field theory; the
coupling then does not reach the ultraviolet fixed point of asymptotically safe gravity, and the article reads the
geometry as an effective metric of scale-dependent gravity. In Paper 2 the positive α of Bonanno and Reuter gives a
regular black hole, with two horizons above an extremal mass and one double horizon at that mass (M = 3.50274 for
γ = 9/2 in the reports). Paper 3 uses the Hayward metric, in the rescaled radius ρ = r/M:

```text
f(ρ) = 1 − 2ρ² / (ρ³ + γ),        horizons:  ρ_h³ − 2ρ_h² + γ = 0,
```

with two horizons for 0 < γ < 32/27, a double one at ρ = 4/3 for γ = 32/27, and none beyond. In all three papers each
perturbation obeys, with the time dependence e^(−iωt) and spherical harmonics of degree ℓ ⩾ s,

```text
f d/dr( f dψ/dr ) + [ ω² − U(r) ] ψ = 0,        U(r) = f(r) [ ε f′(r)/r + ℓ(ℓ + 1)/r² ],        ε = 1 − s²,
```

with ψ ingoing at the event horizon and outgoing at infinity; for s = 2 the Hayward article defines its spin-2 equation
by this potential, an axial model rather than a consequence of a perturbed action. After the asymptotic behaviour at
both ends has been factored out, the regular part of ψ is expanded in N Chebyshev polynomials and collocated at the N
roots y_i = cos((2i − 1)π/(2N)), which gives

```text
( M0 + i ω M1 + ω² M2 ) a = 0
```

with three real N × N matrices; the procedures export `M1` without the factor i, which multiplies it when the problem is
solved.

## Quick start

The files `matrixassembler.mpl` of `matlab/` are Maple procedures, despite the name of the folder: they assemble the
three matrices for MATLAB. Each writes `data/M0_N.mat`, `data/M1_N.mat` and `data/M2_N.mat`, in the ASCII format of
MATLAB, into the folder of the path it is given, which must exist. To recompute the report of Paper 1 for s = 0, ℓ = 0, γ = 1/2 and M = 1,
whose header lists the resolutions 200, 190 and 180, from the folder of the repository:

```bash
mkdir -p data
```

```maple
read "matlab/Paper1NegAlpha/matrixassembler.mpl":
MatrixAssembler(200, 200, 0, 0, 0.5, 1.0, "/absolute/path/to/ASafeGravity"):
MatrixAssembler(200, 190, 0, 0, 0.5, 1.0, "/absolute/path/to/ASafeGravity"):
MatrixAssembler(200, 180, 0, 0, 0.5, 1.0, "/absolute/path/to/ASafeGravity"):
```

The arguments of the non-extremal procedures (Papers 1 and 2) are the number of digits `d`, the number of Chebyshev
modes `N` (the size of the matrices), the spin `s`, the multipole `L` (that is ℓ), the parameter `g` (that is γ), the
mass `M` and the path `p`; they find the event horizon by `fsolve`, from r = 2.6 in Paper 1 and r = 9.03 in Paper 2,
and print it. The extremal procedure of Paper 2, `matlab/Paper2PosAlpha/extreme/matrixassembler.mpl`, takes only
`(d, N, s, L, p)`: it sets γ = 9/2 and solves for the extremal mass and the double horizon (searched between 4 and 5),
which it prints. The eigenvalue problem is then solved in MATLAB with the Advanpix Multiprecision Computing Toolbox, as
in the articles:

```matlab
addpath('/path/to/advanpix');                     % Multiprecision Computing Toolbox
for N = [200, 190, 180]
    mp.Digits(N);
    M0 = mp.read(sprintf('data/M0_%d.mat', N));
    M1 = mp.read(sprintf('data/M1_%d.mat', N));
    M2 = mp.read(sprintf('data/M2_%d.mat', N));
    omega = polyeig(M0, mp('1i')*M1, M2);         % (M0 + i omega M1 + omega^2 M2) a = 0
    plot(real(omega), imag(omega), '.'), hold on
end
```

The quasinormal frequencies are the eigenvalues common to the three resolutions, which the reports list. For Paper 3,
the scripts of the next section do all of this, from the assembly of the matrices to the report.

## The codes of Paper 3

For Paper 3 the repository holds the whole chain that wrote the 53 reports of `qnmdata/Paper3HaywardBH/`: Maple
assembles the three matrices at each resolution, MATLAB solves the quadratic eigenvalue problem in multiprecision
arithmetic and writes the report, and a Python driver runs the two in turn for every set of parameters. The folder
`matlab/Paper3HaywardBH/` holds one copy of each file; the non-extremal and the near-extremal black holes share those of
`nonextreme/`.

| File in `matlab/Paper3HaywardBH/` | Content |
|---|---|
| `nonextreme/matrixassembler.mpl` | Maple procedure `MatrixAssembler(d, N, s, L, g, p)` for 0 ⩽ γ < 32/27; it takes as event horizon ρ_h the real root of ρ³ − 2ρ² + γ = 0 nearest to a tabulated guess, and prints it |
| `extreme/matrixassembler.mpl` | Maple procedure `MatrixAssembler(d, N, s, L, p)` for the extremal black hole, whose coefficients carry γ = 32/27 and ρ_h = 4/3 |
| `nonextreme/chaseQNMs_template.m`, `extreme/chaseQNMs_template.m` | MATLAB script, with the Advanpix toolbox: solves the problem at the three resolutions, keeps the common eigenvalues, writes the report and plots the spectra |
| `nonextreme/run_qnm_computation.py`, `extreme/run_qnm_computation.py` | Python driver: runs Maple, then MATLAB, for every value of γ (non-extremal) or of ℓ (extremal) |
| `utils.py`, `requirements.txt`, `.env.example` | helper functions and Python packages of both drivers, and an example of their configuration file `.env` |

The extremal MATLAB script and driver differ from the non-extremal ones only in the parameter γ, which they do not
take, in the folders they run through, and in the names and headers of the reports. For each set of parameters, the
driver creates the folders `data/`, `results/` and `shots/`, copies the template as `chaseQNMs.m`, and writes for each
resolution N a short Maple script that reads `matrixassembler.mpl` and calls `MatrixAssembler(N, N, s, L, g, path)`
(without `g` in the extremal case); command-line Maple runs the three resolutions in parallel, each writing
`data/M0_N.mat`, `data/M1_N.mat` and `data/M2_N.mat`. The driver then writes the parameters into `matlab_params.mat`
and starts MATLAB on `chaseQNMs.m`, which for each resolution sets the precision to N digits (`mp.Digits(N)`), reads
the matrices with `mp.read` and computes the 2N eigenvalues of `polyeig(M0, mp('1i')*M1, M2)`. Since the first argument
of the procedures is the number of digits, each resolution is thus assembled and solved with N digits. The script keeps
the eigenvalues of the first resolution that lie within 10⁻⁴ of an eigenvalue of each of the other two, and lists those
with a positive real part and the purely imaginary ones (real part below 10⁻¹⁶ in absolute value), by decreasing
imaginary part, in `results/report_QNM_s{s}_L{ℓ}_g{γ}_{date}_{time}.txt`, with γ written to five significant digits
(`g0.14815` for 4/27), or in `results/report_QNM_s{s}_L{ℓ}_extreme_{date}_{time}.txt` for the extremal black hole:
these are the names and the format of the stored reports. It also saves the three spectra in `results/qnms.mat` and
plots them, with the common eigenvalues circled, in `shots/chaseQNMs.png` (with copies in PDF, EPS and FIG formats).

The drivers read the parameters from the names of the folders. The non-extremal driver runs in a folder `s{s}/L{ℓ}/`
whose subfolders `g0`, `g1` and `gX_Y` (for γ = X/Y) are the values of γ to compute, and the extremal driver in a folder
`s{s}/` whose subfolders `L{ℓ}` are the multipoles. The reports of Paper 3 come from the following runs.

| Case | Scripts | Spin s, multipole ℓ and parameter γ | Reports |
|---|---|---|---|
| non-extremal | `nonextreme/` | s = 0, ℓ = 0 and 1: γ = 0, 4/27, 8/27, 12/27, 20/27, 24/27, 1 and 28/27; s = 1: γ = 0, 1, 11/10 and 59/50 for ℓ = 1, γ = 1, 11/10 and 59/50 for ℓ = 2, γ = 1 for ℓ = 3; s = 2: γ = 0, 1/2, 1 and 59/50 for ℓ = 2, γ = 1/2, 1 and 59/50 for ℓ = 3, γ = 1/2 for ℓ = 4 | 33 |
| near-extremal | `nonextreme/` | s = 0, ℓ = 0 and 1, γ = 59/50 | 2 |
| extremal | `extreme/` | γ = 32/27; ℓ = 0 to 5 for s = 0, 1 to 6 for s = 1, 2 to 7 for s = 2 | 18 |

All runs use the resolutions 300, 290 and 280, set by `RESOLUTIONS` in `.env`, except the second report of s = 1,
ℓ = 1, γ = 1, computed with `RESOLUTIONS=400,390,380` and renamed by hand with the suffix `_hires`. The driver passes γ
to Maple as the shortest decimal form of a double-precision number (`g := 0.14814814814814814` for 4/27), so that the
matrices of the reports are those of this decimal, which differs from the fraction by less than 10⁻¹⁶ (the values
11/10, 59/50, 1/2, 0 and 1 are exact).

To run a case again, gather the files of the case and the shared ones in a working folder named after the parameters,
create the subfolders, and write the configuration file `.env` from `.env.example` with the paths of Maple and of the
Advanpix toolbox. For the non-extremal scalar monopole at γ = 4/27 and γ = 1, from the folder of the repository:

```bash
mkdir -p work/s0/L0/g4_27 work/s0/L0/g1
cp matlab/Paper3HaywardBH/nonextreme/* matlab/Paper3HaywardBH/utils.py work/s0/L0/
cp matlab/Paper3HaywardBH/.env.example work/s0/L0/.env     # then set MAPLE_EXECUTABLE and ADVANPIX_PATH
python3 -m pip install -r matlab/Paper3HaywardBH/requirements.txt
cd work/s0/L0 && python3 run_qnm_computation.py
```

For the extremal black hole, the files of `extreme/` and the shared ones go into a folder such as `work/s1/` with the
subfolders `L1`, `L2` and so on, where the driver is run in the same way. The driver first empties the folders `data/`, `results/` and
`shots/` of every subfolder, so that the results of an earlier run are lost. It reports its progress to a Telegram bot
when `BOT_TOKEN` and `CHAT_ID` are set in `.env`, and sends nothing when they are left empty, as in `.env.example`. The
option `MATLAB_PARALLEL` of `.env` has no effect, since the driver does not pass it on to MATLAB, which always solves the
three resolutions one after the other. At the resolutions of the reports, one set of parameters took about a quarter of
an hour.

## The computed spectra

The folder `qnmdata/` holds one plain-text report per set of parameters:

| Folder | Content | Reports |
|---|---|---|
| `Paper1NegAlpha/s{0,1,2}/g{0_5,4_5}/` | Planck stars, γ = 1/2 and 9/2, M = 1, ten multipoles per spin (ℓ = 0 to 9 for s = 0, 1 to 10 for s = 1, 2 to 11 for s = 2) | 60 |
| `Paper2PosAlpha/extreme/s{0,1,2}/L*/` | extremal Bonanno–Reuter black hole, M = 3.50274 | 18 |
| `Paper2PosAlpha/nonextreme/s{0,1,2}/L*/M*/` | non-extremal Bonanno–Reuter black holes, several masses up to M = 1000, close to Schwarzschild | 55 |
| `Paper3HaywardBH/extreme/s{0,1,2}/L*/` | extremal Hayward black hole, γ = 32/27 | 18 |
| `Paper3HaywardBH/nearextreme/s0/L{0,1}/g59_50/` | Hayward black hole close to extremality, γ = 59/50 | 2 |
| `Paper3HaywardBH/nonextreme/s{0,1,2}/L*/g*/` | non-extremal Hayward black holes, several values of γ | 33 |

In the names, `s` and `L` give the spin and the multipole. The value of γ is written `g`, with the underscore as a
decimal point in Paper 1 (`g0_5` is γ = 0.5, `g4_5` is γ = 4.5) and as a fraction bar in Paper 3 (`g4_27` is
γ = 4/27, `g59_50` is γ = 59/50, `g0` and `g1` are 0 and 1). The mass is written `M`, with the underscore as a decimal
point (`M5_0` is M = 5, `M3_503` is M ≈ 3.503). The folder `Paper2PosAlpha/nonextreme/s0/L4/M5_0/` also holds the report
of M = 8, and `Paper3HaywardBH/nonextreme/s1/L1/g1/` a second report computed at the higher resolutions 400, 390 and 380,
which gives the same least damped mode.

Each report starts with the date, the parameters and the three resolutions (200, 190 and 180 in Paper 1, 400, 390 and
380 in Paper 2, 300, 290 and 280 in Paper 3), then lists the modes common to the three within the tolerance 10⁻⁴, by
decreasing imaginary part: the modes with a positive real part, then the purely imaginary ones with the gaps between
consecutive ones. The frequencies are ω in Planck units for the mass of the header in Papers 1 and 2, and the
dimensionless Mω in Paper 3. The Hayward article keeps the purely imaginary entries as unclassified spectral
candidates, not established quasinormal modes, since a reproducible sequence of eigenvalues on the imaginary axis need
not be a pole of the Green function.

## Contents

| Path | Content |
|---|---|
| `maple/` | five Maple 2024 worksheets with the derivations: Paper 1, Paper 2 (non-extremal and extremal) and Paper 3 (non-extremal and extremal) |
| `matlab/` | the five Maple procedures `MatrixAssembler` of Paper 1, Paper 2 (non-extremal and extremal) and Paper 3 (non-extremal and extremal); for Paper 3 also the MATLAB script that solves the eigenvalue problem and writes the reports, and the Python driver with its configuration example |
| `qnmdata/` | the 186 reports of computed spectra described above |
| `docs/figures/` | the script `plot_readme_figure.py` that draws the figure of this page from the reports of Paper 3, with the style file beside it |

## Requirements

Maple 2024 or later, in which the worksheets are saved; the procedures use only `LinearAlgebra`, `ChebyshevT`, `fsolve`
(Papers 1 and 2), `RootOf` with `allvalues` (Paper 3) and `ExportMatrix`, and those of Paper 3 were run with Maple 2022.
To solve the eigenvalue problems, MATLAB with the Advanpix Multiprecision Computing Toolbox (a commercial toolbox); the
MATLAB script of Paper 3 defines functions in the middle of the script, which MATLAB accepts from release R2024a on. The
driver of Paper 3 needs Python 3 with SciPy, python-dotenv and requests (`matlab/Paper3HaywardBH/requirements.txt`), and
calls Maple and MATLAB from the command line. The figure of this page needs Python 3 with NumPy and Matplotlib, and is
redrawn from `docs/figures/` by `python3 plot_readme_figure.py`.

## How to cite

If you use these codes or data, please cite the article of the paper they belong to.

Paper 1. D. Batic, D. Dutykh and F. Scardigli, Spectral analysis of quasinormal modes of Planck stars, The European
Physical Journal C 86, 165 (2026).

https://doi.org/10.1140/epjc/s10052-026-15430-8

```bibtex
@article{Batic2026,
  author = {Batic, Davide and Dutykh, Denys and Scardigli, Fabio},
  title = {Spectral Analysis of Quasinormal Modes of {{Planck}} Stars},
  journal = {The European Physical Journal C},
  year = 2026,
  month = feb,
  volume = {86},
  number = {2},
  pages = {165},
  doi = {10.1140/epjc/s10052-026-15430-8},
  eprint = {2602.19833},
  archiveprefix = {arXiv}
}
```

Paper 2. D. Batic, D. Dutykh and F. Scardigli, Quasinormal modes of Bonanno–Reuter black holes via the spectral method,
Physical Review D 114, 024023 (2026).

https://doi.org/10.1103/sgt6-kwfl

```bibtex
@article{Batic2026f,
  author = {Batic, Davide and Dutykh, Denys and Scardigli, Fabio},
  title = {Quasinormal Modes of {{Bonanno-Reuter}} Black Holes via the Spectral Method},
  journal = {Physical Review D},
  year = 2026,
  month = jul,
  volume = {114},
  number = {2},
  pages = {024023},
  doi = {10.1103/sgt6-kwfl},
  eprint = {2607.10199},
  archiveprefix = {arXiv}
}
```

Paper 3. D. Batic and D. Dutykh, Quasinormal modes of Hayward black holes from Schwarzschild to extremality: a
comparative spectral analysis, Physics of the Dark Universe 54, 102467 (2026).

https://doi.org/10.1016/j.dark.2026.102467

```bibtex
@article{Batic2026h,
  author = {Batic, Davide and Dutykh, Denys},
  title = {Quasinormal Modes of {{Hayward}} Black Holes from {{Schwarzschild}} to Extremality: {{A}} Comparative Spectral Analysis},
  journal = {Physics of the Dark Universe},
  year = 2026,
  month = dec,
  volume = {54},
  pages = {102467},
  doi = {10.1016/j.dark.2026.102467}
}
```

The button "Cite this repository" in the right-hand column of this page gives the first article, read from
`CITATION.cff`, which lists the other two among its references. Free copies of the three articles are on arXiv:

https://arxiv.org/abs/2602.19833

https://arxiv.org/abs/2607.10199

https://arxiv.org/abs/2609.33267

## Licence

GNU General Public License, version 3 (see `LICENSE`).

## Authors

Davide Batic and Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi) and
Fabio Scardigli (Dipartimento di Matematica, Politecnico di Milano). The articles are listed, with the other studies of
quasinormal modes of the same group, in the [Quasi-Normal Modes Hall of Fame](https://qnms.denys-dutykh.com/).

https://www.denys-dutykh.com/
