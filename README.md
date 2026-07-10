# Asymptotically Safe Gravity

![Asymptotically Safe Gravity](pics/asg.png)

Collection of codes and supplemental materials for our series of three papers on the computation of quasi-normal modes (QNMs) of black holes in asymptotically safe gravity.

## Authors

- Dr. Davide Batic (Mathematics Department, Khalifa University of Science and Technology, PO Box 127788, Abu Dhabi, United Arab Emirates)
- Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, PO Box 127788, Abu Dhabi, United Arab Emirates)
- Dr. Fabio Scardigli (Dipartimento di Matematica, Politecnico di Milano, Piazza Leonardo da Vinci 32, 20133 Milano, Italy)

## Publications

This repository accompanies a series of three articles. The directory names
`Paper1NegAlpha/`, `Paper2PosAlpha/`, and `Paper3HaywardBH/` correspond to the
papers listed below.

1. **Paper 1 (published).** D. Batic, D. Dutykh, and F. Scardigli,
   *Spectral analysis of quasinormal modes of Planck stars*,
   Eur. Phys. J. C **86**, 165 (2026).
   [doi:10.1140/epjc/s10052-026-15430-8](https://doi.org/10.1140/epjc/s10052-026-15430-8)
2. **Paper 2 (published).** D. Batic, D. Dutykh, and F. Scardigli,
   *Quasinormal modes of Bonanno–Reuter black holes via the spectral method*,
   Phys. Rev. D **114**, 024023 (2026).
   [doi:10.1103/sgt6-kwfl](https://doi.org/10.1103/sgt6-kwfl)
3. **Paper 3 (in preparation).** D. Batic and D. Dutykh,
   *Quasinormal modes of Hayward black holes via the spectral method*
   (manuscript in preparation).

### How to cite

If you use these codes or data, please cite the relevant paper. BibTeX for the
published papers:

```bibtex
@article{Batic2026PlanckStars,
  author  = {Batic, Davide and Dutykh, Denys and Scardigli, Fabio},
  title   = {Spectral analysis of quasinormal modes of {P}lanck stars},
  journal = {The European Physical Journal C},
  volume  = {86},
  pages   = {165},
  year    = {2026},
  doi     = {10.1140/epjc/s10052-026-15430-8},
}

@article{Batic2026BonannoReuter,
  author  = {Batic, Davide and Dutykh, Denys and Scardigli, Fabio},
  title   = {Quasinormal modes of {B}onanno--{R}euter black holes via the spectral method},
  journal = {Physical Review D},
  volume  = {114},
  pages   = {024023},
  year    = {2026},
  doi     = {10.1103/sgt6-kwfl},
}
```

## Repository Structure

This repository contains source codes and computed data for QNM calculations in different programming environments, as well as the resulting datasets. The main directories are:

- `matlab/`: Matlab scripts and functions for QNM computations, organized by paper:
  - `Paper1NegAlpha/`, `Paper2PosAlpha/`
    - Example: `matrixassembler.mpl` in `Paper1NegAlpha/`
- `maple/`: Maple worksheets and code for symbolic and numeric QNM calculations, organized by paper:
  - `Paper1NegAlpha/`, `Paper2PosAlpha/`, `Paper3HaywardBH/`
    - Example: `01-SUPPLEMENTAL MATERIAL-PLANCK.mw` in `Paper1NegAlpha/`
- `qnmdata/`: Computed QNM data files, organized by paper and further subdivided by physical parameters (e.g., spin, coupling):
  - `Paper1NegAlpha/`, `Paper2PosAlpha/`, `Paper3HaywardBH/`
    - Example: `s0/g0_5/report_QNM_s0_L0_g0.5_M1_0.txt` in `Paper1NegAlpha/`

Each paper directory contains the relevant code and/or data for the corresponding publication. See the subdirectory structure for details on specific parameter sets and results.

## Usage

Browse the respective directories for codes and data associated with each paper. For questions or to reproduce results, refer to the documentation and scripts in the relevant subfolders.

<table>
<tr>
<td><img src="pics/code.jpg" alt="Code" width="400"/></td>
<td><img src="pics/bhcode.png" alt="Black Hole Code" width="400"/></td>
</tr>
</table>
