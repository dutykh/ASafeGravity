# Asymptotically Safe Gravity QNM Repository

Collection of codes and supplemental materials for our papers on the computation of quasi-normal modes (QNMs) of black holes in asymptotically safe gravity.

## Repository Structure

This repository contains source codes and computed data for QNM calculations in different programming environments, as well as the resulting datasets. The main directories are:

- `matlab/` — Matlab scripts and functions for QNM computations, organized by paper:
  - `Paper1NegAlpha/`, `Paper2PosAlpha/`, `Paper3HaywardBH/`
    - Example: `matrixassembler.mpl` in `Paper1NegAlpha/`
- `maple/` — Maple worksheets and code for symbolic and numeric QNM calculations, organized by paper:
  - `Paper1NegAlpha/`, `Paper2PosAlpha/`, `Paper3HaywardBH/`
    - Example: `01-SUPPLEMENTAL MATERIAL-PLANCK.mw` in `Paper1NegAlpha/`
- `qnmdata/` — Computed QNM data files, organized by paper and further subdivided by physical parameters (e.g., spin, coupling):
  - `Paper1NegAlpha/`, `Paper2PosAlpha/`, `Paper3HaywardBH/`
    - Example: `s0/g0_5/report_QNM_s0_L0_g0.5_M1_0.txt` in `Paper1NegAlpha/`

Each paper directory contains the relevant code and/or data for the corresponding publication. See the subdirectory structure for details on specific parameter sets and results.

## Usage

Browse the respective directories for codes and data associated with each paper. For questions or to reproduce results, refer to the documentation and scripts in the relevant subfolders.
