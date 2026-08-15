# SARS-CoV-2 Transmission Bottleneck Estimation

**Convergent Estimation of the SARS-CoV-2 Transmission Bottleneck via Joint Beta-Binomial and Presence-Absence Likelihood Modeling**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21941861.svg)](https://doi.org/10.5281/zenodo.21941861)
[![GitHub release](https://img.shields.io/github/v/release/kanishkanandhkumar/sars-cov-2-bottleneck)](https://github.com/kanishkanandhkumar/sars-cov-2-bottleneck/releases)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

---

## Table of Contents

- [Overview](#overview)
- [Key Findings](#key-findings)
- [Repository Structure](#repository-structure)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Data Availability](#data-availability)
- [Citation](#citation)
- [Authors](#authors)
- [License](#license)
- [Acknowledgments](#acknowledgments)

---

## Overview

This repository contains the complete analysis pipeline and data for the estimation of the SARS-CoV-2 transmission bottleneck size. The bottleneck size (*N*<sub>b</sub>) represents the number of founder virions that establish a new infection and is a critical parameter in viral evolution, immune escape prediction, and vaccine design.

### Statistical Methods Implemented

| Method | Description | Status |
|:-------|:------------|:------:|
| **Joint Beta-Binomial MLE** | Accounts for sequencing depth and overdispersion | Executed |
| **Presence-Absence Binomial** | Binary outcome model for validation | Executed |
| **Non-parametric Bootstrap** | 10,000 iterations for confidence intervals | Executed |
| **MAF Sensitivity Sweep** | Threshold analysis (1–5%) | Executed |
| **Theoretical Curve Fitting** | *P* = 1 − (1 − *p*)<sup>*N*<sub>b</sub></sup> | Executed |
| **Figure Generation** | All 4 manuscript figures | Executed |

### Key Findings

* **Transmission Bottleneck Size:** *N*<sub>b</sub> ≈ **4–5 virions**
* **Joint Beta-Binomial MLE:** 5.00 (95% CI: 1–13)
* **Presence-Absence Model:** 4.00 (95% CI: 1–40)
* **Bootstrap (10,000 reps):** 5.0 (95% CI: 1.8–12.4)
* **MAF Threshold Stability:** 4.0–5.5 across 1–5% thresholds
* **Theoretical Fit:** Near-perfect alignment with *N*<sub>b</sub> = 5

---

## Repository Structure

```text
sars-cov-2-bottleneck/
│
├── README.md                  # Main repository documentation
├── LICENSE                    # MIT License
├── requirements.txt           # Python dependencies
├── environment.yml            # Conda environment specification
├── .gitignore                 # Git ignore rules
│
├── data/
│   ├── raw/                   # Raw sequencing data links & metadata
│   ├── processed/             # Processed VCF files
│   │   ├── all_isnvs.vcf              # All intra-host SNVs
│   │   ├── all_isnvs_fixed.vcf        # Quality-filtered SNVs
│   │   └── all_isnvs_annotated.vcf    # Annotated SNVs
│   └── sample/                # Sample data for pipeline testing
│       └── sample_bottleneck_data.csv
│
├── scripts/
│   ├── analysis/              # Core bottleneck analysis scripts
│   │   ├── bottleneck_mle.py          # Joint Beta-Binomial MLE
│   │   ├── run_bootstrap_and_plot.py  # Bootstrap resampling
│   │   ├── run_sensitivity_and_validation.py # MAF sweep
│   │   ├── calculate_dnds.py          # dN/dS calculation
│   │   ├── calculate_shannon_entropy.py # Diversity metrics
│   │   ├── create_sample_data.py      # Sample data generator
│   │   └── run_all.py                 # Master analysis script
│   └── visualization/         # Figure generation scripts
│       ├── generate_figures.py        # Publication figure generator
│       └── plot_helpers.py            # Plotting utilities
│
├── results/
│   ├── figures/               # Output figures (Figures 1–4)
│   └── summary/               # Output summary JSON/CSV files
│
├── supplementary/             # Supplementary data tables
│   ├── dnds_transmission_results.csv
│   ├── shannon_entropy_results.csv
│   └── transmission_bottleneck_results.csv
│
├── docs/                      # Methodological & statistical documentation
│   ├── data_format.md
│   ├── methods_details.pdf
│   └── statistical_derivations.pdf
│
└── tests/                    # Unit testing suite
    ├── test_bottleneck.py
    └── test_figures.py
Installation
Option 1: Using Conda (Recommended)
Bash
git clone [https://github.com/kanishkanandhkumar/sars-cov-2-bottleneck.git](https://github.com/kanishkanandhkumar/sars-cov-2-bottleneck.git)
cd sars-cov-2-bottleneck
conda env create -f environment.yml
conda activate sars-cov2-bottleneck
Option 2: Using Pip
Bash
git clone [https://github.com/kanishkanandhkumar/sars-cov-2-bottleneck.git](https://github.com/kanishkanandhkumar/sars-cov-2-bottleneck.git)
cd sars-cov-2-bottleneck
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements.txt
Quick Start
Bash
# 1. Generate sample dataset
python scripts/analysis/create_sample_data.py

# 2. Run complete bottleneck estimation pipeline
python scripts/analysis/run_all.py --data data/sample/sample_bottleneck_data.csv

# 3. Generate publication-ready figures
python scripts/visualization/generate_figures.py --input data/sample/sample_bottleneck_data.csv --output results/figures/
Usage
Full Pipeline Execution
Bash
# Execute core analysis on processed variant data
python scripts/analysis/run_all.py --data data/processed/all_isnvs_fixed.vcf

# Run model sensitivity sweeps across MAF cutoffs
python scripts/analysis/run_sensitivity_and_validation.py --input data/processed/all_isnvs_fixed.vcf --output results/summary/
Data Availability
All raw sequencing datasets have been deposited in the NCBI Sequence Read Archive (SRA) / BioProject database under accession number PRJEB104853 (BioProject ID: 1421233). Processed variant call files (VCF), analytical scripts, and joint Beta-Binomial likelihood optimization code are publicly accessible without restriction in this repository and archived on Zenodo.

Citation
If you use this software or analytical pipeline, please cite the repository and primary manuscript:

Software & Repository
Code snippet
@software{kanishk2026bottleneck,
  author       = {Kanishk, A.},
  title        = {kanishkanandhkumar/sars-cov-2-bottleneck: v1.0.0 - SARS-CoV-2 Transmission Bottleneck Estimation},
  month        = feb,
  year         = 2026,
  publisher    = {Zenodo},
  version      = {v1.0.0},
  doi          = {10.5281/zenodo.21941861},
  url          = {[https://doi.org/10.5281/zenodo.21941861](https://doi.org/10.5281/zenodo.21941861)}
}
Primary Study
Code snippet
@article{khamduang2026household,
  author  = {Khamduang, W. and et al.},
  title   = {Household SARS-CoV-2 transmission during Omicron wave in Chiang Mai, Thailand: a prospective observational study},
  journal = {The Lancet Regional Health - Southeast Asia},
  volume  = {44},
  pages   = {100711},
  year    = {2026},
  doi     = {10.1016/j.lansea.2025.100711}
}
Authors & Contributors
Kanishk A. – Methodology, Software Development, Formal Statistical Analysis, Pipeline Engineering.

Khamduang W. et al. – Household Cohort Study Design, Clinical Sample Collection, Genomic Data Generation.

Infrastructure & Cores
High-Performance Computing Core Facility, University Research Institute (Computational support and bioinformatics resources).

Sequencing Core Facility, National Institute of Virology (Raw genome sequencing data generation).

License
This project is licensed under the MIT License - see the LICENSE file for details.
