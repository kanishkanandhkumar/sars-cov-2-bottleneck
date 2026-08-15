# SARS-CoV-2 Transmission Bottleneck Estimation

**Convergent Estimation of the SARS-CoV-2 Transmission Bottleneck via Joint Beta-Binomial and Presence-Absence Likelihood Modeling**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

## 📋 Overview

This repository contains the complete analysis pipeline and data for the estimation of the SARS-CoV-2 transmission bottleneck size using joint Beta-Binomial Maximum Likelihood Estimation (MLE) and Presence-Absence binomial modeling.

**Key Findings:**
- Transmission bottleneck size: **N_b ≈ 4-5 virions**
- Convergent estimates from two independent statistical models
- Robust across MAF thresholds (1-5%)
- Validated through bootstrapping and sensitivity analysis

## 📁 Repository Structure
publication_repo/
├── README.md # This file
├── LICENSE # MIT License
├── requirements.txt # Python dependencies
├── environment.yml # Conda environment specification
│
├── data/
│ ├── raw/ # Raw sequencing data (accessions provided)
│ └── processed/ # Processed VCF files and variant calls
│
├── scripts/
│ ├── analysis/ # Core analysis scripts
│ │ ├── bottleneck_mle.py # Joint Beta-Binomial MLE
│ │ ├── presence_absence.py # Presence-Absence binomial model
│ │ ├── bootstrap_analysis.py # Non-parametric bootstrap
│ │ └── sensitivity_sweep.py # MAF threshold sensitivity
│ └── visualization/ # Figure generation scripts
│ ├── figure1_cohort_flow.py
│ ├── figure2_bottleneck_est.py
│ ├── figure3_sensitivity.py
│ └── figure4_theoretical_fit.py
│
├── results/
│ ├── figures/ # Publication-ready figures (PNG, PDF)
│ │ ├── Fig1_flowchart.png
│ │ ├── Fig2_bottleneck.png
│ │ ├── Fig3_sensitivity.png
│ │ └── Fig4_theoretical.png
│ └── tables/ # Summary tables (CSV)
│ ├── table_s1_pair_metrics.csv
│ └── table_s2_sensitivity.csv
│
├── supplementary/ # Supplementary materials
│ ├── supplementary_tables.pdf
│ └── additional_analyses/
│
└── docs/ # Additional documentation
├── methods_details.pdf
└── statistical_derivations.pdf

## 🔧 Installation

### Option 1: Using Conda (Recommended)

# Clone the repository
git clone https://github.com/kanishkanandhkumar/Project-1.git
cd Project-1

# Create and activate conda environment
conda env create -f environment.yml
conda activate sars-cov2-bottleneck 

### Option 2: Using pip

# Clone the repository
git clone https://github.com/kanishkanandhkumar/Project-1.git
cd Project-1

# Install dependencies
pip install -r requirements.txt

### Usage
##Reproduce All Analyses

python scripts/analysis/run_all.py

##Run Individual Analyses

# Run bottleneck estimation
python scripts/analysis/bottleneck_mle.py --input data/processed/variants.vcf --output results/tables/

# Run bootstrap analysis
python scripts/analysis/bootstrap_analysis.py --iterations 10000

# Run sensitivity sweep
python scripts/analysis/sensitivity_sweep.py --thresholds 0.01 0.02 0.03 0.04 0.05

# Generate all figures
python scripts/visualization/generate_all_figures.py

### Data Availability
Raw sequencing reads: NCBI SRA BioProject PRJEB104853

Processed variant calls: Available in data/processed/

Analysis scripts: This repository

###Citation

Kanishk K., Doe J., Smith J., Taylor A.R. (2024). 
Convergent Estimation of the SARS-CoV-2 Transmission Bottleneck 
via Joint Beta-Binomial and Presence-Absence Likelihood Modeling.
Journal of Virology. [DOI: XXXX]

###License
This project is licensed under the MIT License - see the LICENSE file for details.

## Authors
Kanishk K. - Conceptualization, Methodology, Software, Analysis - @kanishkanandhkumar

Jane Doe - Data Curation, Investigation

John Smith - Resources, Supervision

Alex R. Taylor - Supervision, Funding

##Contact
For questions or collaborations, please contact: kanishkanandhkumar@gmail.com
 
