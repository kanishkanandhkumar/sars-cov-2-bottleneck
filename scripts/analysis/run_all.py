#!/usr/bin/env python3
"""
Master Script to Reproduce All Analyses

This script runs the complete analysis pipeline for SARS-CoV-2
transmission bottleneck estimation.

Usage:
    python scripts/analysis/run_all.py --data data/processed/all_isnvs.vcf

Author: Kanishk K.
Date: 2024
"""

import os
import sys
import subprocess
import argparse
import logging
from pathlib import Path

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Get project root
PROJECT_ROOT = Path(__file__).parent.parent.parent
SCRIPTS_DIR = PROJECT_ROOT / "scripts" / "analysis"
RESULTS_DIR = PROJECT_ROOT / "results"


def run_command(cmd, description):
    """Run a command and log output."""
    logger.info(f"Running: {description}")
    logger.info(f"Command: {' '.join(cmd)}")
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        logger.error(f"Error running {description}")
        logger.error(f"STDERR: {result.stderr}")
        sys.exit(1)
    
    logger.info(f"Successfully completed: {description}")
    if result.stdout:
        logger.info(f"Output: {result.stdout[:500]}...")
    
    return result


def main():
    parser = argparse.ArgumentParser(description='Run complete bottleneck analysis')
    parser.add_argument('--data', type=str, required=True,
                       help='Path to input VCF file')
    parser.add_argument('--iterations', type=int, default=10000,
                       help='Number of bootstrap iterations')
    parser.add_argument('--output', type=str, default=str(RESULTS_DIR),
                       help='Output directory')
    args = parser.parse_args()
    
    logger.info("=" * 60)
    logger.info("SARS-CoV-2 Transmission Bottleneck Analysis Pipeline")
    logger.info("=" * 60)
    
    # Step 1: Run bottleneck MLE
    run_command([
        sys.executable,
        str(SCRIPTS_DIR / "bottleneck_mle.py"),
        "--input", args.data,
        "--output", str(RESULTS_DIR / "summary"),
        "--model", "both"
    ], "Joint Beta-Binomial MLE and Presence-Absence models")
    
    # Step 2: Run bootstrap analysis
    run_command([
        sys.executable,
        str(SCRIPTS_DIR / "run_bootstrap_and_plot.py"),
        "--input", args.data,
        "--iterations", str(args.iterations),
        "--output", str(RESULTS_DIR / "summary")
    ], "Bootstrap resampling analysis")
    
    # Step 3: Run sensitivity sweep
    run_command([
        sys.executable,
        str(SCRIPTS_DIR / "run_sensitivity_and_validation.py"),
        "--input", args.data,
        "--output", str(RESULTS_DIR / "summary")
    ], "MAF threshold sensitivity sweep")
    
    # Step 4: Generate figures
    run_command([
        sys.executable,
        str(PROJECT_ROOT / "scripts" / "visualization" / "generate_figures.py"),
        "--input", args.data,
        "--output", str(RESULTS_DIR / "figures")
    ], "Generate manuscript figures")
    
    logger.info("=" * 60)
    logger.info("Analysis complete!")
    logger.info(f"Results saved to: {RESULTS_DIR}")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
