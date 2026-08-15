#!/usr/bin/env python3
"""
Bottleneck Estimation Script for SARS-CoV-2 Transmission

This script implements the Joint Beta-Binomial MLE and Presence-Absence
binomial models for estimating the transmission bottleneck size.

Author: Kanishk K.
Date: 2024
"""

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import betabinom, binom
from scipy.special import beta, betaln
import argparse
import logging
from typing import Tuple, Dict, List

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class BottleneckEstimator:
    """Main class for bottleneck estimation."""
    
    def __init__(self):
        self.n_b_range = np.arange(1, 501)
    
    def beta_binomial_likelihood(self, N_b: int, donor_freq: float, 
                                 recipient_depth: int, recipient_count: int) -> float:
        """
        Calculate Beta-Binomial likelihood for a single variant.
        
        Args:
            N_b: Bottleneck size
            donor_freq: Donor minor allele frequency (p)
            recipient_depth: Sequencing depth at this locus (n)
            recipient_count: Minor allele count in recipient (k)
        
        Returns:
            Log-likelihood value
        """
        # Beta distribution parameters
        alpha = N_b * donor_freq
        beta_param = N_b * (1 - donor_freq)
        
        # Avoid numerical issues
        if alpha <= 0 or beta_param <= 0:
            return -np.inf
        
        # Beta-Binomial log-likelihood
        try:
            # Using scipy's betabinom
            ll = betabinom.logpmf(recipient_count, recipient_depth, alpha, beta_param)
            if np.isnan(ll) or np.isinf(ll):
                return -np.inf
            return ll
        except Exception as e:
            logger.warning(f"Error in likelihood calculation: {e}")
            return -np.inf
    
    def joint_likelihood(self, N_b: int, data: pd.DataFrame) -> float:
        """
        Calculate joint log-likelihood across all variants.
        
        Args:
            N_b: Bottleneck size
            data: DataFrame with columns ['donor_freq', 'recipient_depth', 'recipient_count']
        
        Returns:
            Total log-likelihood
        """
        total_ll = 0.0
        for _, row in data.iterrows():
            ll = self.beta_binomial_likelihood(
                N_b, 
                row['donor_freq'], 
                row['recipient_depth'], 
                row['recipient_count']
            )
            if np.isinf(ll):
                return -np.inf
            total_ll += ll
        return total_ll
    
    def optimize_bottleneck(self, data: pd.DataFrame) -> Dict:
        """
        Find MLE for bottleneck size.
        
        Args:
            data: DataFrame with variant data
        
        Returns:
            Dictionary with estimate and confidence intervals
        """
        logger.info(f"Optimizing bottleneck with {len(data)} variants")
        
        # Evaluate likelihood over search space
        likelihoods = []
        for N_b in self.n_b_range:
            ll = self.joint_likelihood(N_b, data)
            likelihoods.append(ll)
        
        likelihoods = np.array(likelihoods)
        
        # Find MLE
        max_idx = np.argmax(likelihoods)
        N_b_hat = self.n_b_range[max_idx]
        max_ll = likelihoods[max_idx]
        
        logger.info(f"MLE estimate: N_b = {N_b_hat}, Log-likelihood = {max_ll:.4f}")
        
        # Profile likelihood confidence intervals
        # 95% CI: ΔLL <= 1.92 for 1 parameter
        threshold = max_ll - 1.92
        
        # Find CI bounds
        valid_idx = np.where(likelihoods >= threshold)[0]
        if len(valid_idx) > 0:
            ci_lower = self.n_b_range[valid_idx[0]]
            ci_upper = self.n_b_range[valid_idx[-1]]
        else:
            ci_lower = 1
            ci_upper = self.n_b_range[-1]
        
        return {
            'n_b_hat': N_b_hat,
            'max_ll': max_ll,
            'ci_lower': ci_lower,
            'ci_upper': ci_upper,
            'likelihoods': likelihoods
        }
    
    def presence_absence_likelihood(self, N_b: int, data: pd.DataFrame) -> float:
        """
        Calculate presence-absence binomial likelihood.
        
        Args:
            N_b: Bottleneck size
            data: DataFrame with columns ['donor_freq', 'transmitted']
        
        Returns:
            Log-likelihood
        """
        total_ll = 0.0
        for _, row in data.iterrows():
            p = row['donor_freq']
            y = row['transmitted']  # 1 if transmitted, 0 if absent
            
            # P(transmission) = 1 - (1-p)^N_b
            p_trans = 1 - (1 - p) ** N_b
            
            # Binomial log-likelihood
            if y == 1:
                ll = np.log(p_trans) if p_trans > 0 else -np.inf
            else:
                ll = np.log(1 - p_trans) if p_trans < 1 else -np.inf
            
            if np.isinf(ll):
                return -np.inf
            total_ll += ll
        
        return total_ll
    
    def optimize_presence_absence(self, data: pd.DataFrame) -> Dict:
        """Optimize presence-absence model."""
        logger.info(f"Optimizing presence-absence model with {len(data)} variants")
        
        likelihoods = []
        for N_b in self.n_b_range:
            ll = self.presence_absence_likelihood(N_b, data)
            likelihoods.append(ll)
        
        likelihoods = np.array(likelihoods)
        
        # Find MLE
        max_idx = np.argmax(likelihoods)
        N_b_hat = self.n_b_range[max_idx]
        max_ll = likelihoods[max_idx]
        
        logger.info(f"Presence-absence MLE: N_b = {N_b_hat}, Log-likelihood = {max_ll:.4f}")
        
        # Profile likelihood CI
        threshold = max_ll - 1.92
        valid_idx = np.where(likelihoods >= threshold)[0]
        if len(valid_idx) > 0:
            ci_lower = self.n_b_range[valid_idx[0]]
            ci_upper = self.n_b_range[valid_idx[-1]]
        else:
            ci_lower = 1
            ci_upper = self.n_b_range[-1]
        
        return {
            'n_b_hat': N_b_hat,
            'max_ll': max_ll,
            'ci_lower': ci_lower,
            'ci_upper': ci_upper,
            'likelihoods': likelihoods
        }


def main():
    parser = argparse.ArgumentParser(description='Estimate SARS-CoV-2 transmission bottleneck')
    parser.add_argument('--input', type=str, required=True, help='Input CSV file with variant data')
    parser.add_argument('--output', type=str, default='results/tables/', help='Output directory')
    parser.add_argument('--model', type=str, choices=['joint', 'pa', 'both'], default='both',
                       help='Model to run: joint, pa (presence-absence), or both')
    args = parser.parse_args()
    
    # Load data
    logger.info(f"Loading data from {args.input}")
    data = pd.read_csv(args.input)
    
    # Initialize estimator
    estimator = BottleneckEstimator()
    
    results = {}
    
    # Run joint model
    if args.model in ['joint', 'both']:
        logger.info("Running Joint Beta-Binomial MLE...")
        joint_results = estimator.optimize_bottleneck(data)
        results['joint'] = joint_results
        logger.info(f"Joint MLE: N_b = {joint_results['n_b_hat']} "
                   f"(95% CI: {joint_results['ci_lower']}-{joint_results['ci_upper']})")
    
    # Run presence-absence model
    if args.model in ['pa', 'both']:
        logger.info("Running Presence-Absence model...")
        pa_results = estimator.optimize_presence_absence(data)
        results['pa'] = pa_results
        logger.info(f"PA MLE: N_b = {pa_results['n_b_hat']} "
                   f"(95% CI: {pa_results['ci_lower']}-{pa_results['ci_upper']})")
    
    # Save results
    import os
    os.makedirs(args.output, exist_ok=True)
    
    for model, result in results.items():
        output_file = f"{args.output}/{model}_results.json"
        import json
        with open(output_file, 'w') as f:
            json.dump(result, f, indent=2)
        logger.info(f"Results saved to {output_file}")


if __name__ == "__main__":
    main()
