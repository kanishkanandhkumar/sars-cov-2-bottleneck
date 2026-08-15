#!/usr/bin/env python3
"""
Create sample data for testing and demonstration.

This script generates a small sample dataset that mimics the structure
of the real data for users to test the pipeline.

Author: Kanishk K.
Date: 2024
"""

import pandas as pd
import numpy as np
import os
from pathlib import Path

def create_sample_data(output_dir="data/sample"):
    """Create sample dataset."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Set random seed for reproducibility
    np.random.seed(42)
    
    # Simulate 9 transmission pairs with 3-5 variants each
    data = []
    
    pairs = [
        ("Pair_01", 2.0),
        ("Pair_02", 5.0),
        ("Pair_03", 4.0),
        ("Pair_04", 8.0),
        ("Pair_05", 3.0),
        ("Pair_06", 6.0),
        ("Pair_07", 12.0),
        ("Pair_08", 4.0),
        ("Pair_09", 5.0),
    ]
    
    for pair_id, true_nb in pairs:
        n_variants = np.random.randint(2, 6)
        for i in range(n_variants):
            # Donor frequency (uniform between 0.01 and 0.50)
            donor_freq = np.random.uniform(0.01, 0.50)
            
            # Recipient depth (1000-5000)
            recipient_depth = np.random.randint(1000, 5001)
            
            # Probability of transmission
            p_trans = 1 - (1 - donor_freq) ** true_nb
            
            # Transmitted or not
            transmitted = 1 if np.random.random() < p_trans else 0
            
            # If transmitted, count is binomial from depth
            if transmitted:
                recipient_count = np.random.binomial(recipient_depth, donor_freq)
                # Ensure at least 1 read
                recipient_count = max(1, recipient_count)
            else:
                recipient_count = 0
            
            data.append({
                'pair_id': pair_id,
                'variant_id': f"{pair_id}_var_{i+1}",
                'donor_freq': round(donor_freq, 4),
                'recipient_depth': recipient_depth,
                'recipient_count': recipient_count,
                'transmitted': transmitted,
                'true_nb': true_nb
            })
    
    df = pd.DataFrame(data)
    
    # Save as CSV
    output_file = Path(output_dir) / "sample_bottleneck_data.csv"
    df.to_csv(output_file, index=False)
    print(f"Sample data saved to {output_file}")
    print(f"Total variants: {len(df)}")
    print("\nData summary:")
    print(df.groupby('pair_id').agg({
        'true_nb': 'first',
        'variant_id': 'count',
        'transmitted': 'sum'
    }).rename(columns={'variant_id': 'n_variants', 'transmitted': 'n_transmitted'}))
    
    return df

if __name__ == "__main__":
    create_sample_data()
