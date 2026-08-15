import os, glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

# 1. Load transmission pair results
res_df = pd.read_csv("transmission_bottleneck_results.csv")
pairs = res_df.dropna(subset=['Nb_joint_mle'])

def compute_sample_entropy(csv_path):
    if not os.path.exists(csv_path):
        return np.nan
    df = pd.read_csv(csv_path)
    if df.empty or 'ALT_FREQ' not in df.columns:
        return 0.0
    
    # Calculate H per iSNV position: - [p*log2(p) + (1-p)*log2(1-p)]
    entropies = []
    for _, row in df.iterrows():
        p = row['ALT_FREQ']
        if 0 < p < 1:
            h = -(p * np.log2(p) + (1 - p) * np.log2(1 - p))
            entropies.append(h)
    return np.mean(entropies) if entropies else 0.0

# 2. Compute entropy for donors and recipients in informative pairs
data = []
for _, row in pairs.iterrows():
    d_run, r_run = row['donor_run'], row['recipient_run']
    d_h = compute_sample_entropy(f"results/{d_run}_isnv_summary.csv")
    r_h = compute_sample_entropy(f"results/{r_run}_isnv_summary.csv")
    
    data.append({
        'Donor': d_run,
        'Recipient': r_run,
        'Donor_Entropy': d_h,
        'Recipient_Entropy': r_h,
        'Nb': row['Nb_joint_mle']
    })

entropy_df = pd.DataFrame(data)
entropy_df.to_csv("shannon_entropy_results.csv", index=False)

print("\n" + "="*50)
print("       SHANNON ENTROPY COMPARISON (PER PAIR)       ")
print("="*50)
print(entropy_df[['Donor', 'Recipient', 'Donor_Entropy', 'Recipient_Entropy']].to_string(index=False))
print("="*50)

# 3. Statistical Test (Paired Wilcoxon Signed-Rank Test)
stat, p_val = wilcoxon(entropy_df['Donor_Entropy'], entropy_df['Recipient_Entropy'])
print(f"\nMean Donor Entropy:    {entropy_df['Donor_Entropy'].mean():.4f}")
print(f"Mean Recipient Entropy: {entropy_df['Recipient_Entropy'].mean():.4f}")
print(f"Wilcoxon p-value:       {p_val:.4f}")

# 4. Visualization
plt.figure(figsize=(7, 5))
for idx, row in entropy_df.iterrows():
    plt.plot(['Donor', 'Recipient'], [row['Donor_Entropy'], row['Recipient_Entropy']], 
             marker='o', color='gray', alpha=0.6, linewidth=1.5)

plt.errorbar(['Donor', 'Recipient'], 
             [entropy_df['Donor_Entropy'].mean(), entropy_df['Recipient_Entropy'].mean()],
             yerr=[entropy_df['Donor_Entropy'].std(), entropy_df['Recipient_Entropy'].std()],
             fmt='-o', color='#d95f02', linewidth=3, capsize=5, label='Mean ± SD')

plt.ylabel('Mean Shannon Entropy ($H$)')
plt.title(f'Intra-Host Diversity Shift (p = {p_val:.4f})')
plt.grid(axis='y', linestyle='--', alpha=0.5)
plt.legend()
plt.tight_layout()
plt.savefig('shannon_entropy_shift.png', dpi=300)
print("\nSaved entropy shift plot to 'shannon_entropy_shift.png'.")
