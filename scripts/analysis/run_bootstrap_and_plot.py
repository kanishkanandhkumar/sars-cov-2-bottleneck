import os, glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 1. Load summary results
if not os.path.exists("transmission_bottleneck_results.csv"):
    print("Error: transmission_bottleneck_results.csv not found.")
    exit(1)

res_df = pd.read_csv("transmission_bottleneck_results.csv")
informative = res_df.dropna(subset=['Nb_joint_mle']).copy()

print("\n" + "="*58)
print("             PER-PAIR BOTTLENECK SUMMARY TABLE           ")
print("="*58)
summary_table = informative[['donor_run', 'recipient_run', 'num_donor_isnvs', 'Nb_joint_mle', 'Nb_presence_absence']]
summary_table.columns = ['Donor', 'Recipient', 'Donor iSNVs', 'Nb (Joint MLE)', 'Nb (Pres-Abs)']
print(summary_table.to_string(index=False))
print("="*58 + "\n")

# 2. Perform 1,000 Bootstrap Resamples (Pair-level)
n_boot = 1000
np.random.seed(42)

mle_vals = informative['Nb_joint_mle'].values
pa_vals = informative['Nb_presence_absence'].values
n_pairs = len(mle_vals)

boot_medians_mle = []
boot_medians_pa = []

for _ in range(n_boot):
    idx = np.random.choice(n_pairs, size=n_pairs, replace=True)
    boot_medians_mle.append(np.median(mle_vals[idx]))
    boot_medians_pa.append(np.median(pa_vals[idx]))

ci_mle = np.percentile(boot_medians_mle, [2.5, 97.5])
ci_pa = np.percentile(boot_medians_pa, [2.5, 97.5])

print("="*58)
print("     95% CONFIDENCE INTERVALS (1,000 BOOTSTRAP RESAMPLES)    ")
print("="*58)
print(f"Joint Beta-Binomial MLE Median Nb: {np.median(mle_vals):.2f} (95% CI: {ci_mle[0]:.2f} - {ci_mle[1]:.2f})")
print(f"Presence-Absence Model Median Nb:   {np.median(pa_vals):.2f} (95% CI: {ci_pa[0]:.2f} - {ci_pa[1]:.2f})")
print("="*58 + "\n")

# 3. Generate Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

# Panel A: Per-pair values
pair_labels = [f"P{i+1}\n({row['Donor iSNVs']} iSNVs)" for i, row in summary_table.iterrows()]
x = np.arange(len(pair_labels))
width = 0.35

ax1.bar(x - width/2, summary_table['Nb (Joint MLE)'], width, label='Joint MLE', color='#2b5c8f')
ax1.bar(x + width/2, summary_table['Nb (Pres-Abs)'], width, label='Presence-Absence', color='#d95f02')
ax1.set_ylabel('Bottleneck Size ($N_b$)')
ax1.set_xlabel('Transmission Pairs')
ax1.set_title('Per-Pair Transmission Bottleneck Estimates')
ax1.set_xticks(x)
ax1.set_xticklabels(pair_labels)
ax1.legend()
ax1.grid(axis='y', linestyle='--', alpha=0.5)

# Panel B: Bootstrap distribution histogram
ax2.hist(boot_medians_mle, bins=12, alpha=0.7, color='#2b5c8f', label='Joint MLE', edgecolor='black')
ax2.hist(boot_medians_pa, bins=12, alpha=0.7, color='#d95f02', label='Presence-Absence', edgecolor='black')
ax2.axvline(np.median(mle_vals), color='#1b365d', linestyle='--', linewidth=2, label=f'MLE Median ({np.median(mle_vals):.1f})')
ax2.axvline(np.median(pa_vals), color='#a63603', linestyle='--', linewidth=2, label=f'PA Median ({np.median(pa_vals):.1f})')
ax2.set_xlabel('Median Bottleneck Size ($N_b$)')
ax2.set_ylabel('Bootstrap Count')
ax2.set_title('1,000 Bootstrap Resampling Distributions')
ax2.legend()
ax2.grid(axis='y', linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig('bottleneck_distribution.png', dpi=300)
print("Saved distribution plot to 'bottleneck_distribution.png'.")
