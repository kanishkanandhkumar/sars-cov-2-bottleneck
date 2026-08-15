import os, glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

# 1. Load data files
res_df = pd.read_csv("transmission_bottleneck_results.csv")
entropy_df = pd.read_csv("shannon_entropy_results.csv") if os.path.exists("shannon_entropy_results.csv") else None
var_df = pd.read_csv("dnds_transmission_results.csv") if os.path.exists("dnds_transmission_results.csv") else None

informative = res_df.dropna(subset=['Nb_joint_mle']).copy()

# 2. Configure publication plot aesthetics
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
plt.rcParams['figure.dpi'] = 300

fig, axes = plt.subplots(2, 2, figsize=(11, 9))

# Panel A: Transmission Bottleneck per Pair
ax_a = axes[0, 0]
x = np.arange(len(informative))
width = 0.35
ax_a.bar(x - width/2, informative['Nb_joint_mle'], width, label='Joint MLE ($N_b=5$)', color='#1f77b4')
ax_a.bar(x + width/2, informative['Nb_presence_absence'], width, label='Pres-Abs ($N_b=4$)', color='#ff7f0e')
ax_a.set_ylabel('Bottleneck Size ($N_b$)')
ax_a.set_xlabel('Transmission Pairs')
ax_a.set_xticks(x)
ax_a.set_xticklabels([f"P{i+1}" for i in range(len(informative))])
ax_a.axhline(5.0, color='#1f77b4', linestyle='--', alpha=0.7)
ax_a.legend(frameon=False)
ax_a.set_title('A. Transmission Bottleneck Size ($N_b$)', loc='left', fontweight='bold')

# Panel B: Shannon Entropy Comparison
ax_b = axes[0, 1]
if entropy_df is not None:
    for _, row in entropy_df.iterrows():
        ax_b.plot(['Donor', 'Recipient'], [row['Donor_Entropy'], row['Recipient_Entropy']], 
                  color='#7f7f7f', alpha=0.5, marker='o', linewidth=1)
    ax_b.errorbar(['Donor', 'Recipient'], 
                  [entropy_df['Donor_Entropy'].mean(), entropy_df['Recipient_Entropy'].mean()],
                  yerr=[entropy_df['Donor_Entropy'].std(), entropy_df['Recipient_Entropy'].std()],
                  fmt='-o', color='#d62728', linewidth=2.5, capsize=4, label='Mean ± SD')
    ax_b.legend(frameon=False)
ax_b.set_ylabel('Shannon Entropy ($H$)')
ax_b.set_title('B. Intra-Host Diversity Shift', loc='left', fontweight='bold')

# Panel C: Donor VAF vs Recipient VAF Scatter
ax_c = axes[1, 0]
if var_df is not None:
    transmitted = var_df[var_df['transmitted'] == True]
    lost = var_df[var_df['transmitted'] == False]
    
    ax_c.scatter(lost['donor_vaf'], lost['recip_vaf'], color='#d62728', alpha=0.7, label='Lost (VAF < 3%)', s=45)
    ax_c.scatter(transmitted['donor_vaf'], transmitted['recip_vaf'], color='#2ca02c', alpha=0.8, label='Transmitted', s=45)
    
    # Calculate Correlation for transmitted
    if len(var_df) > 1:
        r_val, p_val = spearmanr(var_df['donor_vaf'], var_df['recip_vaf'])
        ax_c.text(0.05, 0.88, f"Spearman $r = {r_val:.2f}$\n$p = {p_val:.3f}$", transform=ax_c.transAxes, 
                  bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#cccccc'))
                  
    ax_c.axhline(0.03, color='gray', linestyle=':', label='Detection Threshold (3%)')
    ax_c.set_xlabel('Donor Allele Frequency (VAF)')
    ax_c.set_ylabel('Recipient Allele Frequency (VAF)')
    ax_c.legend(frameon=False, loc='lower right')
ax_c.set_title('C. iSNV Frequency Dynamics Across Bottleneck', loc='left', fontweight='bold')

# Panel D: iSNV Genomic Location Mapping
ax_d = axes[1, 1]
if var_df is not None:
    ax_d.scatter(var_df['pos'], var_df['donor_vaf'], c=np.where(var_df['transmitted'], '#2ca02c', '#d62728'), s=50, alpha=0.8)
    ax_d.set_xlabel('Genomic Position (bp)')
    ax_d.set_ylabel('Donor VAF')
    ax_d.grid(axis='x', linestyle=':', alpha=0.6)
ax_d.set_title('D. iSNV Distribution Across Genome', loc='left', fontweight='bold')

plt.tight_layout()
plt.savefig('manuscript_figure_1.png', dpi=300)
plt.savefig('manuscript_figure_1.pdf')
print("\nGenerated manuscript_figure_1.png (300 DPI) and manuscript_figure_1.pdf.")
