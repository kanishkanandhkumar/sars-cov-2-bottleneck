import os, glob, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact

# 1. Load transmission pair results
res_df = pd.read_csv("transmission_bottleneck_results.csv")
pairs = res_df.dropna(subset=['Nb_joint_mle'])

# 2. Helper function to classify variants as Non-synonymous (N) or Synonymous (S)
def classify_variant(row):
    for col in ['EFFECT', 'CONSEQUENCE', 'MUTATION_TYPE', 'ANN', 'TYPE', 'ANNOTATION', 'GENE_EFFECT']:
        if col in row and pd.notna(row[col]):
            val = str(row[col]).lower()
            if any(term in val for term in ['missense', 'nonsynonymous', 'non_synonymous', 'protein_altering', 'stop_gained']):
                return 'N'
            elif any(term in val for term in ['synonymous', 'silent']):
                return 'S'
    
    if 'AA_CHANGE' in row and pd.notna(row['AA_CHANGE']):
        match = re.search(r'([A-Z])\d+([A-Z])', str(row['AA_CHANGE']), re.IGNORECASE)
        if match:
            return 'S' if match.group(1).upper() == match.group(2).upper() else 'N'
            
    return 'N'

# 3. Process variants across all informative pairs
variants_data = []

for _, row in pairs.iterrows():
    d_run, r_run = row['donor_run'], row['recipient_run']
    d_file, r_file = f"results/{d_run}_isnv_summary.csv", f"results/{r_run}_isnv_summary.csv"
    
    if not os.path.exists(d_file) or not os.path.exists(r_file):
        continue
        
    df_d, df_r = pd.read_csv(d_file), pd.read_csv(r_file)
    donor_isnvs = df_d[(df_d['ALT_FREQ'] >= 0.03) & (df_d['ALT_FREQ'] < 0.50)] if not df_d.empty else pd.DataFrame()
    
    for _, var in donor_isnvs.iterrows():
        pos, pd_freq = var['POS'], var['ALT_FREQ']
        effect = classify_variant(var)
        
        match = df_r[df_r['POS'] == pos] if not df_r.empty else pd.DataFrame()
        pr_freq = float(match.iloc[0]['ALT_FREQ']) if not match.empty else 0.0
        
        variants_data.append({
            'donor_run': d_run,
            'recipient_run': r_run,
            'pos': pos,
            'donor_vaf': pd_freq,
            'recip_vaf': pr_freq,
            'effect': effect,
            'transmitted': pr_freq >= 0.03
        })

var_df = pd.DataFrame(variants_data)

if var_df.empty:
    print("Error: No donor iSNVs available to calculate dN/dS.")
    exit(1)

# 4. Calculate N/S Ratios and Selection Significance
n_trans_N = len(var_df[(var_df['transmitted'] == True) & (var_df['effect'] == 'N')])
n_trans_S = len(var_df[(var_df['transmitted'] == True) & (var_df['effect'] == 'S')])
n_lost_N  = len(var_df[(var_df['transmitted'] == False) & (var_df['effect'] == 'N')])
n_lost_S  = len(var_df[(var_df['transmitted'] == False) & (var_df['effect'] == 'S')])

tot_N, tot_S = n_trans_N + n_lost_N, n_trans_S + n_lost_S

# Fisher's Exact Test Contingency Table: [[Transmitted_N, Lost_N], [Transmitted_S, Lost_S]]
odds_ratio, p_value = fisher_exact([[n_trans_N, n_lost_N], [n_trans_S, n_lost_S]])

print("\n" + "="*58)
print("       dN/dS (N/S RATIO) TRANSMISSION SELECTION ANALYSIS   ")
print("="*58)
print(f"Total Donor iSNVs Analyzed:     {len(var_df)}")
print(f"  - Non-synonymous (N):        {tot_N}")
print(f"  - Synonymous (S):            {tot_S}")
print(f"  - Overall N/S Ratio:          {tot_N / max(1, tot_S):.2f}")
print("-" * 58)
print(f"Transmitted iSNVs (VAF >= 3%):  {n_trans_N + n_trans_S}")
print(f"  - Transmitted N / S:          {n_trans_N} / {n_trans_S} (Ratio: {n_trans_N / max(1, n_trans_S):.2f})")
print(f"Lost iSNVs (VAF < 3%):          {n_lost_N + n_lost_S}")
print(f"  - Lost N / S:                 {n_lost_N} / {n_lost_S} (Ratio: {n_lost_N / max(1, n_lost_S):.2f})")
print("="*58)
print("FISHER'S EXACT TEST FOR SELECTION:")
print(f"  - Odds Ratio (Transmitted vs Lost N/S): {odds_ratio:.2f}")
print(f"  - p-value:                             {p_value:.4f}")
print("="*58 + "\n")

var_df.to_csv("dnds_transmission_results.csv", index=False)

# 5. Visualization
plt.figure(figsize=(8, 5))
categories = ['All Donor iSNVs', 'Transmitted iSNVs', 'Lost iSNVs']
x = np.arange(len(categories))
width = 0.35

plt.bar(x - width/2, [tot_N, n_trans_N, n_lost_N], width, label='Non-synonymous (N)', color='#d95f02')
plt.bar(x + width/2, [tot_S, n_trans_S, n_lost_S], width, label='Synonymous (S)', color='#2b5c8f')

plt.ylabel('Number of iSNVs')
plt.title(f'iSNV Selection Analysis Across Bottleneck (Fisher p = {p_value:.4f})')
plt.xticks(x, categories)
plt.legend()
plt.grid(axis='y', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig('dnds_selection_plot.png', dpi=300)
print("Saved selection plot to 'dnds_selection_plot.png'.")
