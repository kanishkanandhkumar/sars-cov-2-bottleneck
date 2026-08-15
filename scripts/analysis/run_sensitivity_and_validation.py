import os, glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binom

# 1. Load primary variant data
if not os.path.exists("dnds_transmission_results.csv"):
    print("Error: dnds_transmission_results.csv missing.")
    exit(1)

var_df = pd.read_csv("dnds_transmission_results.csv")
res_df = pd.read_csv("transmission_bottleneck_results.csv")
informative = res_df.dropna(subset=['Nb_joint_mle'])

# 2. Sensitivity Analysis across VAF thresholds (1% to 5%)
thresholds = [0.01, 0.02, 0.03, 0.04, 0.05]
sens_results = []

for tau in thresholds:
    # Recalculate presence-absence Nb across pairs for cutoff tau
    nb_pa_list = []
    for _, p in informative.iterrows():
        d_file = f"results/{p['donor_run']}_isnv_summary.csv"
        r_file = f"results/{p['recipient_run']}_isnv_summary.csv"
        if not os.path.exists(d_file) or not os.path.exists(r_file):
            continue
        df_d, df_r = pd.read_csv(d_file), pd.read_csv(r_file)
        donor_isnvs = df_d[(df_d['ALT_FREQ'] >= tau) & (df_d['ALT_FREQ'] < 0.50)] if not df_d.empty else pd.DataFrame()
        if donor_isnvs.empty:
            continue
            
        best_nb, best_ll = 1, -np.inf
        for Nb in range(1, 201):
            ll = 0.0
            for _, row in donor_isnvs.iterrows():
                pos, pd_val = row['POS'], row['ALT_FREQ']
                match = df_r[df_r['POS'] == pos] if not df_r.empty else pd.DataFrame()
                pr_val = float(match.iloc[0]['ALT_FREQ']) if not match.empty else 0.0
                p_trans = np.clip(1.0 - (1.0 - pd_val)**Nb, 1e-6, 1.0 - 1e-6)
                ll += np.log(p_trans if pr_val >= tau else (1.0 - p_trans))
            if ll > best_ll:
                best_ll, best_nb = ll, Nb
        nb_pa_list.append(best_nb)
        
    sens_results.append({
        'VAF_Threshold': f"{int(tau*100)}%",
        'Median_Nb': np.median(nb_pa_list) if nb_pa_list else np.nan,
        'Mean_Nb': np.mean(nb_pa_list) if nb_pa_list else np.nan,
        'Num_Pairs': len(nb_pa_list)
    })

sens_df = pd.DataFrame(sens_results)

print("\n" + "="*50)
print("     SENSITIVITY ANALYSIS ACROSS VAF CUTOFFS     ")
print("="*50)
print(sens_df.to_string(index=False))
print("="*50 + "\n")

# 3. Figure Generation (Supplementary Fig 1)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

# Panel A: Cutoff Sensitivity
ax1.plot(thresholds, sens_df['Median_Nb'], marker='o', linewidth=2, color='#1f77b4', label='Median $N_b$')
ax1.axhline(4.0, linestyle='--', color='gray', alpha=0.7, label='Primary Median ($N_b=4$)')
ax1.set_xlabel('iSNV Frequency Cutoff ($\tau$)')
ax1.set_ylabel('Estimated Bottleneck Size ($N_b$)')
ax1.set_title('A. Sensitivity to iSNV Detection Threshold')
ax1.set_xticks(thresholds)
ax1.set_xticklabels([f"{int(t*100)}%" for t in thresholds])
ax1.legend(frameon=False)
ax1.grid(True, linestyle=':', alpha=0.6)

# Panel B: Empirical vs Theoretical Transmission Probability Curve
p_grid = np.linspace(0.01, 0.50, 100)
ax2.plot(p_grid, 1.0 - (1.0 - p_grid)**5, color='#2ca02c', linewidth=2.5, label='Theoretical Model ($N_b = 5$)')
ax2.plot(p_grid, 1.0 - (1.0 - p_grid)**1, color='gray', linestyle=':', label='Theoretical Model ($N_b = 1$)')

# Binned empirical transmission rates
bins = np.linspace(0.03, 0.50, 6)
var_df['vaf_bin'] = pd.cut(var_df['donor_vaf'], bins)
binned = var_df.groupby('vaf_bin', observed=False)['transmitted'].agg(['mean', 'count', lambda x: np.mean(var_df.loc[x.index, 'donor_vaf'])]).reset_index()

ax2.scatter(binned['<lambda_0>'], binned['mean'], color='#d62728', s=70, zorder=5, label='Empirical Data Binned')
ax2.set_xlabel('Donor Allele Frequency ($p_D$)')
ax2.set_ylabel('Transmission Probability')
ax2.set_title('B. Empirical vs Theoretical Transmission Probability')
ax2.legend(frameon=False)
ax2.grid(True, linestyle=':', alpha=0.6)

plt.tight_layout()
plt.savefig('manuscript_supp_fig1.png', dpi=300)
plt.savefig('manuscript_supp_fig1.pdf')
print("Generated manuscript_supp_fig1.png and manuscript_supp_fig1.pdf.")
