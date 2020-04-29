# beta_vs_pval.py
# plot beta vs p value

import argparse
from cyvcf2 import VCF
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import make_violin_plot # to find beta
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

def build_args():
    parser = argparse.ArgumentParser(description='Graph p values vs beta values')
    parser.add_argument('--rnafile', default='Data/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt')
    parser.add_argument('--vcffile', default='allSVs_1MB.vcf.gz')
    return parser.parse_args()

def get_samples(args):
    df = pd.read_csv(args.rnafile, sep='\t', index_col=0)
    pop_df = pd.read_csv('Data/integrated_call_samples_v3.20130502.ALL.panel', sep='\t', index_col=0)
    geuvadis_samples = df.iloc[:,3:].keys().values
    svs_samples = pop_df.index.values
    samples = np.intersect1d(geuvadis_samples, svs_samples)
    return samples

def calc_p(record, df, gene, samples, sample_ord):
    affected_samples = np.array([])
    for sample in samples:
        gt = record.genotypes[sample_ord.index(sample)][0] + record.genotypes[sample_ord.index(sample)][1]
        if gt > 0:
            affected_samples = np.append(affected_samples, sample)

    unaffected_samples = np.setdiff1d(samples, affected_samples)

    dat = df.loc[gene].iloc[3:].to_numpy()
    aff_dat = df.loc[gene][affected_samples].to_numpy()
    unaff_dat = df.loc[gene][unaffected_samples].to_numpy()
    aff_dat = aff_dat-np.mean(dat)
    unaff_dat = unaff_dat-np.mean(dat)
    aff_dat = aff_dat/np.std(dat)
    unaff_dat = unaff_dat/np.std(dat)

    _, p = mannwhitneyu(aff_dat, unaff_dat, alternative='two-sided')
    return p

args = build_args()
keep_samples = get_samples(args)
qtls = VCF(args.vcffile, samples=keep_samples.tolist())
sample_order = qtls.samples

df = pd.read_csv(args.rnafile, sep='\t', index_col=0)
df = df[keep_samples]

pvals = []
betavals = []
for record in qtls:
    gene = record.INFO.get('gene')
    p = calc_p(record, df, gene, keep_samples, sample_order)
    beta = make_violin_plot.find_beta(record, df, keep_samples, sample_order, verbose=False)
    pvals.append(p)
    betavals.append(abs(beta))

plt.scatter(pvals, betavals)
model = LinearRegression()
model.fit(np.log10(np.array(pvals)).reshape(-1, 1), np.array(betavals).reshape(-1, 1))
beta_preds = model.predict(np.log10(pvals).reshape(-1,1))
plt.plot(pvals, beta_preds, color='orange')
print('Coefficient of determination: %.2f'
      % r2_score(betavals, beta_preds))
plt.xscale('log')
plt.xlabel('P values (log)')
plt.ylabel('Beta values (absolute value)')
plt.title('P vs Beta')
plt.show()
