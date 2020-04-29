# make_quantile_plot.py
# Makes quantiles of beta values and plots statistics on the quantiles
# uses VCF file with pvalues and exonic from determine_coding_vs_noncoding.py

import make_violin_plot # to find beta
import argparse
import numpy as np
import pandas as pd
from cyvcf2 import VCF, Writer
import matplotlib.pyplot as plt
import ast
import mygene

def build_args():
    parser = argparse.ArgumentParser(description='Make quantiles of beta values')
    parser.add_argument('--rnafile', default='Data/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt')
    parser.add_argument('--vcffile', default='allSVs_withexonic.vcf.gz')
    parser.add_argument('--exonfile', default='exons.csv')
    parser.add_argument('--savequantiles', const=True, default=False, nargs='?')
    return parser.parse_args()

def get_gene_type(record, exon_df):
    gene = record.INFO.get('gene')
    gene = gene.split('.')[0]
    if pd.isna(exon_df.loc[gene,'ensembl.type_of_gene']):
        dict_list = ast.literal_eval(gene['ensembl'])
        type = dict_list[0]['type_of_gene']
    else:
        type = exon_df.loc[gene,'ensembl.type_of_gene']
    pseudos = set(['polymorphic_pseudogene', 'processed_pseudogene', 'unprocessed_pseudogene', 'transcribed_unitary_pseudogene', 'transcribed_processed_pseudogene', 'transcribed_unprocessed_pseudogene', 'unitary_pseudogene'])
    igs = set(['IG_C_gene', 'IG_V_gene', 'IG_V_pseudogene', 'IG_C_pseudogene'])
    if type in pseudos:
        type = 'pseudogene'
    if type in igs:
        type = 'IG Gene'
    return type

args = build_args()

keep_samples = make_violin_plot.get_samples(args)
qtls = VCF(args.vcffile, samples=keep_samples.tolist())
sample_order = qtls.samples

df = pd.read_csv(args.rnafile, sep='\t', index_col=0)
df = df[keep_samples]

beta_df = pd.DataFrame(columns=['Beta', 'Exonic', 'Type', 'ID', 'GeneType'])

exon_df = pd.read_csv(args.exonfile, index_col=0)

print('Finding beta values...')
for record in qtls:
    exonic = record.INFO.get('exonic')
    if exonic == 'y':
        coding = True
    else:
        coding = False

    beta = make_violin_plot.find_beta(record, df, keep_samples, sample_order, verbose=False)
    print('{}\t{}\t{}'.format(record.INFO['gene'], record.ID, beta))
    if 'SVLEN' in record.INFO:
        size = record.INFO.get('SVLEN')
    else:
        size = record.end - record.start
    type = get_gene_type(record, exon_df)
    beta_df = beta_df.append({'Beta': abs(beta), 'Beta_Sign': 0 if beta == 0 else beta/abs(beta), 'Exonic': coding, 'Type': record.INFO.get('SVTYPE'), 'ID': record.ID, 'Size': size, 'GeneType': type}, ignore_index=True)

q = [0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
quantiles = beta_df.quantile(q=q)
quantile_names = ['(0,0.5]','(0.5,0.6]','(0.6,0.7]','(0.7,0.8]','(0.8,0.9]','(0.9,1.0]']
sv_types = {}
sv_sizes = {}
gene_types = {}
for i in range(len(q)-1):
    q_up = quantiles.loc[q[i+1], 'Beta']
    q_low = quantiles.loc[q[i], 'Beta']
    quant = beta_df[np.logical_and(beta_df['Beta'] > q_low, beta_df['Beta'] <= q_up)]
    sv_types[quantile_names[i]] = quant['Type'].value_counts()
    sv_sizes[quantile_names[i]] = pd.cut(quant['Size'], bins=np.logspace(0,7,8), include_lowest=True, labels=['(10^0,10^1]','(10^1,10^2]','(100^2,10^3]','(100^3,10^4]','(100^4,10^5]','(100^5,10^6]','(100^6,10^7]']).value_counts()
    gene_types[quantile_names[i]] = quant['GeneType'].value_counts()

ax = pd.DataFrame(sv_types).transpose().plot(kind='bar', stacked=True)
ax.set_xlabel('Quantile')
ax.set_ylabel('SV/Gene Pair Count')
ax.set_title('Beta Value Quantiles for SV Types')
plt.show()

ax = pd.DataFrame(sv_sizes).transpose().plot(kind='bar', stacked=True)
ax.set_xlabel('Quantile')
ax.set_ylabel('SV/Gene Pair Count')
ax.set_title('Beta Value Quantiles for SV Sizes')
plt.show()

ax = pd.DataFrame(gene_types).transpose().plot(kind='bar', stacked=True)
ax.set_xlabel('Quantile')
ax.set_ylabel('SV/Gene Pair Count')
ax.set_title('Beta Value Quantiles for Gene Type')
plt.show()

# plt.plot(quantile_names, enrichment)
# plt.xlabel('Quantile')
# plt.ylabel('Enrichment in Exon Regions')
# plt.title('Exonic Enrichment')
# plt.show()
qtls.close()

if (args.savequantiles):
# save each quantile in a separate VCF
    qtls = VCF(args.vcffile, samples=keep_samples.tolist())
    sample_order = qtls.samples
    qtls.add_info_to_header({'ID': 'beta', 'Description': 'beta value',
        'Type':'Character', 'Number': '1'})
    writer_list = [Writer('quantile_vcfs/'+quantile_names[i]+'.vcf.gz', qtls) for i in range(len(quantile_names))]
    for record in qtls:
        id = record.ID
        beta = beta_df[beta_df['ID'] == id].iloc[0]['Beta']
        beta_sign = beta_df[beta_df['ID'] == id].iloc[0]['Beta_Sign']
        record.INFO['beta'] = beta * beta_sign
        qtil = [i for i in range(len(quantiles)-1) if quantiles.loc[q[i],'Beta'] < beta]
        if len(qtil) == 0:
            continue
        writer_list[qtil[-1]].write_record(record)
    for w in writer_list:
        w.close()
    qtls.close()
