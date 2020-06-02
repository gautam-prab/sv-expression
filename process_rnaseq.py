# process RNA seq data with Pandas

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from cyvcf2 import VCF, Writer
import pickle
from scipy.stats import mannwhitneyu

def build_args():
    parser = argparse.ArgumentParser(description='Find closeby structural variants to genes with differential expression')
    parser.add_argument('--p', default=0.05/128065, type=float, help='P-value cutoff') # 12816 for 100kb tolerance, 128065 for 1Mb
    parser.add_argument('--rnafile', default='Data/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt', help='TSV of RNA seq data for samples of interest.')
    parser.add_argument('--vcffile', default='Data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz', help='VCF file containing genotypes for variants of interest. Sample set must be a subset of samples in RNA file.')
    parser.add_argument('--outfile', default=None, help='Output file, a VCF which will contain only the variants deemed significant')
    parser.add_argument('--tolerance', default=1000000, type=int, help='cis region near genes to consider (default: 1,000,000 bases)')
    parser.add_argument('--maxchrom', default=None, help='set to 22 if you want to do autosomes only')
    parser.add_argument('--maf', default=0, type=float, help='Minor allele frequency cutoff; only consider those with higher MAF')
    parser.add_argument('--pfile', default='pvalues.csv', help='File to save p-values to for p-value cutoff analysis')
    return parser.parse_args()

# find closeby structural variants to an ORF
def find_closeby_svs(vcf, chrom, start, end, tol):
    return vcf('{}:{}-{}'.format(chrom, max(0, start-tol), end+tol))

def create_dataframe(args):
    df1 = pd.read_csv(args.rnafile, sep='\t', index_col=0)
    geuvadis_samples = df1.iloc[:,3:].keys().values
    return df1, geuvadis_samples

args = build_args()

df, geuvadis_samples = create_dataframe(args)
if args.maxchrom:
    df = df[pd.to_numeric(df['Chr'], errors='coerce') <= int(args.maxchrom)]
fig = plt.figure()

vcf = VCF(args.vcffile, gts012=True)
samples = np.asarray(vcf.samples)
if args.outfile:
    vcf.add_info_to_header({'ID': 'gene', 'Description': 'overlapping gene', 'Type':'Character', 'Number': '1'})
    w = Writer(args.outfile, vcf)

pvals = []

for i in range(len(df)):
    gene = df.index.values.tolist()[i]
    pos = int(df.iloc[i,2])
    chrom = df.iloc[i,1]
    # print('{}\t{}\t{}'.format(gene, chrom, pos))

    vcf_reader = find_closeby_svs(vcf, chrom, pos, pos, args.tolerance)

    for record in vcf_reader:
        gts = record.gt_types
        affected_samples = samples[np.logical_or(gts == 1, gts == 2)]
        unaffected_samples = samples[gts == 0]
        called_samples = samples[gts != 3]

        if len(affected_samples)/len(called_samples) > args.maf and len(affected_samples)/len(called_samples) < (1-args.maf): # allele frequency > 5%

            dat = df.iloc[i,3:].to_numpy()
            aff_dat = df.iloc[i][affected_samples].to_numpy()
            unaff_dat = df.iloc[i][unaffected_samples].to_numpy()
            aff_dat = aff_dat-np.mean(dat)
            unaff_dat = unaff_dat-np.mean(dat)
            aff_dat = aff_dat/np.std(dat)
            unaff_dat = unaff_dat/np.std(dat)

            _, p = mannwhitneyu(aff_dat, unaff_dat, alternative='two-sided')
            pvals.append(str(p))
            if p < args.p:
                # print('\tSV ID {} is significant with p={}'.format(record.ID, p))
                if args.outfile:
                    record.INFO['gene'] = gene
                    w.write_record(record)

    vcf_reader.close()



if args.outfile:
    w.close()

vcf.close()

with open(args.pfile, 'w') as savefile:
    print(','.join(pvals), file=savefile)
