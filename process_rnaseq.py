# process RNA seq data with Pandas

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from cyvcf2 import VCF, Writer
import pickle
from scipy.stats import mannwhitneyu
from pyliftover import LiftOver

def build_args():
    parser = argparse.ArgumentParser(description='Find closeby structural variants to genes with differential expression')
    parser.add_argument('--p', default=0.05/128065, type=float, help='P-value cutoff') # 12816 for 100kb tolerance, 128065 for 1Mb
    parser.add_argument('--rnafile', default='Data/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt', help='TSV of RNA seq data for samples of interest')
    parser.add_argument('--vcffile', default='Data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz', help='VCF file containing genotypes for variants of interest')
    parser.add_argument('--outfile', const='SVs.vcf.gz', default=None, nargs='?', help='Output file, a VCF which will contain only the variants deemed significant')
    parser.add_argument('--tolerance', default=100000, type=int, help='cis region near genes to consider (in # of bases)')
    parser.add_argument('--maxchrom', default=None, help='set to 22 if you want to do autosomes only')
    parser.add_argument('--maf', default=0.05, type=float, help='Minor allele frequency cutoff; only consider those with higher MAF')
    parser.add_argument('--pfile', default='pvalues.csv', help='File to save p-values to for p-value cutoff analysis')
    parser.add_argument('--reference', default='hg19', help='Reference genome to use for VCF file. Enter either hg19 or hg38.')
    return parser.parse_args()

# find closeby structural variants to an ORF
def find_closeby_svs(vcf, chrom, start, end, tol=100000):
    return vcf('{}:{}-{}'.format(chrom, max(0, start-tol), end+tol))

def create_dataframe(args):
    df1 = pd.read_csv(args.rnafile, sep='\t', index_col=0)

    # get rid of samples that aren't in both datasets
    pop_df = pd.read_csv('Data/integrated_call_samples_v3.20130502.ALL.panel', sep='\t', index_col=0)
    geuvadis_samples = df1.iloc[:,3:].keys().values
    svs_samples = pop_df.index.values
    samples = np.intersect1d(geuvadis_samples, svs_samples)
    exclude_samples = np.setdiff1d(geuvadis_samples, samples)
    del pop_df # don't need this anymore
    df1 = df1.drop(columns=exclude_samples) # get rid of extraneous samples
    return df1, samples

args = build_args()

df, samples = create_dataframe(args)
if args.maxchrom:
    df = df[pd.to_numeric(df['Chr'], errors='coerce') <= int(args.maxchrom)]
    print(df)
fig = plt.figure()

vcf = VCF(args.vcffile)
sample_ord = vcf.samples
if args.outfile:
    vcf.add_info_to_header({'ID': 'gene', 'Description': 'overlapping gene', 'Type':'Character', 'Number': '1'})
    w = Writer(args.outfile, vcf)

pvals = []

if args.reference.lower() == 'hg38':
    lo = LiftOver('hg19', 'hg38')

for i in range(len(df)):
    gene = df.index.values.tolist()[i]
    pos = int(df.iloc[i,2])
    chrom = df.iloc[i,1]
    print('{}\t{}\t{}'.format(gene, chrom, pos))

    if args.reference.lower() == 'hg38':
        pos = lo.convert_coordinate('chr' + chrom, pos)
        if len(pos) > 0:
            pos = pos[0][1]
        else:
            continue

    vcf_reader = find_closeby_svs(vcf, chrom, pos, pos, tol=args.tolerance)
    flag = False

    for record in vcf_reader:
        affected_samples = np.array([])
        for sample in samples:
            gt = record.genotypes[sample_ord.index(sample)][0] + record.genotypes[sample_ord.index(sample)][1]
            if gt > 0:
                affected_samples = np.append(affected_samples, sample)
        if len(affected_samples)/len(samples) > args.maf and len(affected_samples)/len(samples) < (1-args.maf): # allele frequency > 5%
            unaffected_samples = np.setdiff1d(samples, affected_samples)

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
                print('\tSV ID {} is significant with p={}'.format(record.ID, p))
                if args.outfile:
                    record.INFO['gene'] = gene
                    w.write_record(record)

    vcf_reader.close()



if args.outfile:
    w.close()

vcf.close()

with open(args.pfile, 'w') as savefile:
    print(','.join(pvals), file=savefile)
