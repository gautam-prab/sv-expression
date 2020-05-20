import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cyvcf2 import VCF
from scipy.stats import mannwhitneyu, linregress
import argparse

def build_args():
    parser = argparse.ArgumentParser(description='Plot the genotype-phenotype regression')
    parser.add_argument('inputfile', type=str, help='VCF (with gene annotations) of variants to make boxplots of')
    parser.add_argument('--rnafile', default='Data/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt')
    return parser.parse_args()

def resolve_alts(gt, alts, alt_dict):
    ret = 0
    if int(gt[0]) != 0: # not reference allele
        alt = alts[int(gt[0]) - 1]
        ret += alt_dict[alt]
    else: # reference allele is copy number 1
        ret += 1
    if int(gt[1]) != 0: # not reference allele
        alt = alts[int(gt[1]) - 1]
        ret += alt_dict[alt]
    else:
        ret += 1
    return ret

def main(args):
    rna_df = pd.read_csv(args.rnafile, sep='\t', index_col=0)
    rna_samples = rna_df.iloc[:,3:].keys().values

    vcf_reader = VCF(args.inputfile)
    vcf_samples = vcf_reader.samples

    samples = np.intersect1d(rna_samples, vcf_samples)

    alt_dict = {'<CN0>' : 0, '<CN2>' : 2, '<CN3>' : 3, '<CN4>' : 4, '<CN5>' : 5, '<CN6>' : 6, '<CN7>' : 7, '<CN8>' : 8, '<CN9>' : 9}

    for v in vcf_reader:
        gene = v.INFO.get('gene')
        phenotypes = rna_df.loc[gene].iloc[3:]
        phenotypes = (phenotypes - np.mean(phenotypes)) / np.std(phenotypes)

        gts = np.zeros(len(samples))
        phens = np.zeros(len(samples))

        alts = v.ALT
        for idx, sample in enumerate(samples):
            gt = v.genotypes[vcf_samples.index(sample)]
            if v.INFO.get('SVTYPE') == 'CNV' and alts[0] in alt_dict:
                gts[idx] = resolve_alts(gt, alts, alt_dict) # resolve copy numbers
            else:
                gts[idx] = gt[0] + gt[1]
            phens[idx] = phenotypes[sample]

        data = []
        copy_nums = np.unique(gts)
        cp_dict = {}
        for idx, num in enumerate(copy_nums):
            cp_dict[num] = idx + 1
            data.append(phens[gts == num])
        plt.boxplot(data, showfliers=False)
        line_x = np.arange(1, len(cp_dict) + 1) # this is where the boxplots will be
        gts_as_x = [cp_dict[gt] for gt in gts.astype(int)] # as x coordinates
        slope, intercept, r_value, p_value, std_err = linregress(gts_as_x, phens)
        line_y = line_x * slope + intercept
        plt.plot(line_x,line_y,linestyle='--')

        # add datapoints with 'jitter'
        gts_jitter = np.random.normal(np.zeros(len(gts)), 0.04)
        plt.scatter(gts_as_x + gts_jitter, phens, c=gts_as_x, cmap='rainbow', marker='.', alpha=0.4)

        plt.xlabel('Genotype/Copy Number')
        plt.ylabel('Phenotype')
        plt.xticks(line_x, labels=copy_nums)
        plt.title('Effect of {} at chr{}:{}-{} on Expression of {}'.format(v.INFO['SVTYPE'], v.CHROM, v.start, v.end, gene))
        plt.show()

if __name__ == '__main__':
    args = build_args()
    main(args)
