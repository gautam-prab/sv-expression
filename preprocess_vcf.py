# preprocess_vcf.py
# used to filter and create a VCF with per-genotype samples

from cyvcf2 import VCF, Writer
import pandas as pd
from snphwe import snphwe
import argparse

def build_args():
    parser = argparse.ArgumentParser(description='Filter vcf by % called and HWE')
    parser.add_argument('--vcffile', required=True, help='VCF to filter')
    parser.add_argument('--statsfile', required=True, help='Genotype Frequencies')
    parser.add_argument('--outfile', required=False, default='', help='Output for filtered VCF')
    parser.add_argument('--uncalled', default=0.4, type=float, help='Fraction uncalled required to reject (default 0.4)')
    parser.add_argument('--hwe', default=0.0001, type=float, help='HWE p value to reject (default 0.0001)')
    parser.add_argument('--nsamples', default=2504, type=int, help='# of samples in VCF')
    return parser.parse_args()

def make_vcfs():
    genotypes = pd.read_csv('Data/cancer_data/cosmicOrganoid_SV_genotypes_1KGP_all.txt', sep='\t', index_col=0)
    filenames = ['48T_refined_cosmicVariants_allPlatforms_paragraphFormat_Nov12.original.vcf', '51T_refined_cosmicVariants_allPlatforms_paragraphFormat_Nov12.original.vcf', 'SKBR3_refined_cosmicVariants_allPlatforms_paragraphFormat_Nov12.original.vcf']
    for filename in filenames:
        old_vcf = pd.read_csv('Data/cancer_data/'+filename, sep='\t', comment='#', names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'N/A'])
        old_vcf = old_vcf.drop('N/A', axis='columns')
        old_vcf['FORMAT'] = 'GT'

        geno_interest = genotypes.loc[old_vcf['ID']].set_index(old_vcf.index)
        new_vcf = pd.concat([old_vcf, geno_interest], axis=1)
        new_vcf.to_csv('Data/cancer_data/'+filename[:-4]+'.per_sample.vcf', index=False, sep='\t')
        # I manually pasted in the metadata from the VCFs


def filter_vcfs(args):
    N = args.nsamples
    vcf_reader = VCF(args.vcffile)
    if (args.outfile != ''):
        vcf_writer = Writer(args.outfile, vcf_reader)
    stats = pd.read_csv(args.statsfile, sep='\t', index_col=0)
    overall_count = 0
    uncalled_count = 0
    hwe_count = 0
    for v in vcf_reader:
        overall_count += 1
        # check % uncalled
        unc = stats.loc[v.ID, 'uncalled']
        if (unc > args.uncalled): # too many uncalled
            uncalled_count += 1
            continue

        # check HW equilibrium
        het = stats.loc[v.ID, 'het']
        alt = stats.loc[v.ID, 'homAlt']
        ref = stats.loc[v.ID, 'homRef']
        called_N = round((het + alt + ref) * N)
        n_het = round(het * called_N)
        n_alt = round(alt * called_N)
        n_ref = called_N - n_het - n_alt
        p = snphwe(n_het, n_alt, n_ref)
        if (p < args.hwe): # does not follow HWE expectation
            hwe_count += 1
            continue

        if (args.outfile != ''):
            vcf_writer.write_record(v)

    if (args.outfile != ''):
        vcf_writer.close()
    vcf_reader.close()

    print('Saw %d variants overall' % overall_count)
    print('Filtered out %d due to too many uncalled' % uncalled_count)
    print('Filtered out %d due to breaking HWE assumption' % hwe_count)
    print('Left with %d variants' % (overall_count - uncalled_count - hwe_count))


if __name__ == '__main__':
    args = build_args()
    filter_vcfs(args)
