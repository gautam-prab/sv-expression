# preprocess_vcf.py
# used to filter and create a VCF with per-genotype samples

from cyvcf2 import VCF, Writer
import pandas as pd
from snphwe import snphwe
import argparse

def build_args():
    parser = argparse.ArgumentParser(description='Filter vcf by \% called and HWE')
    parser.add_argument('vcffile', help='VCF to filter')
    parser.add_argument('--outfile', required=False, default='', help='Output for filtered VCF')
    parser.add_argument('--uncalled', default=0.4, type=float, help='Filter out variant with > this fraction of uncalled genotypes (default 0.4)')
    parser.add_argument('--hwe', default=0.0001, type=float, help='HWE p value to reject (default 0.0001)')
    parser.add_argument('--maf', default=10, type=int, help='Require at least this # of het or homAlt genotypes (default 10)')
    parser.add_argument('--ref', default=10, type=int, help='Require at least this # of homRef genotypes (default 10)')
    return parser.parse_args()

def make_vcfs(vcf_filenames, gt_file): # input a list of VCF filenames and a per-sample genotypes file
    genotypes = pd.read_csv(gt_file, sep='\t', index_col=0)
    filenames = vcf_filenames
    for filename in filenames:
        old_vcf = pd.read_csv(filename, sep='\t', comment='#', names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'N/A'])
        old_vcf = old_vcf.drop('N/A', axis='columns')
        old_vcf['FORMAT'] = 'GT'

        geno_interest = genotypes.loc[old_vcf['ID']].set_index(old_vcf.index)
        new_vcf = pd.concat([old_vcf, geno_interest], axis=1)
        new_vcf.to_csv(filename[:-4]+'.per_sample.vcf', index=False, sep='\t')
        # I manually pasted in the metadata from the VCFs


# this requires per-sample VCFs
def filter_vcfs(args):
    vcf_reader = VCF(args.vcffile, gts012=True) # sets gt: 0 = homRef, 1 = het, 2 = homAlt, 3 = unknown
    nsamples = len(vcf_reader.samples)
    if (args.outfile != ''):
        vcf_writer = Writer(args.outfile, vcf_reader)
    overall_count = 0
    uncalled_count = 0
    hwe_count = 0
    maf_count = 0
    for v in vcf_reader:
        overall_count += 1

        gts = v.gt_types
        ref = len(gts[gts == 0])
        het = len(gts[gts == 1])
        alt = len(gts[gts == 2])
        unc = len(gts[gts == 3])

        # check number uncalled
        if ((unc / nsamples) > args.uncalled): # too many uncalled
            uncalled_count += 1
            continue

        # check HW equilibrium
        called_N = het + alt + ref
        n_het = round(het * called_N / nsamples)
        n_alt = round(alt * called_N / nsamples)
        n_ref = called_N - n_het - n_alt
        p = snphwe(n_het, n_alt, n_ref)
        if (p < args.hwe): # does not follow HWE expectation
            hwe_count += 1
            continue

        # check number alt and ref
        if ((het + alt) < args.maf or ref < args.ref):
            maf_count += 1
            continue

        if (args.outfile != ''):
            vcf_writer.write_record(v)

    if (args.outfile != ''):
        vcf_writer.close()
    vcf_reader.close()

    print('Saw %d variants overall' % overall_count)
    print('Filtered out %d due to too many uncalled' % uncalled_count)
    print('Filtered out %d due to breaking HWE assumption' % hwe_count)
    print('Filtered out %d due to not enough Alt or Ref genotpyes' % maf_count)
    print('Left with %d variants' % (overall_count - uncalled_count - hwe_count - maf_count))


if __name__ == '__main__':
    args = build_args()
    filter_vcfs(args)
