# process_vcf.py
# Reads VCF of structural variant calls into pandas for basic statistics

from cyvcf2 import VCF
import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt
import argparse

def build_args():
    parser = argparse.ArgumentParser(description='Find basic summary statistics for samples in a vcf')
    parser.add_argument('--vcffile', default='Data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz')
    parser.add_argument('--lengthfile', default='length_statistics.dict') #output file
    parser.add_argument('--samplefile', default='samples_counts.csv') #output file
    return parser.parse_args()

def countVariants(vcf_reader):
    '''this code counts the overall number of variants of each type'''
    dict = {}
    last_chrom = ''
    for record in vcf_reader:
        if record.CHROM != last_chrom:
            last_chrom = record.CHROM
            print('Chromosome #'+str(last_chrom))
        type = record.INFO['SVTYPE']
        if type in dict:
            dict[type] += 1
        else:
            dict[type] = 1

    return dict

def iterateSamples(vcf_reader, sv_types, outfile='samples_counts.csv'):
    '''this creates a table of genotypes for each sample'''
    samples = vcf_reader.samples
    sample_df = pd.DataFrame(index=samples, columns=sv_types)
    sample_df = sample_df.fillna(0) # fill with 0's

    last_chrom = '' # for printouts
    for record in vcf_reader:

        if record.CHROM != last_chrom:
            last_chrom = record.CHROM
            print('Chromosome #'+str(last_chrom))

        type = record.INFO['SVTYPE']
        for sample in samples:
            if record.genotypes[samples.index(sample)][0] + record.genotypes[samples.index(sample)][1] > 0:
                sample_df.at[sample, type] += 1

    sample_df.to_csv(outfile)

def lengthStatistics(vcf_reader, outfile='length_statistics.dict'):
    dict = {}
    last_chrom = ''
    for record in vcf_reader:
        if record.CHROM != last_chrom:
            last_chrom = record.CHROM
            print('Chromosome #'+str(last_chrom))
        type = record.INFO['SVTYPE']
        if record.INFO['SVLEN'] is not None:
            len = record.INFO['SVLEN']
            if len < 0:
                len = -len
        else:
            len = record.end - record.start

        if type not in dict:
            dict[type] = [len]
        else:
            dict[type].append(len)

        with open(outfile, 'wb') as file:
            pickle.dump(dict, file)

def main():
    args = build_args()
    count_dict = countVariants(VCF(args.vcffile))
    print(count_dict)
    lengthStatistics(VCF(args.vcffile), outfile=args.lengthfile)
    iterateSamples(VCF(args.vcffile), list(count_dict.keys()), outfile=args.samplefile)


if(__name__ == '__main__'):
    main()
