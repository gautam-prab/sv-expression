# process_vcf.py
# Reads VCF of structural variant calls into pandas for basic statistics

import vcf
import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt

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

    print(dict)
    print('Done!')

def iterateSamples(vcf_reader, outfile='samples_counts.csv'):
    '''this creates a table of genotypes for each sample'''
    sv_types = ['ALU', 'LINE1', 'SVA', 'INS', 'DEL', 'DEL_ALU', 'DEL_LINE1', 'DEL_SVA', 'DEL_HERV', 'DUP', 'INV', 'CNV']
    sample_df = pd.DataFrame(index=vcf_reader.samples, columns=sv_types)
    sample_df = sample_df.fillna(0) # fill with 0's

    last_chrom = '' # for printouts
    for record in vcf_reader:
        if record.CHROM == 'X':
            break

        if record.CHROM != last_chrom:
            last_chrom = record.CHROM
            print('Chromosome #'+str(last_chrom))

        type = record.INFO['SVTYPE']
        for sample in record.samples:
            gt = sample['GT']
            if gt != '0|0' and gt != '.' and gt != '0':
                sample_df.at[sample.sample, type] += 1

    sample_df.to_csv(outfile)

def lengthStatistics(vcf_reader):
    dict = {}
    last_chrom = ''
    for record in vcf_reader:
        if record.CHROM != last_chrom:
            last_chrom = record.CHROM
            print('Chromosome #'+str(last_chrom))
        type = record.INFO['SVTYPE']
        if 'SVLEN' in record.INFO:
            len = int(record.INFO['SVLEN'][0])
        else:
            len = int(record.INFO['END'] - record.POS)

        if type not in dict:
            dict[type] = [len]
        else:
            dict[type].append(len)

        with open('length_statistics.dict', 'wb') as file:
            pickle.dump(dict, file)

def mainLength():
    print('Reading vcf file...')
    vcf_path = 'Data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf'
    vcf_reader = vcf.Reader(open(vcf_path, 'r'))
    iterateSamples(vcf_reader, outfile='samples_counts_nonsex.csv')

if(__name__ == '__main__'):
    mainLength()
