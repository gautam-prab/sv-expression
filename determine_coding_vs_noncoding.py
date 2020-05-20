# determine_coding_vs_noncoding.py
# read a VCF and use mygene to determine exonic regions
# uses VCF from process_vcf.py with gene identification

import argparse
import numpy as np
import pandas as pd
from cyvcf2 import VCF, Writer
import mygene
from pyliftover import LiftOver
from intervaltree import IntervalTree, Interval
import ast

def build_args():
    parser = argparse.ArgumentParser(description='Read a VCF and use mygene to determine exonic regions')
    parser.add_argument('--vcffile', default='allSVs.vcf.gz')
    parser.add_argument('--fromscratch', default=False, const=True, nargs='?') # query genes to get exons from mygene
    parser.add_argument('--exonfile', default='exons.csv')
    parser.add_argument('--outvcf', default=None)
    return parser.parse_args()

def get_samples(args):
    rna_df = pd.read_csv('Data/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt', sep='\t', index_col=0)
    pop_df = pd.read_csv('Data/integrated_call_samples_v3.20130502.ALL.panel', sep='\t', index_col=0)
    geuvadis_samples = rna_df.iloc[:,3:].keys().values
    svs_samples = pop_df.index.values
    samples = np.intersect1d(geuvadis_samples, svs_samples)
    del rna_df
    del pop_df
    return samples

# check if SV in record is within the exons in gene
def check_exonic(args, record, gene, lo):
    chrom = record.CHROM
    sv_start = record.start
    sv_end = record.end

    if np.array(pd.notna(gene['exons_hg19'])).any():
        t = IntervalTree()
        if args.fromscratch:
            dict_list = gene['exons_hg19']
        else:
            dict_list = ast.literal_eval(gene['exons_hg19'])
        for transcript in dict_list:
            count = 0
            for exon in transcript['position']:
                count +=1
                if exon[0] > exon[1]:
                    t[exon[1]:exon[0]] = count
                else:
                    t[exon[0]:exon[1]] = count
        if len(t[sv_start:sv_end]) > 0: # exonic
            return True
        else:
            return False
    else:
        gene_start = lo.convert_coordinate('chr' + chrom, int(gene['genomic_pos.start']), gene['genomic_pos.strand'])[0]
        gene_end = lo.convert_coordinate('chr' + chrom, int(gene['genomic_pos.end']), gene['genomic_pos.strand'])[0]
        t = IntervalTree()
        if gene_start[1] > gene_end[1]:
            t[gene_end[1]:gene_start[1]] = 'a'
        else:
            t[gene_start[1]:gene_end[1]] = 'a'
        return len(t[sv_start:sv_end]) > 0

# when a gene doesn't say if it's protein-coding, determine that
def fix_df(args, gene):
    if args.fromscratch:
        dict_list = gene['ensembl']
    else:
        dict_list = ast.literal_eval(gene['ensembl'])
    gene['ensembl.type_of_gene'] = dict_list[0]['type_of_gene']
    if args.fromscratch:
        dict_list = gene['genomic_pos']
    else:
        dict_list = ast.literal_eval(gene['genomic_pos'])
    gene['genomic_pos.start'] = dict_list[0]['start']
    gene['genomic_pos.end'] = dict_list[0]['end']
    gene['genomic_pos.strand'] = dict_list[0]['strand']
    return gene

args = build_args()

keep_samples = get_samples(args)
vcf = VCF(args.vcffile, samples=keep_samples.tolist())
sample_order = vcf.samples

lo = LiftOver('hg38', 'hg19')

# get all unique gene names in VCF
list = []
for record in vcf:
    gene = record.INFO.get('gene')
    gene = gene.split('.')[0]
    if gene not in list:
        list.append(gene)
vcf.close()

# query each gene name to MyGene.info to get exon positions and type
if args.fromscratch:
    mg = mygene.MyGeneInfo()
    df = mg.getgenes(list, fields='genomic_pos,ensembl.type_of_gene,exons_hg19.position,exons.position', as_dataframe=True)
    print('Found {} out of {} genes'.format(len(df[df['notfound'].isna()]),len(df)))
    print(df['ensembl.type_of_gene'].value_counts())
else:
    df = pd.read_csv(args.exonfile, index_col = 0)

# create VCF reader and writer
vcf = VCF(args.vcffile, samples=keep_samples.tolist())
if args.outvcf:
    vcf.add_info_to_header({'ID': 'exonic', 'Description': 'is it exonic?', 'Type':'Character', 'Number': '1'})
    w = Writer(args.outvcf, vcf)

# check exon for each record
for record in vcf:
    gene = record.INFO.get('gene')
    gene = gene.split('.')[0]
    if df.loc[gene,'notfound'] == True:
        continue
    if pd.isna(df.loc[gene,'ensembl.type_of_gene']) or pd.isna(df.loc[gene,'genomic_pos.start']):
        df.loc[gene] = fix_df(args, df.loc[gene])

    if args.outvcf:
        if check_exonic(args, record, df.loc[gene], lo):
            record.INFO['exonic'] = 'y'
        else:
            record.INFO['exonic'] = 'n'
        w.write_record(record)
    else:
        if check_exonic(args, record, df.loc[gene], lo):
            print('{} is exonic'.format(record.ID))
        else:
            print('{} is not exonic'.format(record.ID))

vcf.close()

if args.fromscratch:
    df.to_csv(args.exonfile)

if args.outvcf:
    w.close()
