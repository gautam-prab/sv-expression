# enrichment_permutation_test.py
# do a permutation test for the enrichment of variants to a specific feature
# feature file should be in .bed format (or --enhancers if using Dragon db enhancers)

import argparse
import numpy as np
import pandas as pd
from cyvcf2 import VCF
import matplotlib.pyplot as plt
import mygene
import ast
import pybedtools
from intervaltree import IntervalTree
from pyliftover import LiftOver
import os, glob

def build_args():
    parser = argparse.ArgumentParser(description='Do enrichment permutation tests')
    parser.add_argument('--rnafile', default='Data/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt')
    parser.add_argument('--vcffolder', default='quantile_vcfs/')
    parser.add_argument('--exons', default=None, nargs='?') # for method 1 exons
    parser.add_argument('--enhancers', default=False, const=True, nargs='?')
    parser.add_argument('--bedfile', default=None)
    parser.add_argument('--genedistance', default=1000000, type=int)
    parser.add_argument('--flankingdist', default=0, type=int)
    parser.add_argument('--genebed', default=None)
    return parser.parse_args()
args = build_args()

def build_interval_list(vcf, exon_df, tol=1000000):
    sv_list = []
    genes = []
    chrs = {}
    for record in vcf:
        sv_list.append(('chr'+record.CHROM, record.start, record.end))

        gene = record.INFO.get('gene')
        gene = gene.split('.')[0]
        if gene not in exon_df.index:
            print('could not find gene: {}'.format(gene))
            continue
        genes.append(gene)
        chrs[gene] = 'chr' + record.CHROM
    gene_list = []
    exon_list = []
    lo = LiftOver('hg38', 'hg19')
    for gene in genes:
        gene_list.append(get_gene_coords(exon_df.loc[gene], chrs[gene], tol, lo))
        exon_list.extend(get_exon_coords(exon_df.loc[gene], chrs[gene], lo))
    return sv_list, gene_list, exon_list

def get_gene_coords(gene, chrom, tol, lo):
    if pd.isna(gene['genomic_pos.start']) or pd.isna(gene['genomic_pos.end']) or pd.isna(gene['genomic_pos.strand']):
        start, end, strand = resolve_gene_coords(gene)
    else:
        start = gene['genomic_pos.start']
        end = gene['genomic_pos.end']
        strand = gene['genomic_pos.strand']
    gene_start = lo.convert_coordinate(chrom, int(start), strand)[0]
    gene_end = lo.convert_coordinate(chrom, int(end), strand)[0]
    if gene_start > gene_end: # reverse strand
        return (chrom, max(gene_end[1] - tol, 0), gene_start[1] + tol)
    else:
        return (chrom, max(gene_start[1] - tol, 0), gene_end[1] + tol)

def resolve_gene_coords(gene):
    dict_list = ast.literal_eval(gene['genomic_pos'])
    start = dict_list[0]['start']
    end = dict_list[0]['end']
    strand = dict_list[0]['strand']
    return start, end, strand

def get_exon_coords(gene, chrom, lo):
    exon_list = []
    if not pd.isna(gene['exons_hg19']): # has exons
        dict_list = ast.literal_eval(gene['exons_hg19'])
        transcript = dict_list[0] # just use first transcript to not get weird variants
        for exon in transcript['position']:
            if exon[0] > exon[1]:
                exon_list.append((chrom, max(exon[1], 0), exon[0]))
            else:
                exon_list.append((chrom, max(exon[0], 0), exon[1]))
    else: # no exons (lncRNA, etc)
        exon_list = [get_gene_coords(gene, chrom, 0, lo)]
    return exon_list

def set_gene_tolerance(f, tol):
    f.start = max(0, f.start - tol)
    f.stop = max(f.stop + tol, f.start + 1)
    return f

def build_sv_list(vcf):
    sv_list = []
    for record in vcf:
        sv_list.append(('chr'+record.CHROM, record.start, record.end))
    return sv_list

def read_enhancers():
    temp_df = pd.read_csv('enhancers.csv', usecols=['Chrom', 'Start', 'End'])
    return pybedtools.BedTool().from_dataframe(temp_df)


args = build_args()
if args.genebed is None:
    if args.genedistance == 1000000:
        genebed = 'gene_bedfiles/hg19_nongappedgenes_1MB.bed'
    else:
        genes = pybedtools.BedTool('gene_bedfiles/hg19_genes.bed').each(set_gene_tolerance, args.genedistance)
        gaps = pybedtools.BedTool('gene_bedfiles/hg19_gaps.bed')
        temp = genes.subtract(gaps) # nongapped genome
        temp.saveas('temp.bed')
        genebed = 'temp.bed'
else:
    genes = pybedtools.BedTool(args.genebed).each(set_gene_tolerance, args.genedistance)
    gaps = pybedtools.BedTool('gene_bedfiles/hg19_gaps.bed')
    temp = genes.subtract(gaps) # nongapped genome
    temp.saveas('temp.bed')
    genebed = 'temp.bed'

if (args.enhancers):
    enhancers = read_enhancers().each(set_gene_tolerance, args.flankingdist)
    enhancers.saveas('temp_bedfiles/enhancers_wflank.bed')
    featfile = 'temp_bedfiles/enhancers_wflank.bed'
elif (args.bedfile is not None):
    feats = pybedtools.BedTool(args.bedfile).each(set_gene_tolerance, args.flankingdist)
    feats.saveas('temp_bedfiles/features_wflank.bed')
    featfile = 'temp_bedfiles/features_wflank.bed'
elif (args.exons is not None):
    exon_df = pd.read_csv(args.exons, index_col = 0)

for file in sorted(glob.glob(args.vcffolder + '*.vcf.gz')):
    vcf = VCF(file)

    if (args.exons is not None):
        sv_list, gene_list, exon_list = build_interval_list(vcf, exon_df, args.genedistance)
        svs = pybedtools.BedTool(sv_list)
        genes = pybedtools.BedTool(gene_list)
        gaps = pybedtools.BedTool('gene_bedfiles/hg19_gaps.bed')
        temp = genes.subtract(gaps) # nongapped genome
        temp.saveas('temp.bed')
        genebed = 'temp.bed'
        exons = pybedtools.BedTool(exon_list)

        print('Doing permutation test for file {}...'.format(file))
        actual = len(svs.intersect(exons, u=True))
        print('Actual value: %d' % actual)
        outs = list(svs.randomintersection(exons, iterations=100, shuffle_kwargs={'genome': 'hg19', 'incl': genebed}, intersect_kwargs={'u': True}))
        results = np.array(outs)
        print('Permuted median: %d' % np.median(results))
        print('Percentile: %d' % len(results[results < actual]))
    else:
        sv_list = build_sv_list(vcf)
        svs = pybedtools.BedTool(sv_list)

        print('Doing permutation test for file {}...'.format(file))
        feats = pybedtools.BedTool(featfile)
        actual = len(svs.intersect(feats, u=True))
        print('Actual value: %d' % actual)
        outs = list(svs.randomintersection(feats, iterations=100, shuffle_kwargs={'genome': 'hg19', 'incl': genebed}, intersect_kwargs={'u': True}))
        results = np.array(outs)
        print('Permuted median: %d' % np.median(results))
        print('Percentile: %d' % len(results[results < actual]))
