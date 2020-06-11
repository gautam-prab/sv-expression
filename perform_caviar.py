import subprocess
import argparse
import pandas as pd
import numpy as np
import re

parser = argparse.ArgumentParser(description='Run CAVIAR given SNP file and SV file')
parser.add_argument('--snpfile', type=str, required=True, help='file with SNPs (.vcf.gz)')
parser.add_argument('--svfile', type=str, required=True, help='file with SVs (.vcf.gz)')
parser.add_argument('--rnafile', type=str, required=True, help='file with gene RNA-seq (.txt)')
parser.add_argument('--gene', type=str, required=True, help='ENSEMBL ID to associate')
parser.add_argument('--genename', type=str, help='Standard Gene Name')
parser.add_argument('--plinkpath', type=str, default='binaries/PLINK', help='path to PLINK binary')
parser.add_argument('--caviarpath', type=str, default='binaries/CAVIAR', help='path to CAVIAR binary')
parser.add_argument('--outfolder', type=str, required=True, help='folder to put output in')
args = parser.parse_args()

print('Merging variant files...')
with open(args.outfolder + 'all_variants.vcf', 'w') as file:
    subprocess.run(['bcftools', 'concat', args.snpfile, args.svfile], stdout=file)
subprocess.run(['bgzip', args.outfolder + 'all_variants.vcf'])
with open(args.outfolder + 'all_variants.vcf', 'w') as file:
    subprocess.run(['bcftools', 'sort', args.outfolder + 'all_variants.vcf.gz'], stdout=file)
subprocess.run(['rm', args.outfolder + 'all_variants.vcf.gz'])
subprocess.run(['bgzip', args.outfolder + 'all_variants.vcf'])
subprocess.run(['tabix', '-p', 'vcf', args.outfolder + 'all_variants.vcf.gz'])

print('Calculating p values...')
subprocess.run(['python', 'process_pairs.py', args.outfolder + 'all_variants.vcf.gz', '--gene', args.gene, '--outcsv', args.outfolder + 'all_pairs.csv'])
all_df = pd.read_csv(args.outfolder + 'all_pairs.csv')
all_df = all_df[all_df['p'] != 0]
sorted_df = all_df.sort_values(by='p')
p_val = sorted_df.iloc[1000]['p']
p_val = p_val*1.000001

print('Creating .Z file...')
subprocess.run(['python', 'process_rnaseq.py', '--vcffile', args.outfolder + 'all_variants.vcf.gz', '--rnafile', args.rnafile, '--p', str(p_val), '--outfile', args.outfolder + 'sig_variants.vcf.gz'])
with open(args.outfolder + 'sig_variants.vcf', 'w') as file:
    subprocess.run(['bcftools', 'sort', args.outfolder + 'sig_variants.vcf.gz'], stdout=file)
subprocess.run(['rm', args.outfolder + 'sig_variants.vcf.gz'])
subprocess.run(['bgzip', args.outfolder + 'sig_variants.vcf'])

subprocess.run(['python', 'process_pairs.py', args.outfolder + 'sig_variants.vcf.gz', '--gene', args.gene, '--outcsv', args.outfolder + 'sig_pairs.csv'])
sig_df = pd.read_csv(args.outfolder + 'sig_pairs.csv')
pd.Series(np.array(sig_df['beta'] / sig_df['std err']), index=np.array(sig_df['SV ID'])).to_csv(args.outfolder + 'Z.txt', sep = '\t', header = False)
sv_indices = [str(i) for i,item in enumerate(np.array(sig_df['SV ID'])) if (not item.startswith('var') and not item.startswith('snp'))]

print('Creating .ld file...')
subprocess.run(['gunzip', args.outfolder + 'sig_variants.vcf.gz'])
# edit sig_variants.vcf since vcftools won't like it
with open(args.outfolder + 'sig_variants.vcf', 'r') as infile, open(args.outfolder + 'sig_variants_vcftools.vcf', 'w') as outfile:
    for line in infile:
            if '##fileformat' in line:
                outfile.write('##fileformat=VCFv4.1\n')
            else:
                outfile.write(re.sub('".*?"', '""', line))
subprocess.run('ulimit -n 3000', shell=True)
subprocess.run('ulimit -n 3000 && vcftools --vcf {} --plink --out {}'.format(args.outfolder + 'sig_variants_vcftools.vcf', args.outfolder + 'plink'), shell=True)
subprocess.run([args.plinkpath, '--file', args.outfolder + 'plink', '--r', '--matrix', '--out', args.outfolder + 'plink_out'])
subprocess.run(['mv', args.outfolder + 'plink_out.ld', args.outfolder + 'ld.txt'])

print('Running CAVIAR...')
subprocess.run([args.caviarpath, '-z', args.outfolder + 'Z.txt', '-l', args.outfolder + 'ld.txt', '-o', args.outfolder + 'caviar_out'])
if not args.genename:
    name = args.gene
else:
    name = args.genename
subprocess.run(['python', 'make_manhattan_plot.py', args.outfolder + 'sig_pairs.csv', '--column', 'p', '--title', 'Association P Value for %s' % name, '--legend', 'SV', 'SNP', '--highlight'] + sv_indices)
pd.read_csv(args.outfolder + 'caviar_out_post', sep='\t').to_csv(args.outfolder + 'caviar_out_post.csv', sep=',')
subprocess.run(['python', 'make_manhattan_plot.py', args.outfolder + 'caviar_out_post.csv', '--column', 'Causal_Post._Prob.', '--title', 'CAVIAR Posterior Probability for %s' % name, '--legend', 'SV', 'SNP', '--poslog', '--highlight'] + sv_indices)
