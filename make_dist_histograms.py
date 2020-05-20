# make histograms of the distribution of +/-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from cyvcf2 import VCF, Writer
import pickle
from scipy.stats import mannwhitneyu, linregress

# find closeby structural variants to an ORF
def find_closeby_svs(vcf, chrom, start, end, tol=100000):
    return vcf('{}:{}-{}'.format(chrom, max(0, start-tol), end+tol))

def create_dataframe():
    df1 = pd.read_csv('Data/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt', sep='\t', index_col=0)

    # get rid of samples that aren't in both datasets
    pop_df = pd.read_csv('Data/integrated_call_samples_v3.20130502.ALL.panel', sep='\t', index_col=0)
    geuvadis_samples = df1.iloc[:,3:].keys().values
    svs_samples = pop_df.index.values
    samples = np.intersect1d(geuvadis_samples, svs_samples)
    exclude_samples = np.setdiff1d(geuvadis_samples, samples)
    del pop_df # don't need this anymore
    df1 = df1.drop(columns=exclude_samples) # get rid of extraneous samples
    return df1, samples

fromscratch = False

df, samples = create_dataframe()
fig = plt.figure()

if fromscratch:
    vcf = VCF('Data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz')
    sample_ord = vcf.samples

    bias = -1
    j = 1

    while j < 17:
        i = j + bias
        gene = df.index.values.tolist()[i]
        pos = int(df.iloc[i,2])
        chrom = df.iloc[i,1]
        print('{}\t{}\t{}'.format(gene, chrom, pos))
        if (chrom == 'X' or chrom == 'Y'):
            bias+=1
            continue

        vcf_reader = find_closeby_svs(vcf, chrom, pos, pos, tol=100000)
        flag = False

        for record in vcf_reader:
            affected_samples = np.array([])
            for sample in samples:
                gt = record.genotypes[sample_ord.index(sample)][0] + record.genotypes[sample_ord.index(sample)][1]
                if gt > 0:
                    affected_samples = np.append(affected_samples, sample)
            if len(affected_samples)/len(samples) > 0.05 and len(affected_samples)/len(samples) < 0.95: # allele frequency > 5%
                unaffected_samples = np.setdiff1d(samples, affected_samples)

                dat = df.iloc[i,3:].to_numpy()
                aff_dat = df.iloc[i][affected_samples].to_numpy()
                unaff_dat = df.iloc[i][unaffected_samples].to_numpy()
                aff_dat = aff_dat-np.mean(dat)
                unaff_dat = unaff_dat-np.mean(dat)
                aff_dat = aff_dat/np.std(dat)
                unaff_dat = unaff_dat/np.std(dat)

                _, p = mannwhitneyu(aff_dat, unaff_dat, alternative='two-sided')
                if p < (10 ** -40):
                    flag = True
                    print('\tSV ID {} is significant with p={}'.format(record.ID, p))
                    break
        vcf_reader.close()

        if not flag: # failed to find an SV
            bias += 1
            continue

        ax = fig.add_subplot(4, 4, j)
        ax.hist([aff_dat, unaff_dat], bins=50, stacked=True, label=['Positive', 'Negative'], color=['red', 'blue'])
        ax.set_title('Exp. for Gene {0} with SV {1}\np={2:.3E}'.format(i, record.INFO.get('SVTYPE'), p), fontsize=8)
        ax.set_xlabel('Normalized Expression', fontsize=8)
        ax.set_ylabel('Count', fontsize=8)
        ax.legend(prop={'size': 6})
        plt.xticks(fontsize=6)
        plt.yticks(fontsize=6)
        j += 1

else:
    vcf = VCF('cosmicVariants_filtered.vcf.gz')
    sample_ord = vcf.samples

    j = 1
    while j < 11:
        record = next(vcf)
        gene = record.INFO.get('gene')
        pos = int(df.loc[gene].iloc[2])
        chrom = df.loc[gene].iloc[1]
        print('{}\t{}\t{}'.format(gene, chrom, pos))

        flag = False
        affected_samples = np.array([])
        for sample in samples:
            gt = record.genotypes[sample_ord.index(sample)][0] + record.genotypes[sample_ord.index(sample)][1]
            if gt > 0:
                affected_samples = np.append(affected_samples, sample)
        if len(affected_samples)/len(samples) > 0.05 and len(affected_samples)/len(samples) < 0.95: # allele frequency > 5%
            unaffected_samples = np.setdiff1d(samples, affected_samples)

            dat = df.loc[gene].iloc[3:].to_numpy()
            aff_dat = df.loc[gene][affected_samples].to_numpy()
            unaff_dat = df.loc[gene][unaffected_samples].to_numpy()
            aff_dat = aff_dat-np.mean(dat)
            unaff_dat = unaff_dat-np.mean(dat)
            aff_dat = aff_dat/np.std(dat)
            unaff_dat = unaff_dat/np.std(dat)

            _, p = mannwhitneyu(aff_dat, unaff_dat, alternative='two-sided')
            if p < (6 * 10 ** -6):
                flag = True
                print('\tSV ID {} is significant with p={}'.format(record.ID, p))

        if not flag: # failed to find an SV
            continue

        ax = fig.add_subplot(2,5, j)
        bins = np.histogram(np.hstack([aff_dat, unaff_dat]), bins=50)[1]
        ax.hist(aff_dat, bins=bins, label='Positive', color='red', alpha=0.5)
        ax.hist(unaff_dat, bins=bins, label='Negative', color='blue', alpha=0.5)
        ax.set_title('Exp. for Gene {0}\nwith SV {1}\np={2:.3E}'.format(gene, record.ID, p), fontsize=8)
        ax.set_xlabel('Normalized Expression', fontsize=8)
        ax.set_ylabel('Count', fontsize=8)
        ax.legend(prop={'size': 6})
        plt.xticks(fontsize=6)
        plt.yticks(fontsize=6)

        if (p < 5 * 10 ** -57):
            out_unaff = unaff_dat
            out_aff = aff_dat
            most = pd.DataFrame(aff_dat, columns=['expr'], index=affected_samples, dtype=float)
            most['gt'] = 0
            for sample in samples:
                gt0 = record.genotypes[sample_ord.index(sample)][0]
                gt1 = record.genotypes[sample_ord.index(sample)][1]
                if gt0 > 1 or gt1 > 1:
                    print('found: {}'.format(most.loc[sample, 'expr']))
                    most.loc[sample, 'gt'] = 3
                else:
                    most.loc[sample, 'gt'] = gt0 + gt1
            most.to_csv('affected.csv')
        j += 1


plt.show()



# plt.boxplot([out_unaff, out_aff])
# out_y = np.concatenate((out_unaff, out_aff))
# out_x = np.concatenate((np.zeros(len(out_unaff)), np.ones(len(out_aff))))
# slope, intercept, r_value, p_value, std_err = linregress(out_x.astype(float), out_y.astype(float))
# line_x = [1, 2] # this is where the boxplots will be
# line_y = [intercept, slope+intercept]
# plt.plot(line_x,line_y)
# plt.xlabel('Genotype')
# plt.ylabel('Phenotype')
# plt.xticks([1,2], labels=['0', '1'])
# plt.title('Fitting a linear model')
# plt.show()

vcf.close()
