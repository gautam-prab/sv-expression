from cyvcf2 import VCF
import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Generate a plot of MAFs')
parser.add_argument('inputfile', help='Input file as a VCF')
args = parser.parse_args()

vcf_reader = VCF(args.inputfile, gts012=True)
samples = vcf_reader.samples

mafs = []
for v in vcf_reader:
    gts = v.gt_types
    maf = len(gts[np.logical_or(gts == 1, gts == 2)]) / len(gts[gts != 3])
    if maf > 0.5:
        maf = 1-maf
    if maf < 0.05:
        continue
    mafs.append(maf)

n, x = np.histogram(mafs, bins=10)
bin_centers = 0.5*(x[1:]+x[:-1])
plt.plot(bin_centers, n)
plt.xlabel('Minor Allele Frequency')
plt.xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
plt.ylabel('Variant count')
plt.show()
