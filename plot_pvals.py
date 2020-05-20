# make plot of p values and calculate p-value thresholds

import numpy as np
import matplotlib.pyplot as plt
import argparse

def build_args():
    parser = argparse.ArgumentParser(description='Plot pvals from process_rnaseq.py')
    parser.add_argument('filename', metavar='filename.csv', type=str, help='The csv file containing pvalues')
    return parser.parse_args()

args = build_args()
pvals = np.loadtxt(args.filename, delimiter=',')
bh_sort = np.sort(pvals)
valid_k = np.where(bh_sort <= np.arange(1,len(pvals)+1)*0.05/len(pvals))[0]
if len(valid_k) != 0:
    bh_kstar = valid_k[bh_sort[valid_k].argmax()] + 1
else:
    bh_kstar = 0
bhp_1 = 0.05*bh_kstar/len(pvals)
valid_k = np.where(bh_sort <= np.arange(1,len(pvals)+1)*0.1/len(pvals))[0]
if len(valid_k) != 0:
    bh_kstar = valid_k[bh_sort[valid_k].argmax()] + 1
else:
    bh_kstar = 0
bhp_2 = 0.1*bh_kstar/len(pvals)

print('You have %d pvals overall' % len(pvals))
print('Bonferroni p: %e' % (0.05/len(pvals)))
print('\tCaptures %d pairs' % len(pvals[pvals < 0.05/len(pvals)]))
print('BH 0.05 p: %e' % bhp_1)
print('\tCaptures %d pairs' % len(pvals[pvals < bhp_1]))
print('BH 0.1 p: %e' % bhp_2)
print('\tCaptures %d pairs' % len(pvals[pvals < bhp_2]))

hist, bins, _ = plt.hist(pvals, bins=50)
plt.xlabel('P value')
plt.ylabel('Frequency')
plt.title('P values')

plt.show()

ax = plt.gca()
logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
plt.hist(pvals, bins=logbins)
plt.xlabel('P value')
plt.axvline(0.05/len(pvals), color='red',linestyle='dashed')
plt.text(x=0.05/len(pvals), y=0.96, s='Bonferroni', fontsize=8, transform=ax.get_xaxis_transform())
plt.axvline(bhp_2, color='red',linestyle='dashed')
plt.text(x=bhp_2, y=0.96, s='BH FDR = 0.1', fontsize=8, transform=ax.get_xaxis_transform())
plt.xscale('log')
plt.ylabel('Frequency')
plt.title('P values (log scale)')
plt.show()
