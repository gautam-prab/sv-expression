# make plot of p values and calculate p-value thresholds

import numpy as np
import matplotlib.pyplot as plt

pvals = np.loadtxt('pvalues_1MB.csv', delimiter=',')
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

logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
plt.hist(pvals, bins=logbins)
plt.xlabel('P value')
plt.axvline(0.05/len(pvals), color='red',linestyle='dashed')
plt.text(x=0.05/len(pvals), y=10000, s='Bonferroni', fontsize=8)
plt.axvline(bhp_1, color='red',linestyle='dashed')
plt.text(x=bhp_1, y=9000, s='BH FDR = 0.05', fontsize=8)
plt.xscale('log')
plt.ylabel('Frequency')
plt.title('P values (log scale)')
plt.show()
