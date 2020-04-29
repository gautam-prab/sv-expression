import matplotlib.pyplot as plt
import numpy as np

quantile_names = ['(0,0.5]','(0.5,0.6]','(0.6,0.7]','(0.7,0.8]','(0.8,0.9]','(0.9,1.0]']
enrichment = np.array([79, 23, 17, 31, 23, 50])
null = np.array([25, 2, 2, 6, 4, 17])
enrichment = enrichment / null

plt.plot(quantile_names, enrichment)
plt.xlabel('Quantile')
plt.ylabel('Enrichment within Exons')
plt.title('Enrichment within Exons of eGenes')
plt.show()
