# make_svsize_plot.py
# Makes a plot of the sizes of SVs
# by default uses length_statistics.dict, which is generated from process_vcf.py

import numpy as np
import pickle
import matplotlib.pyplot as plt

def makeplot(dict):
    dict['ME DEL'] = dict['DEL_ALU']
    dict['ME DEL'] += dict['DEL_LINE1']
    dict['ME DEL'] += dict['DEL_SVA']
    dict['ME DEL'] += dict['DEL_HERV']
    del dict['DEL_ALU']
    del dict['DEL_LINE1']
    del dict['DEL_SVA']
    del dict['DEL_HERV']

    dict['ME INS'] = dict['ALU']
    dict['ME INS'] += dict['LINE1']
    dict['ME INS'] += dict['SVA']
    del dict['ALU']
    del dict['LINE1']
    del dict['SVA']

    colors = {'DEL': 'blue', 'ME DEL': 'orange', 'DUP': 'red', 'INV': 'green', 'CNV': 'gray', 'INS': 'brown', 'ME INS': 'purple'}

    fig, ax = plt.subplots()
    logbins = np.logspace(np.log10(10**2),np.log10(10**7),100)
    for key in dict:
        list = dict[key]
        final_hist, final_bins = np.histogram(list, bins=logbins)
        bin_centers = 0.5*(final_bins[1:]+final_bins[:-1])
        plt.plot(bin_centers,final_hist, label=key, c=colors[key])
    ax.set_xscale('log')
    a=ax.get_xticks().tolist()
    a[2] = '100 bp'
    a[3] = '1 kb'
    a[4] = '10 kb'
    a[5] = '100 kb'
    a[6] = '1 Mb'
    a[7] = '10 Mb'
    ax.set_xticklabels(a)
    ax.legend()
    ax.set_yscale('symlog')
    plt.xlabel('Size')
    plt.ylabel('Number of SVs')
    plt.ylim([0, 10**4])
    plt.title('Size distribution of SVs')
    plt.show()

def main():
    dict = pickle.load(open('length_statistics.dict', 'rb'))
    makeplot(dict)

if(__name__ == '__main__'):
    main()
