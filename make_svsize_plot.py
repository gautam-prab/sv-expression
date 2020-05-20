# make_svsize_plot.py
# Makes a plot of the sizes of SVs
# by default uses length_statistics.dict, which is generated from process_vcf.py

import numpy as np
import pickle
import matplotlib.pyplot as plt
import argparse

def build_args():
    parser = argparse.ArgumentParser(description='Generate a plot of the sizes of SVs')
    parser.add_argument('--dictfile', default='length_statistics.dict')
    parser.add_argument('--mergeMEs', action='store_true', help='Merge Mobile Elements (from 1KGP SV calls)')
    parser.add_argument('--matchcolors', action='store_true', help='Match colors of GTEx nature genetics paper')
    return parser.parse_args()

def makeplot(dict, args):
    if args.mergeMEs:
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

    if args.matchcolors:
        colors = {'DEL': 'blue', 'ME DEL': 'orange', 'DUP': 'red', 'INV': 'green', 'CNV': 'gray', 'INS': 'brown', 'ME INS': 'purple'}
    else:
        col_keys = list(dict)
        col_vals = plt.cm.gist_rainbow(np.linspace(0,1,len(col_keys)))
        colors = {col_keys[i]: col_vals[i] for i in range(len(col_keys))}

    fig, ax = plt.subplots()
    logbins = np.logspace(np.log10(10**2),np.log10(10**10),100)
    for key in dict:
        sizes = dict[key]
        final_hist, final_bins = np.histogram(sizes, bins=logbins)
        bin_centers = 0.5*(final_bins[1:]+final_bins[:-1])
        plt.plot(bin_centers,final_hist, label=key, c=colors[key])
    ax.set_xscale('log')
    a=ax.get_xticks().tolist()
    a[2] = '1 kb'
    a[3] = '100 kb'
    a[4] = '10 Mb'
    a[5] = '1 Gb'
    ax.set_xticklabels(a)
    ax.legend()
    ax.set_yscale('symlog')
    plt.xlabel('Size')
    plt.ylabel('Number of SVs')
    plt.ylim([0, 10**4])
    plt.title('Size distribution of SVs')
    plt.show()

def main():
    args = build_args()
    dict = pickle.load(open(args.dictfile, 'rb'))
    makeplot(dict, args)

if(__name__ == '__main__'):
    main()
