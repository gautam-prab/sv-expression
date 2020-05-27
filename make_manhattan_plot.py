# make_manhattan_plots.py
# make a manhattan plot from variants from process_pairs.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse

def build_args():
    parser = argparse.ArgumentParser(description='Plot a Manhattan plot')
    parser.add_argument('inputfile', type=str, help='CSV to make Manhattan plot of')
    parser.add_argument('--column', default='p', help='Name of column to make plot of')
    parser.add_argument('--highlight', type=int, nargs='*', help='Variant #s (ordering in the CSV, 0-based) to highlight in the plot')
    parser.add_argument('--title', type=str, default='P Values in CSV', help='Title of plot')
    parser.add_argument('--legend', nargs=2, type=str, default=['Special Variants', 'Others'], help='Legend description of highlighted variants, others')
    parser.add_argument('--poslog', action='store_true', help='Use if you want positive log instead of negative log')
    return parser.parse_args()

args = build_args()
df = pd.read_csv(args.inputfile, sep=',')
if args.highlight:
    special = df.loc[args.highlight]
    not_special = df.loc[~df.index.isin(args.highlight)]
    if args.poslog:
        plt.scatter(special.index, np.log10(special[args.column]), c='r', marker='D', label=args.legend[0])
        plt.scatter(not_special.index, np.log10(not_special[args.column]), c='b', label=args.legend[1])
        plt.ylabel('P (log10)')
    else:
        plt.scatter(special.index, -np.log10(special[args.column]), c='r', marker='D', label=args.legend[0])
        plt.scatter(not_special.index, -np.log10(not_special[args.column]), c='b', label=args.legend[1])
        plt.ylabel('P (-log10)')
    plt.legend()
else:
    if args.poslog:
        plt.scatter(df.index, np.log10(df[args.column]), c='b')
        plt.ylabel('P (log10)')
    else:
        plt.scatter(df.index, -np.log10(df[args.column]), c='b')
        plt.ylabel('P (-log10)')
plt.xlabel('Variant along Chromosome')
plt.title(args.title)
plt.show()
