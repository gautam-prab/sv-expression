# make_svpopulation_plot.py
# Makes a plot of each sample and their SVs, separated by superpopulation
# Uses samples_counts.csv, generated from process_vcf.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def makeFigure1b(df, sort_by):
    '''this creates a plot of kinds of structural variants per sample'''
    # first consolidate some datatypes into others
    df['ME DEL'] = df['DEL_ALU']
    df['ME DEL'] += df['DEL_LINE1']
    df['ME DEL'] += df['DEL_SVA']
    df['ME DEL'] += df['DEL_HERV']
    df = df.drop(columns=['DEL_ALU', 'DEL_LINE1', 'DEL_SVA', 'DEL_HERV'])

    df['ME INS'] = df['ALU']
    df['ME INS'] += df['LINE1']
    df['ME INS'] += df['SVA']
    df = df.drop(columns=['ALU', 'LINE1', 'SVA'])

    # now we read superpopulation data
    pop_df = pd.read_csv('Data/integrated_call_samples_v3.20130502.ALL.panel', sep='\t', index_col=0)
    df[sort_by] = pop_df[sort_by]
    df = df.sort_values(by=sort_by, axis=0)

    # plot each sample
    colors = {'DEL': 'blue', 'ME DEL': 'orange', 'DUP': 'red', 'INV': 'green', 'CNV': 'gray', 'INS': 'brown', 'ME INS': 'magenta'}

    fig, ax = plt.subplots()
    for (name, data) in df.iteritems():
        if(name == 'DEL'):
            continue # ensure DELs get plotted last
        if(name == sort_by):
            continue
        ax.scatter(df.index, data, label=name, c=colors[name], s=1)
    name = 'DEL'
    ax.scatter(df.index, df[name], label=name, c=colors[name], s=1)
    ax.legend(loc='upper right')
    plt.xlabel('Donor')
    plt.ylabel('Number of SVs per Donor')
    plt.title('Number of SVs detected by Variant Type (Autosomes)')
    ax.set_xticklabels(['']*len(df))

    # draw vertical lines for populations
    first_person = df.iloc[0].index[0]
    for pop in df[sort_by].unique():
        first_person = df[df[sort_by] == pop].index[0]
        last_person = df[df[sort_by] == pop][-1:].index[0]
        ax.axvline(first_person, color='gray', linestyle='dashed', linewidth=0.5)
        plt.text(x=first_person,y=1750,s=pop,fontsize=8,color='gray')


    plt.show()



def main():
    print('Reading csv file...')
    csv_path = 'samples_counts_nonsex.csv'
    df = pd.read_csv(csv_path, index_col=0)
    makeFigure1b(df, 'super_pop')

if(__name__ == '__main__'):
    main()
