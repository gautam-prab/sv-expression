import pandas as pd
import numpy as np
from os import listdir
from os.path import isfile, join

df = pd.read_csv('samples_counts.csv', index_col = 0)
print(df.index)

onlyfiles = [f for f in listdir('Data/GSE92521_RAW/') if isfile(join('Data/GSE92521_RAW/', f))]
for f in onlyfiles:
    if f[0] == '.':
        continue
    name = f[11:18]
    if name in df:
        print(name)
