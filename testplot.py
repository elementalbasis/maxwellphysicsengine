import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df = pd.read_csv('testdata.tsv', sep = '\t')

for N in set(df['N']):
    t_u = df[df['N'] == N]['t_u']
    L = df[df['N'] == N]['L']
    plt.plot(L, t_u)

plt.show()
