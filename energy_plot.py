import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

filename = 'out.tsv'
df = pd.read_csv(filename, sep = '\t')
m = 1e-3
N = 3000
times = set(df['t'])

t_array = []
K_array = []
for t in times:
    vx = df[df['t'] == t]['vx']
    vy = df[df['t'] == t]['vy']
    vz = df[df['t'] == t]['vz']
    K = 1/2 * m * (vx**2 + vy**2 + vz**2).sum() / N
    t_array.append(t)
    K_array.append(K)

plt.scatter(t_array, K_array)
plt.show()
