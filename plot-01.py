import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('test.tsv', sep = '\t')
t = df['t']
x = df['x']
v = df['v']

plt.plot(t, x)
plt.show()
