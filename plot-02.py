import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv('out.tsv', sep = '\t')
x = df['x']
y = df['y']
z = df['z']

fig = plt.figure()
ax = fig.add_subplot(projection = '3d')

ax.scatter(x, y, z)

plt.show()
