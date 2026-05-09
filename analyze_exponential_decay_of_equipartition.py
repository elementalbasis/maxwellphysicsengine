import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

#filename = "/home/WindowsShared/maxwell/run-2026-05-08_i=001.tsv"
filename = sys.argv[1]
df = pd.read_csv(filename, sep = '\t')
t = df['t']

print(filename)

for s in {'Ex', 'Ey', 'Ez'}:
    y0 = 5000
    y = df[s]
    u = abs(df[s] - y0)
    #plt.scatter(t, u)

    start = 0
    proportion = 0.5
    finish = int(np.round(proportion * len(y)))

    T = t.iloc[start:finish]
    Y = np.log(abs(y.iloc[start:finish] - y0))
    plt.scatter(T, Y)

    a, b = np.polyfit(T, Y, 1)
    print("%s:\ta = %f\tb = %f" % (s, a, b))
    plt.plot(T, a * T + b)

#plt.yscale('log')
plt.show()
