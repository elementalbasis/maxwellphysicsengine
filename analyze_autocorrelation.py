import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.signal import correlate

# Read the file.
filename = sys.argv[1]
df = pd.read_csv(filename, sep = '\t')
t = df['t']
Ex = df['Ex']
Ey = df['Ey']
Ez = df['Ez']

# Set the dot size for plotting.
s = 0.1

# Pick a time index where the gas has already been thermalized.
start = 150000
finish = len(t)
dt = t.iloc[1] - t.iloc[0]
#plt.scatter(t.iloc[start:finish], Ez.iloc[start:finish], s = s)
#plt.show()

t = t.iloc[start:finish]
t = t.to_numpy()

#fig1 = plt.figure()
#fig2 = plt.figure()
fig, axs = plt.subplots(2, 1, constrained_layout = True, figsize = (6, 8))
axs1 = axs[0]
axs2 = axs[1]

for E in [Ex, Ey, Ez]:
    x = E - np.mean(E)
    x = x.iloc[start:finish]
    x = x.to_numpy()

    # Autocorrelation is defined as mean(E[i] * E[i+delta_i]). Therefore, the
    # `correlate()` function returns a (2*N-1)-dimensional array with the
    # autocorrelation at every possible `delta_i`. Here, N is the length of E.
    #
    # Given that this is the case, we take the range of `acf.size//2 : end`
    # because that corresponds to `delta_i = 0 : N`.

    acf = correlate(x, x, mode = 'full')
    acf = acf[acf.size // 2:]
    acf /= acf[0] # Normalize the autocorrelation function.
    #y = np.log(acf)

    # We want to plot the autocorrelation as a function of time lag. Thus, we
    # convert from integer indices into a meaningful time coordinate.
    lags = np.arange(len(acf))
    tau = lags * dt
    #plt.figure(fig1.number)
    axs1.scatter(tau, acf, s = s)
    axs1.set_xlabel(r'$\tau$ [s]')
    axs1.set_ylabel(r'autocorrelation')
    axs1.set_title(r'Autocorrelation vs. time lag')

    # Separately, plot E vs t to see what data we're working with.
    #plt.figure(fig2.number)
    axs2.scatter(t, E.iloc[start:finish], s = s)
    axs2.set_xlabel(r'$t$ [s]')
    axs2.set_ylabel(r'$E$ [J]')
    axs2.set_title(r'Average kinetic energy vs. time')
plt.show()
#newfilename = filename.replace('run', 'figure').replace('tsv', 'png')
#plt.savefig(newfilename)
