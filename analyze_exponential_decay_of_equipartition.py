import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import scipy.optimize as spo

# Read the file.
filename = sys.argv[1]
#print(filename)
df = pd.read_csv(filename, sep = '\t')

# Only use data between 5% and 90% of the time range provided, to avoid unstable
# conditions at the beginning at the end of the simulation.
proportion_start = 0.05
proportion_end = 0.90
start = int(len(df) * proportion_start)
end = int(len(df) * proportion_end)
df = df.iloc[start:end]
t = df['t']

# Set the font to the one used by LaTeX.
#plt.rcParams['text.usetex'] = True
#plt.rcParams['font.family'] = 'serif'
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 14,
    "axes.labelsize": 16,
    "axes.titlesize": 16,
    "legend.fontsize": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
})

# Set the dot size for plotting; create the plotting figures.
s = 0.1
#fig1 = plt.figure()
#fig2 = plt.figure()
#fig3 = plt.figure()
#fig4 = plt.figure()
fig, axs = plt.subplots(2, 2, constrained_layout = True, figsize = (12,8))
axs1 = axs[0,0]
axs2 = axs[0,1]
axs3 = axs[1,0]
axs4 = axs[1,1]

# Calculate the average kinetic energy on each degree of freedom.
y0 = (df['Ex'] + df['Ey'] + df['Ez']) / 3

# Calculate the energy anisotropy, which is the standard deviation of the
# energies on the three degrees of freedom.
A = np.sqrt(((df['Ex'] - y0)**2 + (df['Ey'] - y0)**2 + (df['Ez'] - y0)**2)/3)
U = np.log(A)

# Plot the logarithmic energy anisotropy vs. time.
#plt.figure(fig1.number)
axs1.scatter(t, U, s = s)
axs1.set_xlabel(r'$t$ [s]')
axs1.set_ylabel(r'$\log \sigma_E$')
axs1.set_title(r'Logarithmic energy anisotropy vs. time')

# Perform a linear fit on the logarithmic energy anisotropy.
#m, b = np.polyfit(t, U, 1)
#textstr = rf'Fit: $m = {m:.5f}$, $b = {b:.5f}$'
#axs1.plot(t, m*t + b, color = 'red', label = textstr)
#axs1.legend()

# Plot the energy anisotropy vs. time.
#plt.figure(fig2.number)
axs2.scatter(t, A, s = s)
axs2.set_xlabel(r'$t$ [s]')
axs2.set_ylabel(r'$\sigma_E$ [J]')
axs2.set_title('Energy anisotropy vs. time')

# Perform an exponential fit to the energy anisotropy
def f(t, A_eq, A0, tau):
    return A_eq + A0 * np.exp(-t / tau)
popt, pcov = spo.curve_fit(f, t, A,
                           p0 = [A.iloc[-1], A.iloc[0] - A.iloc[-1], 1.0],
                           bounds = ([0, 0, 0],
                                     [2*A.iloc[-1], 2*A.iloc[0], np.inf]
                               )
                           )
A_eq, A0, tau = popt
u_A_eq = np.sqrt(pcov[0][0])
u_A0 = np.sqrt(pcov[1][1])
u_tau = np.sqrt(pcov[2][2])
A_fit = f(t, A_eq, A0, tau)
axs2.plot(t, A_fit, color = 'red')
#print('%s\t%f\t%f\t%f\t%f\t%f\t%f' % (filename, A_eq, A0, tau, u_A_eq, u_A0, u_tau))
print('%s\t%f\t%f' % (filename, tau, u_tau))

# Plot the average directional kinetic energy vs. time.
#plt.figure(fig3.number)
axs3.scatter(t, y0, s = s)
axs3.set_xlabel(r'$t$ [s]')
axs3.set_ylabel(r'$\overline{E}$ [J]')
axs3.set_title('Average directional kinetic energy vs. time')

# Plot the inidivual energies.
#plt.figure(fig4.number)
axs4.scatter(t, df['Ex'], s = s, label = '$E_x$')
axs4.scatter(t, df['Ey'], s = s, label = '$E_y$')
axs4.scatter(t, df['Ez'], s = s, label = '$E_z$')
axs4.set_xlabel(r'$t$ [s]')
axs4.set_ylabel(r'$E_i$ [J]')
axs4.set_title('Directional kinetic energy vs. time')
axs4.legend(markerscale = 20) # Enlarge the dots in the legend

# Show all graphs for a specified time.
#plt.show()
#plt.show(block = False)
#plt.pause(3)
#plt.close()

# Save the figure.
newfilename = filename.replace('.tsv', '.png').replace('run', 'figure')
fig.savefig(newfilename)

#for s in ['Ex', 'Ey', 'Ez']:
    #y = df[s]
    #u = abs(df[s] - y0)
    #plt.scatter(t, u)

    #start = 0
    #proportion = 0.95
    #finish = int(np.round(proportion * len(y)))

    #T_full = t
    #T = T_full.iloc[start:finish]
    #Y_full = np.log(abs(y - y0))
    #Y = Y_full.iloc[start:finish]

    #label = r'$\log|%s - E_0|$' % {'Ex': 'E_x', 'Ey': 'E_y', 'Ez': 'E_z'}[s]
    #label = r'$%s$' % {'Ex': 'E_x', 'Ey': 'E_y', 'Ez': 'E_z'}[s]
    #plt.scatter(T, Y, s = 3, label = label)
    #plt.scatter(t, y, s = 3, label = label)

    #a, b = np.polyfit(T, Y, 1)
    #a, b = np.polyfit(T, y.iloc[start:finish], 1)
    #tau = 1/a
    #print("%s:\ta = %f\tb = %f\ttau = %f" % (s, a, b, tau))
    #print("%s:\ta = %f\tb = %f" % (s, a, b))
    #plt.plot(T, a * T + b, color = 'red')
    #plt.plot(t, [y0] * len(t), color = 'red', linestyle = 'dashed')
    #plt.scatter(t, u)
    #plt.plot(t, np.exp(a * t + b), color = 'red')

#plt.yscale('log')
#plt.xlabel('$t$ (s)')
#plt.ylabel(r'$\log|E - E_0|$')
#plt.ylabel(r'$E$ (J)')
#plt.title(r'$\log|E - E_0|$ vs. $t$ for each degree of freedom')
#plt.title(r'$E$ vs. $t$ for each degree of freedom')
#plt.legend()
#plt.show()
