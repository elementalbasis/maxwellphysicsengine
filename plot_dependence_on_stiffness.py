import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.optimize as spo

# Read the file
df = pd.read_csv('data.tsv', sep = '\t')
k = df['k']
tau = df['tau']
u_tau = df['u_tau']

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
s = 5
#fig1 = plt.figure()
#fig2 = plt.figure()
#fig3 = plt.figure()
#fig4 = plt.figure()
fig, axs = plt.subplots(1, 2, constrained_layout = True, figsize = (12,8))
axs1 = axs[0]
axs2 = axs[1]

# Plot the data
x = np.log(k)
y = np.log(tau)
#x = 1/k**2
#y = tau
axs1.scatter(x, y, s = s)
axs1.set_xlabel(r'$\ln k$')
axs1.set_ylabel(r'$\ln \tau$')
#axs1.set_title(r'Logarithmic relaxation time vs. logarithmic stiffness')
#axs1.set_xlabel(r'$\frac{1}{k^2}$')
#axs1.set_ylabel(r'$\tau$')

'''
# Plot the data
plt.scatter(k, tau)
plt.xlabel(r'$k$')
plt.ylabel(r'$\tau$')
plt.title(r'Relaxation time vs. stiffness')
'''

# Perform a linear fit
X = []
Y = []
for i in range(len(x)):
    u = x[i]
    if u < 10.75:
    #if True:
        X.append(x[i])
        Y.append(y[i])

m, b = np.polyfit(X, Y, 1)
textstr = rf'Fit: $m = {m:.5f}$, $b = {b:.5f}$'
x_plot = np.arange(min(x), max(x), 1e-3)
axs1.plot(x_plot, m*x_plot+b, label = textstr, color = 'red')
axs1.legend()


def g(x):
    return -x + np.sqrt(1 + x**2)

# Perform an exponential fit
#def f(k, tau0, A):
    #return tau0 + A / k**2
def f(x, a, b, x0, y0):
    return a*g((x-x0)/b) + y0
popt, pcov = spo.curve_fit(f, x, y,
                           #p0 = [tau.iloc[-1], tau.iloc[0] - tau.iloc[-1], 1.0],
                           #p0 = [0.5, 20],
                           p0 = [1, 1, 11, -1],
                           #bounds = ([0, 0, 0],
                           #          [2*tau.iloc[-1], 2*tau.iloc[0], np.inf]
                           #    )
                           )
#tau0, A = popt
#a, b, k0, tau0 = popt
a, b, x0, y0 = popt
u_a = np.sqrt(pcov[0][0])
u_b = np.sqrt(pcov[1][1])
u_x0 = np.sqrt(pcov[2][2])
u_y0 = np.sqrt(pcov[3][3])
print('a =', a, '+/-', u_a)
print('b =', b, '+/-', u_b)
print('x0 =', x0, '+/-', u_x0)
print('y0 =', y0, '+/-', u_y0)
#u_tau0 = np.sqrt(pcov[0][0])
#u_A = np.sqrt(pcov[1][1])
#k_plot = np.arange(min(k), max(k), 1)
#tau_fit = f(k_plot, tau0, A)
#tau_fit = f(k_plot, a, b, k0, tau0)
x_plot = np.arange(min(x), max(x), 1e-4)
y_plot = f(x_plot, a, b, x0, y0)
axs2.scatter(x, y, s = s)
axs2.plot(x_plot, y_plot, color = 'red')
#axs2.set_yscale('log')
#axs2.set_xscale('log')
#print('%s\t%f\t%f\t%f\t%f\t%f\t%f' % (filename, A_eq, A0, tau, u_A_eq, u_A0, u_tau))
#print('%s\t%f\t%f' % (filename, tau, u_tau))

# Show the graph
plt.show()
