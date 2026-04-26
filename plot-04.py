import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

df = pd.read_csv("out.tsv", sep = '\t')

# Get box limits from data, or set manually
L = 1.1
xmin = -L
xmax = L
ymin = -L
ymax = L

times = sorted(df["t"].unique())

fig, ax = plt.subplots()
scat = ax.scatter([], [], s=1)

ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_aspect("equal")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("Gas particles in a box")

def init():
    scat.set_offsets(np.empty((0, 2)))
    return scat,

def update(frame):
    t = times[frame]
    d = df[df["t"] == t]

    xy = d[["x", "y"]].to_numpy()
    scat.set_offsets(xy)

    ax.set_title(f"Gas particles in a box, t = {t:.3f}")
    return scat,

ani = FuncAnimation(
    fig,
    update,
    frames=len(times),
    init_func=init,
    interval=60,
    blit=True,
)

plt.show()
