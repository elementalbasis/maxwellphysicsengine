import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

df = pd.read_csv("moro.tsv", sep="\t", nrows=30000000)

ids = df["id"].drop_duplicates().to_numpy()

if len(ids) != 2:
    raise ValueError(f"Expected exactly 2 particles, but found {len(ids)} UUIDs.")

id1, id2 = ids

df = df[df["id"].isin([id1, id2])].copy()

L = 2
xmin, xmax = -L/2, L/2
ymin, ymax = -L/2, L/2

frames = []

for t, d in df.groupby("t", sort=True):
    xy = d[["x", "y"]].to_numpy()

    colors = np.where(
        d["id"].to_numpy()[:, None] == id1,
        np.array([1, 0, 0, 1.0]),
        np.array([0, 0, 1, 1.0])
    )

    K = 0.5 * (d["vx"]**2 + d["vy"]**2 + d["vz"]**2).sum() * 1e-3

    v1 = float(d.loc[d["id"] == id1, "vx"].iloc[0])
    v2 = float(d.loc[d["id"] == id2, "vx"].iloc[0])

    frames.append((t, xy, colors, K, v1, v2))

fig, ax = plt.subplots(figsize=(6, 6), dpi=150)

scat = ax.scatter([], [], s=20)

info_text = ax.text(
    0.03, 0.97, "",
    transform=ax.transAxes,
    ha="left",
    va="top",
    fontsize=10,
)

ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_aspect("equal", adjustable="box")

ax.set_xticks([])
ax.set_yticks([])
ax.set_xlabel("")
ax.set_ylabel("")
ax.set_title("")

ax.set_position([0.08, 0.08, 0.84, 0.84])

def init():
    scat.set_offsets(np.empty((0, 2)))
    scat.set_facecolors([])
    info_text.set_text("")
    return scat, info_text

def update(frame_index):
    t, xy, colors, K, v1, v2 = frames[frame_index]

    scat.set_offsets(xy)
    scat.set_facecolors(colors)

    info_text.set_text(
        f"t = {t:.3f} seconds\n"
        f"K = {K:.6f} joules\n"
        f"particle 1 vx = {v1:.6f}\n"
        f"particle 2 vx = {v2:.6f}"
    )

    return scat, info_text

ani = FuncAnimation(
    fig,
    update,
    frames=len(frames),
    init_func=init,
    interval=33,
    blit=True,
)

ani.save(
    "collision.mp4",
    fps=60,
    dpi=150,
    bitrate=1800,
    extra_args=["-pix_fmt", "yuv420p"],
)
