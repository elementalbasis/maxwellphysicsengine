# Minimal 3D particle animation (TSV -> MP4), no IDs, no sorting.
# TSV must have columns: t, x, y, z  (multiple rows per time t)

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import animation

# ---- user knobs ----
INPUT_PATH  = "data.tsv"       # tab-separated values with columns: t,x,y,z
OUTPUT_PATH = "particles.mp4"  # MP4 output (requires ffmpeg in PATH)
FPS         = 30
POINT_SIZE  = 10
FIG_SIZE    = 6
BG          = "#111111"
FG          = "#EEEEEE"
DPI         = 120
# --------------------

def load_tsv(path):
    # NOTE: tab is "\t" (backslash-t), not "/t".
    df = pd.read_csv(path, sep="\t")
    if not {"t","x","y","z"}.issubset(df.columns):
        raise ValueError("TSV must have columns: t, x, y, z")
    # No sorting; assuming input already sorted by t
    times = df["t"].to_numpy()
    unique_times = np.unique(times)
    frames = [df[times == t][["x","y","z"]].to_numpy() for t in unique_times]
    return frames, unique_times

frames, times = load_tsv(INPUT_PATH)

# Bounds with cubic aspect so axes have equal scale
all_xyz = np.vstack(frames)                 # shape (sum_i N_i, 3)
mins = np.nanmin(all_xyz, axis=0)
maxs = np.nanmax(all_xyz, axis=0)
span = float(np.max(maxs - mins)) or 1.0
center = (mins + maxs) / 2.0
mins = center - span/2.0
maxs = center + span/2.0

# Figure + 3D axes
fig = plt.figure(figsize=(FIG_SIZE, FIG_SIZE), facecolor=BG)
ax = fig.add_subplot(111, projection="3d", facecolor=BG)
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

ax.set_xlim(mins[0], maxs[0])
ax.set_ylim(mins[1], maxs[1])
ax.set_zlim(mins[2], maxs[2])
ax.tick_params(colors=FG, which="both")
try:
    ax.xaxis.pane.set_edgecolor(FG); ax.yaxis.pane.set_edgecolor(FG); ax.zaxis.pane.set_edgecolor(FG)
    ax.xaxis.pane.set_alpha(0.1);    ax.yaxis.pane.set_alpha(0.1);    ax.zaxis.pane.set_alpha(0.1)
except Exception:
    pass

# First frame
first = frames[0]
scat = ax.scatter(first[:,0], first[:,1], first[:,2], s=POINT_SIZE, depthshade=False)
title = ax.text2D(0.02, 0.95, f"t = {times[0]}", transform=ax.transAxes, color=FG)

def update(i):
    xyz = frames[i]
    scat._offsets3d = (xyz[:,0], xyz[:,1], xyz[:,2])
    title.set_text(f"t = {times[i]}")
    return scat, title

anim = animation.FuncAnimation(fig, update, frames=len(frames), interval=1000/FPS, blit=False)

# Save to MP4 (requires ffmpeg)
writer = animation.FFMpegWriter(fps=FPS, metadata={"artist": "particle_anim"})
anim.save(OUTPUT_PATH, writer=writer, dpi=DPI, savefig_kwargs={"facecolor": BG})
print(f"Saved animation to {OUTPUT_PATH}")

