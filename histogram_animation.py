import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

df = pd.read_csv("out.tsv", sep="\t")
df = df.drop_duplicates()

# Settings
bins = 40
fps = 30

# Precompute time steps
times = np.sort(df["t"].unique())

def make_animation(component):
    # Build frames
    frames = []
    for i, t in enumerate(times):
        d = df[df["t"] == t]
        v = d[component].to_numpy()
        frames.append((i, t, v))

    # Global histogram range
    v_all = df[component].to_numpy()
    vmin, vmax = v_all.min(), v_all.max()
    vmax_abs = max(abs(vmin), abs(vmax))
    hist_range = (-vmax_abs, vmax_abs)

    # Fixed y-axis
    max_count = 0
    for _, _, v in frames:
        counts, _ = np.histogram(v, bins=bins, range=hist_range)
        max_count = max(max_count, counts.max())

    fig, ax = plt.subplots(figsize=(6, 4), dpi=150)

    def init():
        ax.clear()
        ax.set_xlim(*hist_range)
        ax.set_ylim(0, max_count * 1.1)
        ax.set_xlabel(component)
        ax.set_ylabel("Count")
        return []

    def update(frame):
        step, t, v = frames[frame]

        ax.clear()

        counts, edges = np.histogram(v, bins=bins, range=hist_range)

        ax.bar(
            edges[:-1],
            counts,
            width=np.diff(edges),
            align="edge",
        )

        ax.set_xlim(*hist_range)
        ax.set_ylim(0, max_count * 1.1)
        ax.set_xlabel(component)
        ax.set_ylabel("Count")
        ax.set_title(f"Histogram of {component}, step = {step}, t = {t:.3f}")

        return []

    ani = FuncAnimation(
        fig,
        update,
        frames=len(frames),
        init_func=init,
        interval=1000 / fps,
        blit=False,
    )

    ani.save(
        f"{component}_histogram.mp4",
        fps=fps,
        dpi=150,
        bitrate=2500,
        extra_args=["-pix_fmt", "yuv420p"],
    )

    plt.close(fig)  # important: avoid memory buildup


# Generate all three
for comp in ["vx", "vy", "vz"]:
    make_animation(comp)
