import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

INPUT = "out.tsv"
CHUNK_ROWS = 10_000_000

# Box limits
L = 1.1
xmin, xmax = -L, L
ymin, ymax = -L, L

rng = np.random.default_rng(12345)

# Choose red particle ids from the first chunk only
first_chunk = pd.read_csv(INPUT, sep="\t", nrows=CHUNK_ROWS)
all_ids = first_chunk["id"].unique()
red_ids = set(rng.choice(all_ids, size=min(10, len(all_ids)), replace=False))

def make_video(df, chunk_index):
    frames = []

    for t, d in df.groupby("t", sort=True):
        xy = d[["x", "y"]].to_numpy()

        colors = np.full((len(d), 4), [0, 0, 1, 1.0])
        colors[d["id"].isin(red_ids).to_numpy()] = [1, 0, 0, 1.0]

        K = 0.5 * (d["vx"]**2 + d["vy"]**2 + d["vz"]**2).sum() * 1e-3

        frames.append((t, xy, colors, K))

    fig, ax = plt.subplots(figsize=(6, 6), dpi=300)

    scat = ax.scatter([], [], s=1)

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
        info_text.set_text("")
        return scat, info_text

    def update(frame):
        t, xy, colors, K = frames[frame]

        scat.set_offsets(xy)
        scat.set_facecolors(colors)

        info_text.set_text(f"t = {t:.3f} seconds\nK = {K:.1f} joules")

        return scat, info_text

    ani = FuncAnimation(
        fig,
        update,
        frames=len(frames),
        init_func=init,
        interval=33,
        blit=True,
    )

    output = f"gas_chunk_{chunk_index:03d}.mp4"

    ani.save(
        output,
        fps=15,
        dpi=300,
        bitrate=1800,
        extra_args=["-pix_fmt", "yuv420p"],
    )

    plt.close(fig)

    print(f"Saved {output}")


# Process first chunk
make_video(first_chunk, 0)

# Process remaining chunks
reader = pd.read_csv(
    INPUT,
    sep="\t",
    chunksize=CHUNK_ROWS,
    skiprows=range(1, CHUNK_ROWS + 1),
)

for i, chunk in enumerate(reader, start=1):
    make_video(chunk, i)
