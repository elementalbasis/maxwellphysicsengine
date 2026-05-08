import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

INPUT = "out.tsv"
CHUNK_ROWS = 7_200_000

# Box limits
L = 1.1
xmin, xmax = -L, L
ymin, ymax = -L, L

# Particle rendering settings.
# Keep all particles the same size; only alpha changes with z.
PARTICLE_SIZE = 1.0
MIN_PARTICLE_ALPHA = 0.10
MAX_PARTICLE_ALPHA = 1.00

# Red-particle trail settings.
TRAIL_FRAMES = 25
TRAIL_SIZE = PARTICLE_SIZE

rng = np.random.default_rng(12345)


def z_depth_factor(z):
    """
    Map z from roughly [-L, L] to [0, 1].
    Negative z becomes transparent; positive z becomes solid.
    """
    return np.clip((z + L) / (2 * L), 0.0, 1.0)


def alphas_from_z(z):
    depth = z_depth_factor(z)
    return MIN_PARTICLE_ALPHA + depth * (MAX_PARTICLE_ALPHA - MIN_PARTICLE_ALPHA)


# Choose red particle ids from the first chunk only.
first_chunk = pd.read_csv(INPUT, sep="\t", nrows=CHUNK_ROWS)
all_ids = first_chunk["id"].unique()
red_ids = set(rng.choice(all_ids, size=min(10, len(all_ids)), replace=False))


def make_video(df, chunk_index, carryover_red_history=None):
    """
    Render one chunk.

    carryover_red_history is a small list of red-particle snapshots from the
    previous chunk. Each item is (red_xy, red_ids_this_frame). It lets the first
    frames of this chunk draw a trail that began in the previous chunk, without
    keeping the whole previous chunk in memory.
    """
    if carryover_red_history is None:
        carryover_red_history = []

    frames = []
    current_chunk_red_history = []

    for t, d in df.groupby("t", sort=True):
        xy = d[["x", "y"]].to_numpy()
        z = d["z"].to_numpy()
        ids = d["id"].to_numpy()
        red_mask = d["id"].isin(red_ids).to_numpy()

        alphas = alphas_from_z(z)

        colors = np.zeros((len(d), 4))
        colors[:, 2] = 1.0          # blue particles
        colors[:, 3] = alphas       # z-depth transparency
        colors[red_mask, 0] = 1.0   # red particles
        colors[red_mask, 2] = 0.0

        red_xy = xy[red_mask]
        red_ids_this_frame = ids[red_mask]
        red_alphas = alphas[red_mask]

        # Current-frame alpha lookup for red particles.
        # Trails use the current z-depth alpha of the corresponding live red dot,
        # then fade toward 0 as the trail point gets older.
        current_red_alpha_by_id = dict(zip(red_ids_this_frame, red_alphas))

        Kx = 0.5 * (d["vx"] ** 2).sum() * 1e-3
        Ky = 0.5 * (d["vy"] ** 2).sum() * 1e-3
        Kz = 0.5 * (d["vz"] ** 2).sum() * 1e-3
        K = Kx + Ky + Kz

        frames.append((
            t,
            xy,
            colors,
            red_xy,
            red_ids_this_frame,
            red_alphas,
            current_red_alpha_by_id,
            K,
            Kx,
            Ky,
            Kz,
        ))

        # Store only the tiny red-particle snapshot needed for trails.
        current_chunk_red_history.append((red_xy, red_ids_this_frame))

    fig, ax = plt.subplots(figsize=(6, 6), dpi=300)

    # Draw trail first so live particles sit on top of it.
    trail_scat = ax.scatter([], [], s=[])
    scat = ax.scatter([], [], s=PARTICLE_SIZE)

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
        empty_offsets = np.empty((0, 2))
        scat.set_offsets(empty_offsets)
        scat.set_sizes([])
        scat.set_facecolors(np.empty((0, 4)))
        trail_scat.set_offsets(empty_offsets)
        trail_scat.set_sizes([])
        trail_scat.set_facecolors(np.empty((0, 4)))
        info_text.set_text("")
        return trail_scat, scat, info_text

    def update(frame):
        (
            t,
            xy,
            colors,
            red_xy,
            red_ids_this_frame,
            red_alphas,
            current_red_alpha_by_id,
            K,
            Kx,
            Ky,
            Kz,
        ) = frames[frame]

        trail_offsets = []
        trail_sizes = []
        trail_colors = []

        # Include red-particle snapshots from the previous chunk, then the
        # snapshots before this frame in the current chunk. This is the key fix:
        # frame 0 of chunk N can still draw the ending trail from chunk N-1.
        prior_red_history = carryover_red_history + current_chunk_red_history[:frame]
        prior_red_history = prior_red_history[-TRAIL_FRAMES:]

        for history_index, (old_red_xy, old_red_ids) in enumerate(reversed(prior_red_history), start=1):
            age = history_index

            # age = 1 means right behind the dot: alpha is almost the current dot alpha.
            # age = TRAIL_FRAMES means old tail: alpha is close to 0.
            age_fade = (TRAIL_FRAMES - age + 1) / TRAIL_FRAMES

            if len(old_red_xy) == 0:
                continue

            # Use current z-alpha for each corresponding red particle, not the old z-alpha.
            base_alpha = np.array([
                current_red_alpha_by_id.get(pid, 0.0)
                for pid in old_red_ids
            ])
            trail_alpha = base_alpha * age_fade

            keep = trail_alpha > 0.0
            if not np.any(keep):
                continue

            trail_offsets.append(old_red_xy[keep])
            trail_sizes.append(np.full(np.count_nonzero(keep), TRAIL_SIZE))

            old_colors = np.zeros((np.count_nonzero(keep), 4))
            old_colors[:, 0] = 1.0
            old_colors[:, 3] = trail_alpha[keep]
            trail_colors.append(old_colors)

        if trail_offsets:
            trail_scat.set_offsets(np.vstack(trail_offsets))
            trail_scat.set_sizes(np.concatenate(trail_sizes))
            trail_scat.set_facecolors(np.vstack(trail_colors))
        else:
            trail_scat.set_offsets(np.empty((0, 2)))
            trail_scat.set_sizes([])
            trail_scat.set_facecolors(np.empty((0, 4)))

        scat.set_offsets(xy)
        scat.set_sizes(np.full(len(xy), PARTICLE_SIZE))
        scat.set_facecolors(colors)

        info_text.set_text(
            f"t = {t:.3f} s\n"
            f"K = {K/1000:.1f} kJ\n"
            f"Kx = {Kx/1000:.1f} kJ\n"
            f"Ky = {Ky/1000:.1f} kJ\n"
            f"Kz = {Kz/1000:.1f} kJ"
        )

        return trail_scat, scat, info_text

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
        fps=30,
        dpi=300,
        # bitrate=20000,
        extra_args=["-crf", "18", "-pix_fmt", "yuv420p", "-preset", "slow"],
    )

    plt.close(fig)

    print(f"Saved {output}")

    # Return only the last TRAIL_FRAMES red snapshots for the next chunk.
    return current_chunk_red_history[-TRAIL_FRAMES:]


# Process first chunk.
red_history = make_video(first_chunk, 0)

# Process remaining chunks.
reader = pd.read_csv(
    INPUT,
    sep="\t",
    chunksize=CHUNK_ROWS,
    skiprows=range(1, CHUNK_ROWS + 1),
)

for i, chunk in enumerate(reader, start=1):
    red_history = make_video(chunk, i, carryover_red_history=red_history)
