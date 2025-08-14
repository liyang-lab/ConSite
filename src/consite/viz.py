from __future__ import annotations
from pathlib import Path
from typing import Any

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import patheffects as pe

from .utils import Hit


def plot_domain_map(seq_len: int, hits: list[Hit], conserved_idx: list[int], out_png: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 1.8))
    ax.plot([1, seq_len], [0.5, 0.5], color="#444", linewidth=1.2, zorder=1)

    for h in hits:
        x0 = min(int(h.ali_start), int(h.ali_end))
        width = abs(int(h.ali_end) - int(h.ali_start)) + 1
        ax.add_patch(
            Rectangle(
                (x0, 0.25),
                width,
                0.5,
                facecolor="#7fb3d5",
                alpha=0.8,
                linewidth=0,
                zorder=0,
                label=h.family,
            )
        )

    if conserved_idx:
        ax.scatter(
            conserved_idx,
            [0.5] * len(conserved_idx),
            s=30,
            facecolors="none",
            edgecolors="#d62728",
            linewidths=1.5,
            zorder=2,
        )

    ax.set_xlim(1, seq_len)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xlabel("Sequence position")
    ax.set_title("Domain map with conserved sites")
    fig.tight_layout()
    fig.savefig(str(out_png), dpi=200)
    plt.close(fig)


def plot_conservation_track(scores: dict, out_png: Path, title: str = "Conservation (JSD)") -> None:
    jsd = scores["jsd"]
    fig, ax = plt.subplots(figsize=(10, 2.2))
    ax.plot(np.arange(1, len(jsd) + 1), jsd, linewidth=1.5)
    ax.set_xlim(1, len(jsd))
    ax.set_ylim(0, 1)
    ax.set_xlabel("Sequence position (aligned)")
    ax.set_ylabel("JSD (norm)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(str(out_png), dpi=200)
    plt.close(fig)


def plot_alignment_panel(
    seq: str,
    hit: Hit,
    conserved: set[int],    
    out_png: Path,
) -> None:
    """
    Per-hit panel that prioritizes readability:
      • Letters are white with a thin black stroke (always readable).
      • Light cell outlines behind letters (no solid boxes covering text).
      • Conserved sites = hollow red circles slightly above the baseline.
    """
    start, end = int(hit.ali_start), int(hit.ali_end)
    if end < start:
        start, end = end, start

    xs = np.arange(start, end + 1)
    subseq = seq[start - 1 : end]

    # width grows with length so characters don’t cramp
    fig_w = max(10.0, 0.12 * len(xs))
    fig, ax = plt.subplots(figsize=(fig_w, 1.8), constrained_layout=True)
    ax.set_facecolor("white")

    # translucent domain band behind text
    ax.axvspan(start - 0.5, end + 0.5, color="#77b3d5", alpha=0.25, zorder=0)

    # letters: white fill + black stroke so they pop on any background
    text_effect = [pe.withStroke(linewidth=1.0, foreground="black", alpha=0.6)]
    for x, aa in zip(xs, subseq):
        ax.text(
            x,
            0.0,
            aa,
            ha="center",
            va="center",
            fontsize=12,
            family="DejaVu Sans Mono",
            color="white",
            path_effects=text_effect,
            zorder=3,
        )

    # very light cell outlines to help the eye track positions
    for x in xs:
        ax.add_patch(
            Rectangle(
                (x - 0.5, -0.35),
                1.0,
                0.7,
                facecolor="none",
                edgecolor="#cccccc",
                linewidth=0.6,
                zorder=1,
            )
        )

    # conserved markers: hollow red circles slightly above the text baseline
    cons_mask = [p in conserved for p in xs]
    if any(cons_mask):
        xs_cons = xs[cons_mask]
        y_cons = np.full_like(xs_cons, 0.18, dtype=float)
        ms: Any = "o"  # appease type checker; matplotlib accepts string markers
        ax.scatter(
            xs_cons,
            y_cons,
            s=36,
            marker=ms,
            facecolors="none",
            edgecolors="#d62728",
            linewidths=1.5,
            zorder=4,
        )

    # cosmetics
    ax.set_xlim(start - 0.5, end + 0.5)
    ax.set_ylim(-0.6, 0.6)
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ("top", "right", "left", "bottom"):
        ax.spines[spine].set_visible(False)
    ax.set_title(f"{hit.family}  {start}-{end}", fontsize=13, pad=6)

    fig.savefig(str(out_png), dpi=200)
    plt.close(fig)
