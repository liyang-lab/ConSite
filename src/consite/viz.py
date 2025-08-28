from __future__ import annotations
from pathlib import Path
from typing import Any

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import patheffects as pe

from matplotlib.colors import to_rgba


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
    cons_values: list[float] | np.ndarray | None = None,
    *,
    cons_clip: tuple[float, float] = (5.0, 95.0),   # percentile clip range
    cons_gamma: float = 0.75,                       # <1 brightens mids; >1 darkens
    cons_smooth: int = 0,                           # 0=off, else odd window (e.g., 3)
    cons_show_scale: bool = True,                   # draw a tiny grayscale legend
) -> None:
    """
    ...
    NEW controls:
      - cons_clip: percentile clip (low, high) before normalization (forces contrast)
      - cons_gamma: gamma applied after 0-1 normalization (shapes midtones)
      - cons_smooth: optional moving-average window over the span
      - cons_show_scale: draw a mini grayscale bar with numeric endpoints
    """
    start, end = int(hit.ali_start), int(hit.ali_end)
    if end < start:
        start, end = end, start

    xs = np.arange(start, end + 1)
    subseq = seq[start - 1 : end]

    fig_w = max(10.0, 0.12 * len(xs))
    fig, ax = plt.subplots(figsize=(fig_w, 1.8), constrained_layout=True)
    ax.set_facecolor("white")

    # translucent domain band behind text
    ax.axvspan(start - 0.5, end + 0.5, color="#77b3d5", alpha=0.25, zorder=0)

    # ---------- GRADIENT BACKGROUND (improved) ----------
    if cons_values is not None and len(xs) > 0:
        vec = np.asarray(cons_values, dtype=float)

        # Prefer slicing by absolute positions if vector covers full sequence;
        # otherwise accept a vector already matching the hit span.
        if vec.size >= len(seq):
            span_vals = vec[start - 1 : end]
        elif vec.size == len(xs):
            span_vals = vec
        else:
            span_vals = None  # length mismatch -> skip gracefully

        if span_vals is not None and span_vals.size == len(xs):
            v = np.array(span_vals, dtype=float)
            v[~np.isfinite(v)] = 0.0

            # Percentile clip to force usable contrast
            lo_p, hi_p = cons_clip
            lo = float(np.percentile(v, lo_p)) if hi_p > lo_p else float(np.min(v))
            hi = float(np.percentile(v, hi_p)) if hi_p > lo_p else float(np.max(v))
            if hi <= lo:
                lo, hi = float(np.min(v)), float(np.max(v))
            if hi > lo:
                v = (v - lo) / (hi - lo)
            else:
                v = np.zeros_like(v)

            # Optional smoothing (simple moving average)
            if cons_smooth and cons_smooth > 1 and cons_smooth % 2 == 1:
                k = cons_smooth
                pad = k // 2
                v_pad = np.pad(v, (pad, pad), mode="edge")
                kernel = np.ones(k, dtype=float) / k
                v = np.convolve(v_pad, kernel, mode="valid")

            # Gamma shaping for midtone richness
            v = np.clip(v, 0.0, 1.0) ** cons_gamma

            # grayscale color (dark=high)
            def _gray(val: float) -> tuple[float, float, float]:
                g = 1.0 - float(val)
                return (g, g, g)

            # Paint before letters/outlines, above the blue band
            for x, score in zip(xs, v):
                ax.add_patch(
                    Rectangle(
                        (x - 0.5, -0.35),
                        1.0,
                        0.7,
                        facecolor=_gray(score),
                        edgecolor="none",
                        linewidth=0.0,
                        zorder=0.6,
                    )
                )

            # Optional mini scale bar that "says the values" without overlapping
            if cons_show_scale:
                # inset in axes-fraction coords: x, y, w, h
                ax2 = ax.inset_axes([0.012, 0.80, 0.20, 0.12], zorder=5)
                ax2.set_facecolor("white")  # white background (no alpha, avoids tuple error)
                ax2.set_alpha(0.85)  # optional: slight transparency
                for spine in ("top", "right", "left", "bottom"):
                    ax2.spines[spine].set_visible(True)

                # draw a horizontal grayscale bar leaving room for labels below
                grad = np.linspace(0, 1, 256).reshape(1, -1)
                ax2.imshow(grad, aspect="auto", cmap="gray_r", extent=[0, 1, 0.45, 0.95], clip_on=False)

                # no ticks—avoid layout bleed into the main axes
                ax2.set_xticks([]); ax2.set_yticks([])
                ax2.set_xlim(0, 1); ax2.set_ylim(0, 1)

                # label the clipped endpoints and midpoint inside the inset
                ticks_vals = [lo, (lo + hi) / 2.0, hi]
                ax2.text(0.0, 0.20, f"{ticks_vals[0]:.2f}", ha="left",  va="center", fontsize=7, clip_on=False)
                ax2.text(0.5, 0.20, f"{ticks_vals[1]:.2f}", ha="center", va="center", fontsize=7, clip_on=False)
                ax2.text(1.0, 0.20, f"{ticks_vals[2]:.2f}", ha="right", va="center", fontsize=7, clip_on=False)

                # compact title inside the inset
                ax2.text(0.01, 0.98, "bg=JSD", ha="left", va="top", fontsize=7, clip_on=False)

    # ---------- END GRADIENT BACKGROUND ----------

    # letters
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

    # light cell outlines
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

    # conserved markers (hollow red)
    cons_mask = [p in conserved for p in xs]
    if any(cons_mask):
        xs_cons = xs[cons_mask]
        y_cons = np.full_like(xs_cons, 0.18, dtype=float)
        ms: Any = "o"
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
