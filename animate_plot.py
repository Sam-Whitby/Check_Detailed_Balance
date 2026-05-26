#!/usr/bin/env python3
"""
animate_plot.py  --  Animated matplotlib visualisation for DetailedBalanceChecker

Called by animate.wls:
    python3 animate_plot.py <data.json>

Displays three panels:
  Left:   imshow of the lattice at each recorded step.
  Centre: system energy vs step (line grows as animation progresses).
  Right:  model parameters panel:
            - colour matrix of Jpair coupling strengths (when present)
            - scalar parameters (β, lambdaJ, fieldAmp, …)
            - coupling / field formula strings (when provided by the .wl file)
"""

import json
import os
import re
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import ListedColormap, TwoSlopeNorm
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.patches import Patch


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_data(path: str) -> dict:
    with open(path) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Colour maps for the lattice panel
# ---------------------------------------------------------------------------

_TAB10 = [
    "#4c72b0", "#dd8452", "#55a868", "#c44e52", "#8172b3",
    "#937860", "#da8bc3", "#8c8c8c", "#ccb974", "#64b5cd",
]


def make_colormap(n_types: int) -> ListedColormap:
    colors = ["#1a1a1a"] + [_TAB10[i % len(_TAB10)] for i in range(n_types)]
    return ListedColormap(colors)


def make_simple_colormap() -> ListedColormap:
    return ListedColormap(["#dddddd", "#4c72b0"])


# ---------------------------------------------------------------------------
# Parameter parsing
# ---------------------------------------------------------------------------

def _parse_jpair_key(key: str, n_types: int):
    """Return (a, b) (1-indexed) if key is a valid Jpair<a><b> entry, else None."""
    if not key.startswith("Jpair"):
        return None
    suffix = key[5:]
    # Try all splits of the suffix digits into two integers 1..n_types.
    # Reject any split where either part has a leading zero, which would
    # cause ambiguous matches like "1","011" for Jpair1011 (should be 10,11).
    for i in range(1, len(suffix)):
        a_str, b_str = suffix[:i], suffix[i:]
        if (a_str.isdigit() and b_str.isdigit()
                and not a_str.startswith("0") and not b_str.startswith("0")):
            a, b = int(a_str), int(b_str)
            if 1 <= a <= n_types and 1 <= b <= n_types:
                return a, b
    return None


def build_jpair_matrix(params: dict, n_types: int):
    """
    Extract Jpair<a><b> values from params into a symmetric (n_types × n_types)
    numpy array.  Returns (matrix, found) where found=True if any key matched.
    """
    matrix = np.zeros((n_types, n_types))
    found = False
    for key, val in params.items():
        ab = _parse_jpair_key(key, n_types)
        if ab is not None:
            a, b = ab
            matrix[a - 1, b - 1] = val
            matrix[b - 1, a - 1] = val
            found = True
    return matrix, found


def scalar_params(params: dict, n_types: int) -> dict:
    """Return params that are not Jpair coupling entries."""
    out = {}
    for key, val in params.items():
        if _parse_jpair_key(key, n_types) is None:
            out[key] = val
    return out


# ---------------------------------------------------------------------------
# Right-panel drawing helpers
# ---------------------------------------------------------------------------

def draw_coupling_matrix(ax, cbar_ax, matrix: np.ndarray, n_types: int):
    """
    Draw a colour-coded Jpair(a,b) heatmap on *ax* with a colourbar on
    *cbar_ax*.  Diverging RdBu_r colourmap centred at zero.
    """
    # RdBu: blue = positive J (attractive), white = 0, red = negative J (repulsive)
    vabs = max(float(np.abs(matrix).max()), 1e-10)
    norm = TwoSlopeNorm(vmin=-vabs, vcenter=0.0, vmax=vabs)
    im = ax.imshow(matrix, cmap="RdBu", norm=norm, aspect="equal",
                   interpolation="nearest")

    # Show only the first (1) and last (N) tick labels to avoid clutter
    ax.set_xticks([0, n_types - 1])
    ax.set_xticklabels(["1", str(n_types)], fontsize=9)
    ax.set_yticks([0, n_types - 1])
    ax.set_yticklabels(["1", str(n_types)], fontsize=9)
    ax.set_xlabel("Type b", fontsize=8)
    ax.set_ylabel("Type a", fontsize=8)
    ax.set_title("Jpair(a,b)  blue=+  white=0  red=\u2212", fontsize=8, pad=4)

    cbar = plt.colorbar(im, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=7)

    return im


def draw_info_panel(ax, scalars: dict, coupling_formula: str, field_formula: str):
    """
    Draw scalar parameters and optional formula strings as formatted text.
    """
    ax.axis("off")

    lines = []
    for k, v in scalars.items():
        label = "\u03b2" if k == "beta" else k   # β for beta
        lines.append(f"{label} = {v:.5g}")

    if coupling_formula:
        lines.append("")
        lines.append("J(a, b, d\u00b2) =")          # superscript 2
        lines.append(f"  {coupling_formula}")

    if field_formula:
        lines.append("")
        lines.append("f(x, y, L) =")
        lines.append(f"  {field_formula}")

    if not lines:
        lines = ["(no parameters)"]

    ax.text(0.04, 0.97, "\n".join(lines),
            transform=ax.transAxes,
            fontsize=7.5, va="top", ha="left",
            family="monospace",
            bbox=dict(facecolor="#f0f0f0", edgecolor="#cccccc",
                      boxstyle="round,pad=0.4"))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 animate_plot.py <data.json> [--no-params]")
        sys.exit(1)

    no_params_flag = "--no-params" in sys.argv[2:]
    d = load_data(sys.argv[1])

    steps            = d["steps"]
    states           = d["states"]
    energies         = d["energies"]
    n_rows           = d["grid_rows"]
    n_cols           = d["grid_cols"]
    n_types          = max(int(d.get("n_types", 1)), 1)
    fps              = float(d.get("fps", 10.0))
    delay_ms         = 1000.0 / fps
    simple           = bool(d.get("simple", False))
    algo_name        = d.get("algo_file", "")
    n_frames         = len(steps)
    params           = d.get("params", {})
    coupling_formula = d.get("coupling_formula", "")
    field_formula    = d.get("field_formula", "")

    # Parse Jpair coupling matrix and remaining scalar params
    jpair_matrix, has_jpair = build_jpair_matrix(params, n_types)
    scalars = scalar_params(params, n_types)
    has_params = (bool(params) or bool(coupling_formula) or bool(field_formula)) and not no_params_flag

    # Lattice colour map
    if simple:
        cmap = make_simple_colormap()
        vmin, vmax_latt = -0.5, 1.5
    else:
        cmap = make_colormap(n_types)
        vmin, vmax_latt = -0.5, n_types + 0.5

    # ------------------------------------------------------------------ #
    # Figure / axes layout
    # ------------------------------------------------------------------ #
    fig_h = max(5, n_rows + 2)
    fig   = plt.figure(figsize=(18, fig_h), facecolor="#f5f5f5")

    if has_params:
        gs = GridSpec(
            1, 3, figure=fig,
            width_ratios=[1, 1.6, 1.1],
            left=0.06, right=0.97, top=0.92, bottom=0.09, wspace=0.40)

        ax_grid   = fig.add_subplot(gs[0, 0])
        ax_energy = fig.add_subplot(gs[0, 1])

        if has_jpair:
            # Right column: coupling matrix (top) + info text (bottom)
            # The matrix subspec uses a 2-column inner grid so the colourbar
            # gets its own thin column.
            gs_right = GridSpecFromSubplotSpec(
                2, 2,
                subplot_spec=gs[0, 2],
                height_ratios=[max(1.2, n_types * 0.6), 1.0],
                width_ratios=[1, 0.07],
                hspace=0.55, wspace=0.08)

            ax_matrix = fig.add_subplot(gs_right[0, 0])
            ax_cbar   = fig.add_subplot(gs_right[0, 1])
            ax_info   = fig.add_subplot(gs_right[1, :])
        else:
            # No Jpair: just an info text panel
            ax_matrix = None
            ax_cbar   = None
            ax_info   = fig.add_subplot(gs[0, 2])
    else:
        gs = GridSpec(
            1, 2, figure=fig,
            width_ratios=[1, 1.6],
            left=0.06, right=0.97, top=0.92, bottom=0.09, wspace=0.35)
        ax_grid   = fig.add_subplot(gs[0, 0])
        ax_energy = fig.add_subplot(gs[0, 1])
        ax_matrix = ax_cbar = ax_info = None

    mode_tag = (f"  [simple  {fps:.0f} fps]" if simple
                else f"  [{fps:.1f} fps]")
    fig.suptitle(algo_name + mode_tag, fontsize=9, y=0.995, style="italic")

    # ------------------------------------------------------------------ #
    # Lattice panel
    # ------------------------------------------------------------------ #
    def to_grid(state):
        arr = np.array(state, dtype=float).reshape(n_rows, n_cols)
        if simple:
            arr = (arr > 0).astype(float)
        return arr

    im = ax_grid.imshow(
        to_grid(states[0]), cmap=cmap, vmin=vmin, vmax=vmax_latt,
        interpolation="nearest", aspect="equal")

    ax_grid.set_xticks(range(n_cols))
    ax_grid.set_xticklabels(range(1, n_cols + 1), fontsize=7)
    if n_rows > 1:
        ax_grid.set_yticks(range(n_rows))
        ax_grid.set_yticklabels(range(1, n_rows + 1), fontsize=7)
        ax_grid.set_ylabel("Row", fontsize=9)
    else:
        ax_grid.set_yticks([])
    ax_grid.set_xlabel("Site", fontsize=9)
    ax_grid.set_title("Step 0", fontsize=10)

    ax_grid.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax_grid.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax_grid.grid(which="minor", color="white", linewidth=1.5)
    ax_grid.tick_params(which="minor", length=0)

    if simple:
        legend_handles = [
            Patch(facecolor="#dddddd", label="0  (hole)"),
            Patch(facecolor="#4c72b0", label="n>0  (particle)"),
        ]
    else:
        legend_handles = [Patch(facecolor="#1a1a1a", label="0  (empty)")]
        for k in range(1, n_types + 1):
            legend_handles.append(
                Patch(facecolor=_TAB10[(k - 1) % len(_TAB10)],
                      label=f"{k}  (type {k})"))
    ax_grid.legend(handles=legend_handles,
                   loc="upper left", bbox_to_anchor=(1.01, 1),
                   fontsize=8, framealpha=0.9,
                   title="Particle", title_fontsize=8)

    # ------------------------------------------------------------------ #
    # Energy panel
    # ------------------------------------------------------------------ #
    ax_energy.set_facecolor("#fafafa")
    emin = min(energies)
    emax = max(energies)
    epad = max((emax - emin) * 0.15, 0.05)
    ax_energy.set_xlim(0, max(steps) if steps else 1)
    ax_energy.set_ylim(emin - epad, emax + epad)
    ax_energy.set_xlabel("Step", fontsize=9)
    ax_energy.set_ylabel("Energy", fontsize=9)
    ax_energy.set_title("System energy", fontsize=10)
    ax_energy.grid(True, alpha=0.3, linestyle="--")

    e_line, = ax_energy.plot([], [], color="#4c72b0", lw=1.5, label="E(t)")
    e_dot,  = ax_energy.plot([], [], "o", color="#c44e52", ms=6, zorder=5)
    e_text  = ax_energy.text(
        0.97, 0.05, "", transform=ax_energy.transAxes,
        ha="right", va="bottom", fontsize=8, color="#444444")
    ax_energy.legend(fontsize=8, loc="upper right")

    xs: list = []
    ys: list = []

    # ------------------------------------------------------------------ #
    # Right panel: coupling matrix + info
    # ------------------------------------------------------------------ #
    if ax_matrix is not None:
        draw_coupling_matrix(ax_matrix, ax_cbar, jpair_matrix, n_types)

    if ax_info is not None:
        draw_info_panel(ax_info, scalars, coupling_formula, field_formula)
        ax_info.set_title(
            "Parameters" if not has_jpair else "Model",
            fontsize=9, pad=4)

    # ------------------------------------------------------------------ #
    # Animation
    # ------------------------------------------------------------------ #
    # Use frames=None (infinite counter) so FuncAnimation never auto-stops
    # its timer.  On macOS/Tkinter, having the timer stop itself (repeat=False)
    # races with the Tk event loop teardown and causes a SIGSEGV inside
    # mach_vm_allocate_kernel.  Instead we clamp the frame index to the last
    # frame so the display freezes after all frames are shown, and we let the
    # user close the window explicitly.  blit is always False to avoid the
    # extra blitting code path that triggers the same crash on some macOS
    # matplotlib builds.
    def update(frame: int):
        i = min(frame, n_frames - 1)   # freeze on last frame once done

        im.set_data(to_grid(states[i]))
        ax_grid.set_title(f"Step {steps[i]}", fontsize=10)

        if frame < n_frames:
            xs.append(steps[i])
            ys.append(energies[i])
            e_line.set_data(xs, ys)
            e_dot.set_data([xs[-1]], [ys[-1]])
            e_text.set_text(f"E = {energies[i]:.4f}")

        return im, e_line, e_dot, e_text

    ani = animation.FuncAnimation(
        fig, update,
        frames=None,        # infinite — timer never auto-stops
        interval=delay_ms,
        blit=False,         # avoid blit teardown crash on macOS
        cache_frame_data=False)

    # On window close, hard-exit the process immediately.  os._exit(0) bypasses
    # all Python atexit handlers and GC finalizers, which is what triggers the
    # macOS SIGSEGV: Python tries to GC matplotlib Tk objects after Tk has
    # already shut down.
    def on_close(_event):
        os._exit(0)

    fig.canvas.mpl_connect("close_event", on_close)

    try:
        plt.show()
    except Exception:
        pass

    os._exit(0)


if __name__ == "__main__":
    main()
