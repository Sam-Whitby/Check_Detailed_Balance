#!/usr/bin/env python3
"""
DetailedBalanceChecker – Python visualisation
Usage: python3 show_report.py <report.json>
Saves a high-resolution PNG and prints its path to stdout.
"""

import sys, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from collections import defaultdict, deque

GREEN = '#1a7a1a'
RED   = '#8b0000'
BLUE  = '#3a6bb5'

# Panels are skipped above these state-count thresholds to prevent
# matplotlib producing images that exceed the 2^16 pixel-per-axis limit.
MAX_TREE_STATES   = 24   # individual decision-tree panel
MAX_MATRIX_STATES = 50   # full symbolic transition-matrix panel

# ── main ─────────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) < 2:
        print("Usage: show_report.py <report.json>", file=sys.stderr)
        sys.exit(1)

    with open(sys.argv[1]) as f:
        data = json.load(f)

    name      = data["name"]
    passed    = data["pass"]
    kl        = data["kl"]
    states    = data["allStates"]
    sim_freq  = data["simFreq"]
    boltzmann = data["boltzmann"]
    db_pairs  = data["dbPairs"]
    mat_raw   = data["matrix"]
    tree_raw  = data["treeData"]
    alg_code  = data.get("algCode", "")

    n          = len(states)
    code_lines = alg_code.rstrip('\n').split('\n')
    n_code     = len(code_lines)

    show_trees  = n <= MAX_TREE_STATES
    show_matrix = n <= MAX_MATRIX_STATES

    # ── Figure sizing ─────────────────────────────────────────────────────────
    # Each code line at fontsize 8 needs ~0.155 in (including line spacing)
    BANNER_H = 0.55
    TOP_H    = 4.0
    TREE_H   = 7.0  if show_trees  else 0.8   # collapsed placeholder when skipped
    MAT_H    = 2.8  if show_matrix else 0.8
    CODE_H   = max(3.5, n_code * 0.155 + 0.9)
    FIG_H    = TOP_H + TREE_H + MAT_H + CODE_H + BANNER_H

    FIG_W    = 20.0

    fig = plt.figure(figsize=(FIG_W, FIG_H))
    fig.patch.set_facecolor('#f4f6f9')

    # ── Banner ────────────────────────────────────────────────────────────────
    banner_col = GREEN if passed else RED
    verdict    = '✓  DETAILED BALANCE:  PASS' if passed else '✗  DETAILED BALANCE:  FAIL'
    banner_y   = 1.0 - (BANNER_H / 2) / FIG_H
    fig.text(0.5, banner_y,
             f'{name}     {verdict}',
             ha='center', va='center', fontsize=14, fontweight='bold', color='white',
             bbox=dict(boxstyle='round,pad=0.45', facecolor=banner_col, alpha=1.0))

    # ── GridSpec ──────────────────────────────────────────────────────────────
    top_frac = 1.0 - BANNER_H / FIG_H - 0.005
    gs = gridspec.GridSpec(
        4, 3, figure=fig,
        height_ratios=[TOP_H, TREE_H, MAT_H, CODE_H],
        hspace=0.38, wspace=0.28,
        top=top_frac, bottom=0.008, left=0.04, right=0.97)

    ax_sc   = fig.add_subplot(gs[0, :2])
    ax_db   = fig.add_subplot(gs[0, 2])
    ax_mat  = fig.add_subplot(gs[2, :])
    ax_code = fig.add_subplot(gs[3, :])

    # ── Tree row ──────────────────────────────────────────────────────────────
    if show_trees:
        gs_trees = gridspec.GridSpecFromSubplotSpec(
            1, n, subplot_spec=gs[1, :], wspace=0.25)
        ax_tr = [fig.add_subplot(gs_trees[0, i]) for i in range(n)]
        for i in range(n):
            draw_tree(ax_tr[i], states[i],
                      tree_raw[i] if i < len(tree_raw) else [], states)
    else:
        ax_skip = fig.add_subplot(gs[1, :])
        ax_skip.axis('off')
        ax_skip.text(0.5, 0.5,
                     f'Decision-tree panel omitted — {n} states exceeds display threshold '
                     f'({MAX_TREE_STATES}).\nAll states: ' +
                     ', '.join(str(s) for s in states[:40]) +
                     (f' … (+{n-40} more)' if n > 40 else ''),
                     ha='center', va='center', transform=ax_skip.transAxes,
                     fontsize=10, color='#555', wrap=True,
                     bbox=dict(boxstyle='round,pad=0.5', facecolor='#f0f0f0',
                               edgecolor='#aaa', lw=0.8))

    # ── Draw remaining panels ─────────────────────────────────────────────────
    draw_scatter(ax_sc, states, sim_freq, boltzmann, kl)
    draw_db_table(ax_db, db_pairs)
    if show_matrix:
        draw_matrix(ax_mat, states, mat_raw)
    else:
        ax_mat.axis('off')
        ax_mat.text(0.5, 0.5,
                    f'Transition-matrix panel omitted — {n} states exceeds display threshold '
                    f'({MAX_MATRIX_STATES}).',
                    ha='center', va='center', transform=ax_mat.transAxes,
                    fontsize=10, color='#555',
                    bbox=dict(boxstyle='round,pad=0.5', facecolor='#f0f0f0',
                              edgecolor='#aaa', lw=0.8))
    draw_code_panel(ax_code, alg_code)

    # ── Save (with fallback if image exceeds matplotlib's pixel limit) ────────
    out_path = sys.argv[1].replace('.json', '.png')
    try:
        fig.savefig(out_path, dpi=250, bbox_inches='tight',
                    facecolor=fig.get_facecolor())
    except (ValueError, OverflowError) as exc:
        # Image still too large even after panel suppression — fall back to
        # a compact summary-only figure (scatter + DB table + code).
        print(f'  Warning: full figure too large ({exc}); saving compact summary.',
              file=sys.stderr)
        plt.close(fig)
        fig2, axes2 = plt.subplots(1, 2, figsize=(FIG_W, TOP_H + BANNER_H + 0.5))
        fig2.patch.set_facecolor('#f4f6f9')
        fig2.text(0.5, 0.97, f'{name}     {verdict}',
                  ha='center', va='top', fontsize=13, fontweight='bold', color='white',
                  bbox=dict(boxstyle='round,pad=0.4', facecolor=banner_col, alpha=1.0))
        draw_scatter(axes2[0], states, sim_freq, boltzmann, kl)
        draw_db_table(axes2[1], db_pairs)
        fig2.tight_layout(rect=[0, 0, 1, 0.94])
        fig2.savefig(out_path, dpi=150, bbox_inches='tight',
                     facecolor=fig2.get_facecolor())
    print(out_path)


# ── Scatter plot: simulated vs Boltzmann ──────────────────────────────────────

def draw_scatter(ax, states, sim_freq, boltzmann, kl):
    x  = np.array(boltzmann, dtype=float)
    y  = np.array(sim_freq,  dtype=float)
    hi = max(x.max(), y.max()) * 1.18
    lo = 0.0

    # y = x reference line
    ax.plot([lo, hi], [lo, hi], color='#555', linestyle='--', lw=1.8,
            alpha=0.6, label='y = x  (perfect detailed balance)', zorder=1)

    # One point per state
    palette = plt.cm.tab10(np.linspace(0, 0.7, len(states)))
    for i, s in enumerate(states):
        ax.scatter(x[i], y[i], s=160, color=palette[i], zorder=5,
                   edgecolors='white', linewidths=1.0)
        ax.annotate(f'  State {s}', (x[i], y[i]),
                    textcoords='offset points', xytext=(8, 4),
                    fontsize=10, color=palette[i], fontweight='bold')

    # Pearson r
    if len(x) > 1 and x.std() > 0 and y.std() > 0:
        r = float(np.corrcoef(x, y)[0, 1])
        if r > 0.999:
            r_col, r_note = GREEN, '(excellent)'
        elif r > 0.98:
            r_col, r_note = '#886600', '(moderate)'
        else:
            r_col, r_note = RED, '(poor)'
        ax.text(0.04, 0.96,
                f'Pearson  r = {r:.6f}  {r_note}',
                transform=ax.transAxes, fontsize=10, fontweight='bold',
                color=r_col, va='top')

    ax.set_xlabel('Boltzmann probability (theoretical)', fontsize=11)
    ax.set_ylabel('Simulated frequency (MCMC)', fontsize=11)
    ax.set_title('MCMC vs Boltzmann — state frequencies', fontsize=12, fontweight='bold')
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect('equal', adjustable='box')
    ax.legend(fontsize=9, loc='lower right')
    ax.grid(alpha=0.3, zorder=0)
    ax.set_facecolor('#fafbff')

    kl_col = GREEN if kl < 0.02 else RED
    ax.text(0.97, 0.04,
            f'KL divergence = {kl:.6f}',
            transform=ax.transAxes, ha='right', va='bottom',
            fontsize=9, color=kl_col, fontweight='bold')


# ── Decision tree ─────────────────────────────────────────────────────────────

def _bkey(bits):
    return 'root' if not bits else 'b' + ''.join(map(str, bits))

def build_tree_graph(state, paths):
    """Return (nodes, edges, edge_labels) for one starting state."""
    nodes    = {}   # id -> {label, color, round}
    edges    = []   # list of (u, v)
    edge_lbl = {}   # (u, v) -> bit label string

    nodes['root'] = dict(label=f'S={state}', color=BLUE, round=True)

    all_pfx  = set()
    leaf_pfx = set()
    for path in paths:
        bits = tuple(path['bits'])
        leaf_pfx.add(bits)
        for k in range(len(bits)):
            all_pfx.add(bits[:k+1])

    for pfx in all_pfx - leaf_pfx:
        nodes[_bkey(pfx)] = dict(label='?', color='#555', round=True)
    for pfx in leaf_pfx:
        nodes[_bkey(pfx)] = dict(label='↓', color='#999', round=True)

    for pfx in all_pfx:
        parent = _bkey(pfx[:-1])
        child  = _bkey(pfx)
        edges.append((parent, child))
        edge_lbl[(parent, child)] = str(pfx[-1])

    for path in paths:
        bits    = tuple(path['bits'])
        leaf_id = _bkey(bits)
        for i, out in enumerate(path['outcomes']):
            oid   = leaf_id + f'o{i}'
            dest  = out['state']
            prob  = out['probStr']
            color = '#2a9e48' if dest != state else '#d47b10'
            nodes[oid] = dict(label=f'→{dest}\n{prob}', color=color, round=False)
            edges.append((leaf_id, oid))

    return nodes, edges, edge_lbl


def bfs_layout(edges, root):
    children = defaultdict(list)
    for u, v in edges:
        children[u].append(v)
    levels  = defaultdict(list)
    visited = {root}
    q       = deque([(root, 0)])
    while q:
        node, lv = q.popleft()
        levels[lv].append(node)
        for ch in children[node]:
            if ch not in visited:
                visited.add(ch)
                q.append((ch, lv + 1))
    max_lv = max(levels) if levels else 0
    pos    = {}
    for lv, lvnodes in levels.items():
        nn = len(lvnodes)
        for i, nd in enumerate(lvnodes):
            pos[nd] = ((i + 0.5) / nn,
                       1.0 - lv / (max_lv + 1) if max_lv else 0.85)
    return pos


def draw_tree(ax, state, paths, all_states):
    ax.axis('off')
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.05, 1.10)
    ax.set_facecolor('#f7f8fc')
    ax.set_title(f'State {state}', fontsize=11, fontweight='bold', pad=5)

    if not paths:
        ax.text(0.5, 0.5, '(no paths)', ha='center', va='center',
                transform=ax.transAxes, color='gray', fontsize=10)
        return

    nodes, edges, edge_lbl = build_tree_graph(state, paths)
    pos = bfs_layout(edges, 'root')

    # Draw edges first (lower zorder)
    for u, v in edges:
        if u not in pos or v not in pos:
            continue
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        ax.annotate(
            '', xy=(x1, y1), xytext=(x0, y0),
            xycoords='data', textcoords='data',
            arrowprops=dict(arrowstyle='-|>', color='#777',
                            lw=1.1, mutation_scale=12,
                            shrinkA=12, shrinkB=12),
            annotation_clip=False, zorder=1)
        lbl = edge_lbl.get((u, v))
        if lbl is not None:
            ax.text((x0 + x1) / 2 + 0.03, (y0 + y1) / 2,
                    lbl, fontsize=9, color='#334499',
                    fontweight='bold', ha='left', va='center',
                    clip_on=False, zorder=2)

    # Draw nodes as auto-sized text boxes (clip_on=False so label never gets cut)
    for nid, info in nodes.items():
        if nid not in pos:
            continue
        x, y = pos[nid]
        # Outcome nodes (rounded rectangle) carry two-line labels
        bs = 'round,pad=0.38' if info['round'] else 'round,pad=0.32'
        ax.text(x, y, info['label'],
                ha='center', va='center',
                fontsize=8, fontweight='bold', color='white',
                clip_on=False, zorder=5,
                linespacing=1.4,
                bbox=dict(boxstyle=bs,
                          facecolor=info['color'],
                          edgecolor='white', lw=1.2))


# ── Detailed balance pair table ───────────────────────────────────────────────

def draw_db_table(ax, db_pairs):
    ax.axis('off')
    ax.set_title('Detailed Balance\nCheck', fontsize=11, fontweight='bold', pad=4)
    if not db_pairs:
        ax.text(0.5, 0.5, '(no pairs)', ha='center', va='center',
                transform=ax.transAxes)
        return

    rows    = []
    rcolors = []
    for p in db_pairs:
        ok = p['pass']
        rows.append([str(p['i']), str(p['j']),
                     '✓ PASS' if ok else '✗ FAIL'])
        bg = '#d4f0d4' if ok else '#f0d4d4'
        rcolors.append(['#ffffff', '#ffffff', bg])

    tbl = ax.table(cellText=rows,
                   colLabels=['i', 'j', 'Result'],
                   cellColours=rcolors, cellLoc='center', loc='center')
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(11)
    tbl.scale(1.0, 1.8)
    for (r, c), cell in tbl.get_celld().items():
        if r == 0:
            cell.set_facecolor('#c8d8f0')
            cell.set_text_props(fontweight='bold')
        if c == 2 and r > 0:
            cell.set_text_props(fontweight='bold')


# ── Transition matrix (full width, full expressions) ──────────────────────────

def _abbrev(s):
    """Light abbreviation: replace Mathematica escapes, keep as much as useful."""
    s = s.replace('\\[Beta]', 'β').replace('Piecewise', 'Pw')
    s = s.replace('E^(-β', 'e^{-β').replace('E^(β', 'e^{β')
    return s


def draw_matrix(ax, states, mat_raw):
    ax.axis('off')
    ax.set_title('Symbolic Transition Matrix  T[i→j]'
                 '   (Piecewise = Metropolis acceptance prob; β kept symbolic)',
                 fontsize=11, fontweight='bold', pad=5)

    n     = len(states)
    cells = [['' for _ in range(n)] for _ in range(n)]
    for entry in mat_raw:
        fi = states.index(entry['from'])
        ti = states.index(entry['to'])
        cells[fi][ti] = _abbrev(entry['str'])

    tbl = ax.table(
        cellText=cells,
        rowLabels=[f'{s} →' for s in states],
        colLabels=[f'→ {s}' for s in states],
        cellLoc='center', loc='center')
    tbl.auto_set_font_size(True)
    tbl.scale(1.0, 2.2)
    for (r, c), cell in tbl.get_celld().items():
        if r == 0 or c == -1:
            cell.set_facecolor('#c8d8f0')
            cell.set_text_props(fontweight='bold')
        elif r > 0 and c >= 0 and r - 1 == c:
            cell.set_facecolor('#fffde0')


# ── Full algorithm source code ────────────────────────────────────────────────

def draw_code_panel(ax, code):
    ax.axis('off')
    ax.set_facecolor('#fffff8')
    ax.set_title('Algorithm Under Test  (full source)',
                 fontsize=11, fontweight='bold', pad=5)
    # Position at top-left of axes; clip_on=False so tall code is never cut
    ax.text(0.012, 0.988, code,
            transform=ax.transAxes,
            fontsize=8, fontfamily='monospace',
            verticalalignment='top', horizontalalignment='left',
            clip_on=False,
            linespacing=1.45,
            bbox=dict(boxstyle='round,pad=0.55',
                      facecolor='#fffff0', edgecolor='#ccccaa',
                      alpha=0.97, lw=0.8))


if __name__ == '__main__':
    main()
