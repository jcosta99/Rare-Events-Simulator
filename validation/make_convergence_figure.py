#!/usr/bin/env python3

"""Draw the approach of the simulated open SSEP to its large-L closed form.

Derrida's result for the open symmetric exclusion process is exact only to
leading order in the system size, and the engine reports the generating
function per bond, so the combination with a finite limit is ``(L-1) lambda``.
Both panels use it: the left shows the whole curve at several sizes collapsing
onto the limit, the right follows a few counting fields as the size grows.
"""

import glob
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

NNN = Path(__file__).resolve().parent.parent
SIZES = (8, 16, 32, 64, 128)
FIELDS = (-3.0, -1.5, 1.5, 3.0)
TRACKED = (1.0, 2.0, 3.0)

# L is an ordered quantity, so the sizes take a single-hue ramp rather than
# categorical hues; the counting fields on the right are categorical.
RAMP = ('#86b6ef', '#5598e7', '#2a78d6', '#1c5cab', '#0d366b')
SLOTS = ('#2a78d6', '#eb6834', '#1baf7a')
INK, INK_SOFT, GRID_INK = '#0b0b0b', '#52514e', '#d9d8d4'

plt.rcParams.update({
    'figure.facecolor': 'white', 'axes.facecolor': 'white', 'font.size': 9,
    'axes.labelsize': 9, 'legend.fontsize': 8, 'xtick.labelsize': 8.5,
    'ytick.labelsize': 8.5, 'axes.edgecolor': GRID_INK, 'text.color': INK,
    'axes.labelcolor': INK, 'xtick.color': INK_SOFT, 'ytick.color': INK_SOFT,
    'axes.grid': True, 'grid.color': GRID_INK, 'grid.linewidth': 0.6,
    'axes.axisbelow': True, 'legend.frameon': False,
    'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02})


def limit(s, rho_a=0.5, rho_b=0.5):
    """The large-L closed form, already multiplied by (L-1)."""
    omega = (rho_a * (math.exp(s) - 1) + rho_b * (math.exp(-s) - 1)
             + rho_a * rho_b * (math.exp(s) - 1) * (math.exp(-s) - 1))
    if omega >= 0:
        return math.asinh(math.sqrt(omega)) ** 2
    return -math.asin(math.sqrt(-omega)) ** 2


def measured():
    """Every equal-reservoir SSEP run in ../data, as (L, s, (L-1) lambda)."""
    rows = []
    for path in glob.glob(str(NNN / 'data' / 'WASEP_O_*_E0.0_*.csv')):
        meta = json.loads(Path(path.replace('.csv', '.json')).read_text())
        if meta['walkers'] != 2000 or abs(meta['gamma_left_removal'] - 1) > 1e-9:
            continue
        value = float(pd.read_csv(path)['CGFCloningNormalization'].iloc[-1])
        rows.append({'L': meta['length'], 's': meta['bias_s'],
                     'scaled': (meta['length'] - 1) * value})
    return pd.DataFrame(rows)


def main():
    table = measured()
    figure, (left, right) = plt.subplots(1, 2, figsize=(9.2, 3.4))

    fine = np.linspace(-3.2, 3.2, 400)
    left.plot(fine, [limit(s) for s in fine], color=INK, linestyle='--',
              linewidth=1.5, zorder=5, label=r'$L\to\infty$')
    for colour, L in zip(RAMP, SIZES):
        sub = table[table['L'] == L].sort_values('s')
        if sub.empty:
            continue
        left.plot(sub['s'], sub['scaled'], color=colour, marker='o',
                  markersize=4.5, linewidth=1.6, label=f'$L={L}$')
    left.set(xlabel='$s$', ylabel=r'$(L-1)\,\lambda(s)$')
    left.legend(loc='upper center', ncol=2)

    for colour, s in zip(SLOTS, TRACKED):
        sub = table[np.isclose(table['s'], s)].sort_values('L')
        if len(sub) < 3:
            continue
        right.plot(sub['L'], sub['scaled'], color=colour, marker='o',
                   markersize=5, linewidth=1.8, label=f'$s={s:g}$')
        right.axhline(limit(s), color=colour, linestyle=':', linewidth=1.2)
    right.set(xscale='log', xlabel='$L$', ylabel=r'$(L-1)\,\lambda(s)$')
    right.set_xticks(SIZES)
    right.set_xticklabels([str(L) for L in SIZES])
    right.xaxis.set_minor_formatter(plt.NullFormatter())
    right.legend(loc='center right')

    for axes in (left, right):
        for side in ('top', 'right'):
            axes.spines[side].set_visible(False)
    figure.text(0.005, -0.05,
                r'Open SSEP, $E=0$, equal reservoirs, $M=2000$, '
                r'$t=2\times10^4$. Dotted lines are the limit at each $s$.',
                ha='left', fontsize=8, color=INK_SOFT)
    figure.tight_layout()

    out = NNN / 'validation' / 'figures'
    out.mkdir(parents=True, exist_ok=True)
    for suffix in ('pdf', 'png'):
        figure.savefig(out / f'ssep_large_L_convergence.{suffix}', dpi=200)
    print('wrote validation/figures/ssep_large_L_convergence.{pdf,png}')

    for s in TRACKED:
        sub = table[np.isclose(table['s'], s)].sort_values('L')
        gap = (sub['scaled'].iloc[-1] - limit(s)) / abs(limit(s)) * 100
        print(f'  s={s:>5g}: limit {limit(s):.5f}, '
              f'L={sub["L"].iloc[-1]:.0f} gives {sub["scaled"].iloc[-1]:.5f} '
              f'({gap:+.2f}%)')


if __name__ == '__main__':
    main()
