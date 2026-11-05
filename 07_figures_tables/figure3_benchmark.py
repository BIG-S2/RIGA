#!/usr/bin/env python3
"""Benchmark figure (main Figure 3 from 2 November 2026): same data and style as fig2_benchmark_updated_20261027,
panel c retitled 'association methods', descriptive gray text removed. Reads the saved panel sources, no recomputation."""
from pathlib import Path
import os, sys, json, hashlib
os.environ.setdefault('MPLCONFIGDIR', '/tmp/riga_fig3_mpl')
import numpy as np, pandas as pd
REV = Path(__file__).resolve().parents[1]
SRC = REV/'manuscript/integration_20261027/benchmark_assets'
sys.path.insert(0, str(REV/'scripts'))
import fig_plasmode_benchmark_ng_20261014 as F
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.ticker as mticker
single = pd.read_csv(SRC/'figure_panel_a_source.csv'); rv10 = pd.read_csv(SRC/'figure_panel_b_source.csv'); plotext = pd.read_csv(SRC/'figure_panel_c_source.csv')
keys = ['oracle', 'strat_project', 'riga_projectkey', 'pooled_project']
assert set(single.key) == set(keys) and len(plotext) == 18
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans'], 'font.size': 8, 'text.color': F.INK, 'axes.labelcolor': F.INK, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none', 'axes.unicode_minus': False, 'mathtext.fontset': 'dejavusans', 'savefig.facecolor': 'white'})
W, H = 183, 180
fig = plt.figure(figsize=(W/25.4, H/25.4), dpi=300, facecolor='white')
def axmm(x, y, w, h): return fig.add_axes([x/W, y/H, w/W, h/H])
def txt(x, y, s, **kw): return fig.text(x/W, y/H, s, **kw)
a = axmm(17, 101, 66, 55)
F.base_style(a); a.set_yscale('log'); a.set_ylim(.25, 620); a.set_xlim(.1, .65)
a.set_xticks(np.arange(.1, .7, .1)); a.set_yticks([.3, 1, 3, 10, 30, 100, 300]); a.set_yticklabels(['0.3', '1', '3', '10', '30', '100', '300']); a.yaxis.set_minor_locator(mticker.NullLocator()); a.grid(axis='y', color=F.GRID, lw=.5); a.axhline(1, color='#888888', lw=.8, ls=(0, (3, 2)))
for r in single.itertuples():
    color = F.STYLES[r.key][1]
    a.plot([r.power]*2, [r.lo, r.hi], color=color, lw=.95)
    a.plot([r.power-1.96*r.power_se, r.power+1.96*r.power_se], [r.ratio]*2, color=color, lw=.95)
    F.point(a, r.power, r.ratio, r.key)
a.set_xlabel('Power', fontsize=8.5, labelpad=4); a.set_ylabel('Type-I error / nominal', fontsize=8.5, labelpad=4)
b = axmm(130, 101, 49, 55); F.panel_b(b, rv10)
for x, label, title in [(3, 'a', 'Single-variant tests'), (95, 'b', 'Rare-variant tests')]:
    txt(x, 175, label, fontsize=11, fontweight='bold', va='top')
    txt(17 if label == 'a' else 103, 175, title, fontsize=9, fontweight='bold', va='top')
handles = []
for key in keys:
    label, color, marker, size, hollow = F.STYLES[key]
    if key == 'riga_projectkey': label = 'Case–control matched'
    handles.append(Line2D([], [], ls='none', marker=marker, ms=size, mfc='white' if hollow else color, mec=color if hollow else 'white', color=color, label=label))
fig.legend(handles=handles, ncol=4, loc='center', frameon=False, bbox_to_anchor=(.51, 88/H), fontsize=7.6, handletextpad=.4, columnspacing=1.4)
txt(3, 78, 'c', fontsize=11, fontweight='bold', va='top')
txt(17, 78, 'Calibration across association methods', fontsize=9, fontweight='bold', va='top')
armsty = {'retained_joint': (F.BLUE, 's', False, .12), 'retained_legacy_additive': ('#009E73', 'o', True, -.12)}
for ix, (alpha, alabel) in enumerate([(1e-5, r'$\alpha = 10^{-5}$'), (1e-6, r'$\alpha = 10^{-6}$'), (5e-8, r'$\alpha = 5\times10^{-8}$')]):
    ax = axmm(32+52*ix, 25, 43, 38); F.base_style(ax); ax.set_xscale('log'); ax.set_xlim(.3, 400); ax.set_ylim(-.5, 2.5)
    ax.set_xticks([1, 10, 100]); ax.set_xticklabels(['1', '10', '100']); ax.xaxis.set_minor_locator(mticker.NullLocator()); ax.axvline(1, color='#888888', lw=.8, ls=(0, (3, 2))); ax.grid(axis='x', color=F.GRID, lw=.5)
    ax.spines['left'].set_visible(False); ax.set_yticks([2, 1, 0]); ax.set_yticklabels(['GMMAT', 'GENESIS', 'REGENIE'] if ix == 0 else ['', '', '']); ax.tick_params(axis='y', length=0)
    ax.set_title(alabel, fontsize=8, pad=6); ax.set_xlabel('Type-I error / nominal', fontsize=7.6, labelpad=4)
    for ti, tool in enumerate(['gmmat', 'genesis', 'regenie']):
        for arm, (color, marker, hollow, dy) in armsty.items():
            z = plotext[(plotext.tool == tool) & (plotext.arm == arm) & np.isclose(plotext.alpha, alpha, atol=0, rtol=1e-12)].iloc[0]
            y = 2-ti+dy; ax.plot([z.lo, z.hi], [y, y], color=color, lw=1.05)
            ax.plot(z.ratio, y, ls='none', marker=marker, ms=4.8, mfc='white' if hollow else color, mec=color, mew=.9)
handles = [Line2D([], [], ls='none', marker=m, ms=4.8, mfc='white' if h else c, mec=c, label='Joint-unit intercepts' if ar == 'retained_joint' else 'Historical additive coding') for ar, (c, m, h, dy) in armsty.items()]
fig.legend(handles=handles, ncol=2, loc='center', frameon=False, bbox_to_anchor=(.53, 12/H), fontsize=7.6, handletextpad=.4, columnspacing=2)
out = REV/'manuscript/figures/fig3_benchmark_20261102'
for suffix in ('.pdf', '.svg', '.png'): fig.savefig(out.with_suffix(suffix), dpi=300, facecolor='white')
plt.close(fig)
manifest = {'figure': str(out), 'derived_from': 'fig2_benchmark_updated_20261027 (same panel sources)', 'inputs': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in [SRC/'figure_panel_a_source.csv', SRC/'figure_panel_b_source.csv', SRC/'figure_panel_c_source.csv']}, 'changes': ['panel c title: association methods', 'gray subtitles, Nominal label and bottom note removed', 'renumbered as main Figure 3']}
Path(str(out)+'.provenance.json').write_text(json.dumps(manifest, indent=2)); print('Wrote', out)
