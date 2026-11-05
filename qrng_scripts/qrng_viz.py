r"""
## Visualization and Report Styling

**Purpose:** Loads the shared plotting functions and publication-style theme used by the analysis summaries. It does not connect to the Arduino or read any data.

**Usage:** Run this cell first whenever the notebook kernel is started or restarted. There are no settings a normal user needs to change. If a later cell says a plotting function is not defined, run this cell again.

Source: notebook cell 1 (code copied unchanged unless noted below).
"""

import textwrap
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from qrng_common import SENSORS, Fs   # names the notebook cell used from the setup cell

### SUMMARY VISUALIZATION ###
# Publication-style figures for the analysis-cell summaries. Same call signatures the
# analysis cells already use -- only the rendering is journal-style: white ground, serif
# type, boxed axes with inward ticks, and lettered panels.

import textwrap
import matplotlib.ticker as mticker
from matplotlib.colors import LinearSegmentedColormap, to_rgb

VIZ_WIDTH_IN  = 7.16     # double-column figure width (in); ncols=1 uses VIZ_COL_IN
VIZ_COL_IN    = 3.50     # single-column width (in)
VIZ_DPI       = 150      # on-screen dpi
VIZ_SAVE_DPI  = 600      # raster save dpi (ignored for pdf/svg)
VIZ_FORMAT    = 'pdf'    # vector copy saved next to the png: 'pdf', 'svg' or 'png'
VIZ_HATCH     = False    # hatch bar fills so the figure survives grayscale printing
VIZ_SERIF     = ['Times New Roman', 'Nimbus Roman', 'DejaVu Serif']

# Okabe-Ito: the CVD-safe qualitative set used across the sciences. Fixed order, never
# cycled; every series also carries a marker shape so identity is never colour-alone.
_VIZ_SERIES  = ['#0072B2', '#D55E00', '#009E73', '#E69F00',
                '#CC79A7', '#56B4E9', '#F0E442', '#000000']
_VIZ_MARKERS = ['o', 's', '^', 'D', 'v', 'P', 'X', '*']
_VIZ_HATCHES = ['', '///', '\\\\\\', 'xxx', '...', '+++', 'ooo', '**']
VIZ_STATUS   = dict(good='#000000', warning='#D55E00', serious='#D55E00', critical='#B2182B')

_VIZ_INK = dict(ink='#000000', ink2='#333333', muted='#606060',
                rule='#000000', grid='#D9D9D9', surface='#FFFFFF')

_VIZ_RC = {
    'font.family': 'serif', 'font.serif': VIZ_SERIF, 'mathtext.fontset': 'stix',
    'font.size': 14, 'axes.labelsize': 14, 'axes.titlesize': 14,
    'xtick.labelsize': 12, 'ytick.labelsize': 12, 'legend.fontsize': 11,
    'axes.linewidth': 0.6, 'axes.edgecolor': '#000000', 'axes.labelcolor': '#000000',
    'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
    'xtick.major.size': 3.0, 'ytick.major.size': 3.0,
    'xtick.minor.size': 1.6, 'ytick.minor.size': 1.6,
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.minor.width': 0.5, 'ytick.minor.width': 0.5,
    'xtick.color': '#000000', 'ytick.color': '#000000',
    'lines.linewidth': 1.0, 'lines.markersize': 3.5, 'legend.frameon': False,
    'legend.handlelength': 1.6, 'legend.columnspacing': 1.2, 'legend.labelspacing': 0.35,
    'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white',
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
}

_VIZ_DEFAULT_H = dict(stat=1.50, bars=2.25, stacked=2.25, dots=2.15, matrix=2.45,
                      lines=2.35, phasor=2.55, table=2.20, text=1.70, averaging_table=2.10)


def _viz_theme():
    T = dict(_VIZ_INK)
    T['series'] = _VIZ_SERIES
    T['markers'] = _VIZ_MARKERS
    T['hatches'] = _VIZ_HATCHES
    # single-hue sequential (light -> dark) and a blue/white/red diverging pair
    T['seq'] = LinearSegmentedColormap.from_list(
        'paper_seq', ['#F7FBFF', '#D9E7F5', '#A9C9E5', '#6BA3D0', '#3577B1', '#12508C'])
    T['div'] = LinearSegmentedColormap.from_list(
        'paper_div', ['#B2182B', '#EF9A7B', '#FFFFFF', '#8DBBDA', '#2166AC'])
    return T


def _viz_axes_style(ax, T, grid=False, gridaxis='y'):
    ylabel = ax.get_ylabel()
    if len(ylabel) > 24 and '\n' not in ylabel:
        ax.set_ylabel(textwrap.fill(ylabel, width=22, break_long_words=False,
                                     break_on_hyphens=False))
    for sp in ax.spines.values():
        sp.set_visible(True)
        sp.set_linewidth(0.6)
        sp.set_color(T['rule'])
    ax.minorticks_on()
    ax.tick_params(which='both', direction='in', top=True, right=True)
    if grid:
        ax.grid(axis=gridaxis, color=T['grid'], lw=0.4, ls=':')
        ax.set_axisbelow(True)


def _viz_panel_title(ax, p, T, pad=3.5):
    letter = p.get('_letter', '').strip()
    title = p.get('title', '')
    if title:
        spec = ax.get_subplotspec()
        width = ax.figure.get_figwidth() * (spec.colspan.stop - spec.colspan.start) / spec.get_gridspec().ncols
        title = textwrap.fill(title, max(18, int(6.5 * (width - (0.45 if letter else 0)))),
                              break_long_words=False)
    if letter:
        ax.annotate(letter, xy=(0, 1), xycoords='axes fraction',
                    xytext=(0, pad + 17 * title.count('\n')), textcoords='offset points',
                    color=T['ink'], fontsize=16, ha='left', va='bottom')
    if title:
        ax.annotate(title, xy=(0, 1), xycoords='axes fraction',
                    xytext=(29 if letter else 0, pad), textcoords='offset points',
                    color=T['ink'], fontsize=14, ha='left', va='bottom')


def _viz_legend(ax, T, n=2, frame=True, **kw):
    # inside the axes with an opaque hairline frame (journal convention): it stays
    # legible wherever it lands, and the panels below reserve headroom for it
    kw.setdefault('loc', 'best')
    leg = ax.legend(frameon=frame, facecolor='white', edgecolor=T['rule'] if frame else 'none',
                    framealpha=1.0, fancybox=False, borderpad=0.35, labelspacing=0.3,
                    labelcolor=T['ink'], fontsize=11, ncols=min(n, 2),
                    handlelength=1.4, columnspacing=1.0, **kw)
    leg.get_frame().set_linewidth(0.4)
    return leg


def _viz_txt_on(color):
    r, g, b = to_rgb(color)
    return '#000000' if (0.299 * r + 0.587 * g + 0.114 * b) > 0.55 else '#FFFFFF'


def _viz_fmt(v, fmt):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return '--'
    return fmt.format(v)


def _viz_rules(ax, T, ys):
    for y, lw in ys:
        ax.plot([0, 1], [y, y], transform=ax.transAxes, color=T['rule'], lw=lw,
                clip_on=False, solid_capstyle='butt')


def _panel_stat(ax, p, T):
    # booktabs-style key-numbers strip: top rule / labels / mid rule / values + notes / rule
    items = p['items']
    ax.set_axis_off()
    _viz_panel_title(ax, p, T, pad=5)
    ncol = len(items)
    wrap = max(16, int(100 / max(ncol, 1)))          # chars that fit one column
    heads = [textwrap.fill(str(it[0]), max(12, int(0.78 * wrap))) for it in items]
    notes = [textwrap.fill(str(it[2]), wrap) if len(it) > 2 and it[2] else '' for it in items]
    nhead = max(txt.count(chr(10)) + 1 for txt in heads)
    nnote = max([txt.count(chr(10)) + 1 for txt in notes if txt] or [0])
    top, bot = 0.98, 0.04
    u = (top - bot) / (0.85 * nhead + 1.35 + 0.66 * nnote)   # one unit of vertical rhythm
    y_head = top - (0.15 + 0.85 * nhead) * u
    y_bot = y_head - (1.35 + 0.66 * nnote) * u
    _viz_rules(ax, T, [(top, 0.9), (y_head, 0.5), (y_bot, 0.9)])
    w = 1.0 / ncol
    for k, it in enumerate(items):
        value = it[1]
        status = it[3] if len(it) > 3 else None
        flagged = status in ('warning', 'serious', 'critical')
        x = k * w + 0.012
        ax.text(x, top - 0.22 * u, heads[k], transform=ax.transAxes, color=T['ink'],
                fontsize=11, va='top', ha='left', linespacing=1.3)
        ax.text(x, y_head - 0.62 * u, f'{value}{"*" if flagged else ""}',
                transform=ax.transAxes, color=VIZ_STATUS[status] if flagged else T['ink'],
                fontsize=16, va='center', ha='left')
        if notes[k]:
            ax.text(x, y_head - 1.12 * u, notes[k], transform=ax.transAxes,
                    color=T['muted'], fontsize=9.5, va='top', ha='left',
                    style='italic', linespacing=1.3)


def _panel_bars(ax, p, T):
    labels = list(p['labels'])
    series = p.get('series') or [(p.get('name', ''), p['values'])]
    fmt = p.get('fmt', '{:.3g}')
    x = np.arange(len(labels), dtype=float)
    n = len(series)
    w = (0.62 if n == 1 else 0.74 / n)
    for k, (name, vals) in enumerate(series):
        vals = np.asarray(vals, dtype=float)
        off = 0.0 if n == 1 else (k - (n - 1) / 2) * (w + 0.015)
        ax.bar(x + off, vals, width=w, facecolor=T['series'][k % 8], edgecolor=T['rule'],
               lw=0.5, label=(name or None), zorder=3,
               hatch=(T['hatches'][k % 8] if VIZ_HATCH else None))
        if p.get('label_bars', True):
            span = max(abs(np.nanmax(vals)), 1e-30)
            for xi, v in zip(x + off, vals):
                if not np.isfinite(v):
                    continue
                ax.text(xi, v + np.sign(v or 1) * 0.02 * span, _viz_fmt(v, fmt),
                        ha='center', va='bottom' if v >= 0 else 'top',
                        fontsize=11, color=T['ink2'])
    if p.get('ref') is not None:
        ax.axhline(p['ref'], color=T['ink'], lw=0.7, ls=(0, (4, 2.5)), zorder=2)
        ax.annotate(p.get('ref_label', 'ref'), xy=(0.99, p['ref']), xytext=(0, 2.5),
                    xycoords=('axes fraction', 'data'), textcoords='offset points',
                    ha='right', va='bottom', fontsize=11, color=T['ink2'])
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.tick_params(axis='x', which='minor', bottom=False, top=False)
    ax.set_ylabel(p.get('ylabel', ''))
    _viz_panel_title(ax, p, T)
    ax.margins(y=0.32 if n > 1 else 0.20)
    _viz_axes_style(ax, T)
    if n > 1:
        _viz_legend(ax, T, n=n)


def _panel_stacked(ax, p, T):
    labels = list(p['labels'])
    x = np.arange(len(labels), dtype=float)
    bottom = np.zeros(len(labels))
    for k, (name, vals) in enumerate(p['series']):
        vals = np.asarray(vals, dtype=float)
        ax.bar(x, vals, bottom=bottom, width=0.58, label=name, facecolor=T['series'][k % 8],
               edgecolor=T['rule'], lw=0.5, zorder=3,
               hatch=(T['hatches'][k % 8] if VIZ_HATCH else None))
        bottom = bottom + vals
    if p.get('label_bars', True):
        fmt = p.get('fmt', '{:.3g}')
        for xi, tot in zip(x, bottom):
            ax.text(xi, tot * 1.02, _viz_fmt(tot, fmt), ha='center', va='bottom',
                    fontsize=11, color=T['ink2'])
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.tick_params(axis='x', which='minor', bottom=False, top=False)
    ax.set_ylabel(p.get('ylabel', ''))
    _viz_panel_title(ax, p, T)
    ax.margins(y=0.34)
    _viz_axes_style(ax, T)
    _viz_legend(ax, T, n=len(p['series']))


def _panel_dots(ax, p, T):
    # horizontal stem on a log axis: decades of range without a fake bar baseline
    labels = list(p['labels'])
    vals = np.asarray(p['values'], dtype=float)
    y = np.arange(len(labels), dtype=float)[::-1]
    logx = p.get('logx', True)
    base = (np.nanmin(vals[vals > 0]) / 3.0) if logx else 0.0
    for yi, v in zip(y, vals):
        ax.hlines(yi, base, v, color=T['ink'], lw=0.6, zorder=3)
        ax.plot([v], [yi], marker='o', ms=4, color=T['series'][0], mec=T['rule'],
                mew=0.5, zorder=4)
        ax.text(v, yi + 0.22, _viz_fmt(v, p.get('fmt', '{:.3g}')), fontsize=11,
                color=T['ink2'], ha='center', va='bottom')
    if logx:
        ax.set_xscale('log')
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.tick_params(axis='y', which='minor', left=False, right=False)
    ax.set_xlabel(p.get('xlabel', ''))
    _viz_panel_title(ax, p, T)
    ax.set_ylim(-0.75, len(labels) - 0.15)
    _viz_axes_style(ax, T, grid=True, gridaxis='x')


def _panel_lines(ax, p, T):
    # series entries are (name, y) -- sharing p['x'] -- or (name, x, y) per series
    for k, entry in enumerate(p['series']):
        name, y = entry[0], np.asarray(entry[-1], dtype=float)
        x = np.asarray(entry[1] if len(entry) == 3 else p.get('x', np.arange(len(y))),
                       dtype=float)
        ax.plot(x, y, lw=1.0, color=T['series'][k % 8], label=name,
                marker=(T['markers'][k % 8] if len(y) <= 40 else None), ms=3.2,
                mfc='white', mew=0.7, zorder=3)
    if p.get('logx'):
        ax.set_xscale('log')
    if p.get('logy'):
        ax.set_yscale('log')
    ax.set_xlabel(p.get('xlabel', ''))
    ax.set_ylabel(p.get('ylabel', ''))
    _viz_panel_title(ax, p, T)
    ax.margins(y=0.22)
    _viz_axes_style(ax, T)
    if len(p['series']) > 1:
        _viz_legend(ax, T, n=len(p['series']))


def _panel_matrix(ax, p, T):
    M = np.asarray(p['M'], dtype=float)
    labels = list(p.get('labels', SENSORS))
    fmt = p.get('fmt', '{:.3g}')
    if p.get('diverging'):
        lim = p.get('vlim', max(1e-30, np.nanmax(np.abs(M))))
        im = ax.imshow(M, cmap=T['div'], vmin=-lim, vmax=lim)
    else:
        im = ax.imshow(M, cmap=T['seq'], vmin=p.get('vmin', 0.0),
                       vmax=p.get('vmax', np.nanmax(M)))
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, _viz_fmt(M[i, j], fmt), ha='center', va='center', fontsize=10,
                    color=_viz_txt_on(im.cmap(im.norm(M[i, j]))))
    ax.set_aspect('auto')
    short = [str(s).split('_')[-1] for s in labels]
    ax.set_xticks(range(len(labels)), short)
    ax.set_yticks(range(len(labels)), short)
    ax.tick_params(which='both', length=0, top=False, right=False)
    ax.minorticks_off()
    for sp in ax.spines.values():
        sp.set_linewidth(0.6)
        sp.set_color(T['rule'])
    _viz_panel_title(ax, p, T)
    if p.get('note'):
        ax.set_xlabel(p['note'], fontsize=10, color=T['muted'], style='italic')
    cb = ax.figure.colorbar(im, ax=ax, fraction=0.045, pad=0.025)
    cb.ax.tick_params(labelsize=10, length=2, width=0.5, direction='out')
    cb.outline.set_linewidth(0.5)
    if 'e' in fmt:
        cb.ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f'{v:.0e}'))
        cb.ax.yaxis.get_offset_text().set_visible(False)
    if p.get('cbar_label'):
        cb.set_label(p['cbar_label'], fontsize=11)


def _panel_phasor(ax, p, T):
    # complex carrier-band phasors: magnitude and phase in one plane
    vals = p['values']              # {label: complex}
    R = max(abs(v) for v in vals.values()) or 1.0
    for r in np.linspace(R / 3, R, 3):
        ax.add_artist(plt.Circle((0, 0), r, fill=False, color=T['grid'], lw=0.5, ls=':',
                                 zorder=1))
        ax.text(r * 0.72, r * 0.72, f'{r:.3g}', color=T['muted'], fontsize=9.5, zorder=2)
    ax.axhline(0, color=T['grid'], lw=0.5, zorder=1)
    ax.axvline(0, color=T['grid'], lw=0.5, zorder=1)
    for k, (name, c) in enumerate(vals.items()):
        col = T['series'][k % 8]
        ax.plot([0, c.real], [0, c.imag], lw=1.0, color=col, zorder=3,
                label=f'{str(name).split("_")[-1]}: {abs(c):.1f}, {np.degrees(np.angle(c)):+.0f}$^\\circ$')
        ax.plot([c.real], [c.imag], marker=T['markers'][k % 8], ms=3.6, color=col,
                mfc='white', mew=0.7, zorder=4)
    ax.set_xlim(-1.25 * R, 1.25 * R)
    ax.set_ylim(-1.25 * R, 1.25 * R)
    ax.set_aspect('equal')
    ax.set_xlabel(p.get('xlabel', 'In-Phase'))
    ax.set_ylabel(p.get('ylabel', 'Quadrature'))
    _viz_panel_title(ax, p, T)
    _viz_axes_style(ax, T)
    _viz_legend(ax, T, n=1, frame=False, loc='upper left', bbox_to_anchor=(1.01, 1.0))


def _panel_table(ax, p, T):
    # booktabs: top rule, header, mid rule, body, bottom rule. No vertical rules.
    ax.set_axis_off()
    _viz_panel_title(ax, p, T, pad=5)
    cols = list(p['columns'])
    rows = [[('--' if c is None else str(c)) for c in r] for r in p['rows']]
    ncol, nrow = len(cols), len(rows)
    w = 1.0 / ncol
    y_top = 0.97
    slot = min(0.20, (y_top - 0.03) / (nrow + 1))
    y_head = y_top - slot
    y_bot = y_head - nrow * slot
    _viz_rules(ax, T, [(y_top, 0.9), (y_head, 0.5), (y_bot, 0.9)])

    def _x(k):
        return (k * w + 0.008, 'left') if k == 0 else ((k + 1) * w - 0.008, 'right')

    for k, c in enumerate(cols):
        x, ha = _x(k)
        heading = textwrap.fill(str(c), max(7, int(72 / ncol)), break_long_words=False)
        ax.text(x, y_top - 0.5 * slot, heading, transform=ax.transAxes, color=T['ink'],
                fontsize=9.5, va='center', ha=ha, linespacing=1.05)
    for i, r in enumerate(rows):
        yc = y_head - (i + 0.5) * slot
        for k, c in enumerate(r):
            x, ha = _x(k)
            ax.text(x, yc, c, transform=ax.transAxes, color=T['ink'], fontsize=10,
                    va='center', ha=ha)
    if p.get('note'):
        ax.text(0.0, y_bot - 0.55 * slot, p['note'], transform=ax.transAxes,
                color=T['muted'], fontsize=10, va='center', ha='left', style='italic')

def _panel_averaging_table(ax, p, T):
    ax.set_axis_off()
    _viz_panel_title(ax, p, T, pad=5)

    # Columns and rows
    cols = ['Metric', 'Value']
    rows = list(p['rows'])  # list of (metric, value)

    ncol = 2
    nrow = len(rows)
    w = 1.0 / ncol

    # Compact manuscript-style vertical layout
    y_top = 0.97
    slot = min(0.16, (y_top - 0.06) / (nrow + 1))
    y_head = y_top - slot
    y_bot = y_head - nrow * slot

    # Booktabs rules: top, header separator, bottom
    _viz_rules(ax, T, [(y_top, 0.9), (y_head, 0.5), (y_bot, 0.9)])

    def _x(k):
        return (k * w + 0.008, 'left') if k == 0 else ((k + 1) * w - 0.008, 'right')

    # Header
    for k, c in enumerate(cols):
        x, ha = _x(k)
        ax.text(
            x, y_top - 0.5 * slot, str(c),
            transform=ax.transAxes,
            color=T['ink'],
            fontsize=16,
            fontfamily='Times New Roman',
            va='center',
            ha=ha
        )

    # Body
    for i, (metric, value, kind) in enumerate(rows):
        yc = y_head - (i + 0.5) * slot

        # Left column
        x0, ha0 = _x(0)
        ax.text(
            x0, yc, metric,
            transform=ax.transAxes,
            color=T['ink'],
            fontsize=16,
            fontfamily='Times New Roman',
            va='center',
            ha=ha0
        )

        # Right column
        x1, ha1 = _x(1)
        if kind == 'bitstring':
            ax.text(
                x1, yc, value,
                transform=ax.transAxes,
                color=T['ink'],
                fontsize=16,  # emphasized, but still manuscript-like
                fontfamily='Times New Roman',
                va='center',
                ha=ha1
            )
        else:
            ax.text(
                x1, yc, value,
                transform=ax.transAxes,
                color=T['ink'],
                fontsize=16,
                fontfamily='Times New Roman',
                va='center',
                ha=ha1
            )

def _panel_text(ax, p, T):
    ax.set_axis_off()
    _viz_panel_title(ax, p, T, pad=5)
    ax.text(0.0, 0.97, p['body'], transform=ax.transAxes, va='top', ha='left',
            fontsize=11, color=T['ink2'], linespacing=1.45)


_VIZ_PANELS = dict(stat=_panel_stat, bars=_panel_bars, stacked=_panel_stacked,
                   dots=_panel_dots, lines=_panel_lines, matrix=_panel_matrix,
                   phasor=_panel_phasor, table=_panel_table, text=_panel_text, averaging_table=_panel_averaging_table)


def plot_summary(title, panels, subtitle='', ncols=2, save_dir=None, filename=None,
                 fig_number=None, letters=True, show=True, row_gap_in=0.0):
    """Render an analysis cell's summary as a publication-style figure.

    `panels` is a list of dicts; every panel takes optional 'span' (columns, default 1),
    'height' (inches) and 'title'. Panel kinds:
      stat     items=[(label, value, note, status), ...]   status flags the value with *
      bars     labels, values | series=[(name, values)], ylabel, fmt, ref, ref_label
      stacked  labels, series=[(name, values)], ylabel      (parts of one total)
      dots     labels, values, xlabel, fmt, logx=True       (decades of range)
      lines    series=[(name, y)] or [(name, x, y)], x, xlabel, ylabel, logx, logy
      matrix   M, labels, fmt, diverging, vmin/vmax/vlim, cbar_label, note
      phasor   values={label: complex}, xlabel, ylabel      (complex baseband)
      table    columns, rows, note
      text     body
    Panels are lettered (a), (b), ... . row_gap_in adds space between rows.
    With save_dir, writes <filename>.png
    and a VIZ_FORMAT vector copy.
    """
    T = _viz_theme()
    rows, cur, used = [], [], 0
    for p in panels:
        span = min(int(p.get('span', 1)), ncols)
        if used + span > ncols and cur:
            rows.append(cur); cur, used = [], 0
        cur.append((p, span, used)); used += span
    if cur:
        rows.append(cur)
    heights = [max(p.get('height', _VIZ_DEFAULT_H.get(p['kind'], 1.8)) for p, _, _ in r)
               for r in rows]
    width = VIZ_WIDTH_IN if ncols > 1 else VIZ_COL_IN
    gap = max(0.0, float(row_gap_in))
    grid_heights = []
    for height in heights:
        if grid_heights and gap:
            grid_heights.append(gap)
        grid_heights.append(height)

    with plt.rc_context(_VIZ_RC):
        fig = plt.figure(figsize=(width, sum(grid_heights) + 0.12), dpi=VIZ_DPI,
                         layout='constrained')
        fig.get_layout_engine().set(h_pad=0.045, w_pad=0.045, hspace=0.03, wspace=0.03)
        gs = fig.add_gridspec(len(grid_heights), ncols, height_ratios=grid_heights)
        k = 0
        for r, row in enumerate(rows):
            for p, span, c0 in row:
                p['_letter'] = f'({chr(97 + k)}) ' if letters else ''
                k += 1
                grid_row = 2 * r if gap else r
                _VIZ_PANELS[p['kind']](fig.add_subplot(gs[grid_row, c0:c0 + span]), p, T)
        if save_dir is not None:
            stem = Path(filename).stem if filename else \
                ''.join(ch if ch.isalnum() else '_' for ch in title.lower()).strip('_')
            for ext in dict.fromkeys([VIZ_FORMAT, 'png']):
                path = Path(save_dir) / f'{stem}.{ext}'
                fig.savefig(path, dpi=VIZ_SAVE_DPI, bbox_inches='tight', pad_inches=0.02)
                print(f"Saved figure to {path.as_posix()}")
        if show:
            plt.show()
    return fig


# --- per-cell wrappers: pass the numbers the cell already computed -------------------

def plot_averaging_table_summary(a_hat, bits, freq_err, avg_samples,
                                 save_dir=None, fig_number=None):
    n_fp = len(a_hat)

    A = {s: np.array([a_hat[g][s] for g in range(n_fp)]) for s in SENSORS}
    cv = np.array([100 * A[s].std() / abs(A[s].mean()) for s in SENSORS])

    allbits = np.array([bits[g] for g in range(n_fp)])
    uniq, counts = np.unique(allbits, axis=0, return_counts=True)

    modal_bits = uniq[counts.argmax()]
    modal = ' '.join(str(int(x)) for x in modal_bits)
    stability = 100 * counts.max() / n_fp

    panels = [
        dict(
            kind='averaging_table',
            span=2,
            height=2.10,
            title='Fingerprint Repeatability Summary',
            rows=[
                ('Most Observed Bitstring', modal, 'bitstring'),
                ('Distinct Bitstrings Observed', f'{len(uniq)}', 'normal'),
                ('Bit Stability', f'{stability:.0f} %', 'normal'),
                ('Worst Single-Channel Drift', f'{cv.max():.3f} %', 'normal'),
            ]
        )
    ]

    return plot_summary(
        'Fingerprint Repeatability Under Averaging',
        panels,
        subtitle=f'{n_fp} reads of {avg_samples} averaged samples '
                 f'({avg_samples/Fs:.3f} s each), $F_s$ = {Fs:.0f} Hz.',
        ncols=2,
        save_dir=save_dir,
        fig_number=fig_number,
        filename='averaging_summary_table',
        letters=False
    )

def plot_cross_sensor_summary(cov_A, corr_A, norm_cov, a_vec, var_A, var_common, var_n,
                              var_m, n_pairs, freq_err, label='', save_dir=None,
                              fig_number=None):
    pct_indep = 100 * np.clip(var_n, 0, None) / var_A
    ch = [s.split('_')[-1] for s in SENSORS]

    panels = [
        dict(kind='bars', labels=ch, values=a_vec,
             ylabel=r'Carrier Amplitude $a_i$ (mV)', fmt='{:.0f}',
             title='Mean Carrier Amplitude per Diode'),
        dict(kind='stacked', labels=ch,
             series=[
                 (r'Shared Common mode, $a_i^2\,\mathrm{Var}(m)$', var_common),
                 (r'Independent Noise, $\mathrm{Var}(n_i)$', np.clip(var_n, 0, None))
             ],
             ylabel=r'Envelope Variance $\mathrm{Var}(A_i)$ (mV$^2$)', fmt='{:.3f}',
             title='Envelope Variance: Shared Common Mode vs Independent Noise'),
        dict(kind='matrix', M=corr_A, diverging=True, vlim=1.0, fmt='{:.3f}',
             title='Envelope Correlation Between Diodes', cbar_label='Correlation'),
        dict(kind='bars', labels=ch, values=pct_indep,
             ylabel='Share of Envelope Variance (%)', fmt='{:.1f}',
             title='Independent Noise Share per Diode'),
    ]

    return plot_summary(
        f'Cross-Sensor Envelope Variance{f", Window {label}" if label else ""}',
        panels,
        save_dir=save_dir,
        fig_number=fig_number,
        filename=f'cross_sensor_{label or "summary"}'
    )


def plot_noise_floor_summary(summary, records, demod_bw, n_samp, save_dir=None,
                             fig_number=None):
    ok = [s for s in SENSORS if summary.get(s)]
    var = np.array([summary[s]['var'] for s in ok])
    snr = np.array([(summary[s].get('snr_db') or np.nan) for s in ok])
    ch = [s.split('_')[-1] for s in ok]

    panels = [
        dict(kind='dots', labels=ch, values=var,
             xlabel=r'Noise variance (mV$^2$)',
             fmt='{:.2e}',
             title=f'Noise Floor ({demod_bw:.0f} Hz)'),
        dict(kind='bars', labels=ch, values=snr,
             ylabel='SNR (dB)',
             fmt='{:.1f}',
             title='Carrier-to-Noise Ratio per Channel'),
        dict(kind='lines', span=2, height=2.35, logy=True,
             series=[
                 (
                     s.split('_')[-1],
                     [r['center'] for r in sorted(records[s], key=lambda r: r['center'])],
                     [r['density'] for r in sorted(records[s], key=lambda r: r['center'])]
                 )
                 for s in ok
             ],
             xlabel='Band Centre Frequency (Hz)',
             ylabel=r'Noise Power Density (mV$^2$/Hz)',
             title='Noise Density vs Frequency'),
    ]

    return plot_summary(
        'Per-Channel Noise Floor and SNR',
        panels,
        save_dir=save_dir,
        fig_number=fig_number,
        filename='noise_floor',
        row_gap_in=0.55
    )


def plot_dark_test_summary(pick, n_samp, halfbw, lit_a=None, save_dir=None,
                           fig_number=None):
    mag_uv = np.array([pick[s]['mag'] * 1000 for s in SENSORS])
    noise_uv = np.array([pick[s]['noise'] * 1000 for s in SENSORS])

    panels = [
        dict(kind='bars', labels=[s.split('_')[-1] for s in SENSORS],
             series=[
                 (r'Coherent Pickup $|c_i|$', mag_uv),
                 (f'In-Band Noise ($\\pm${halfbw:.0f} Hz)', noise_uv)
             ],
             ylabel=r'Amplitude ($\mu$V)', fmt='{:.0f}',
             title='Coherent Pickup vs the Incoherent Noise Floor'),
        dict(kind='phasor',
             values={s: pick[s]['c'] * 1000 for s in SENSORS},
             xlabel=r'In-Phase Component ($\mu$V)',
             ylabel=r'Quadrature Component ($\mu$V)',
             title='Pickup Phasors at the Carrier Frequency'),
    ]

    return plot_summary(
        'Dark Test: Carrier-Band Pickup With the LEDs Off',
        panels,
        save_dir=save_dir,
        fig_number=fig_number,
        filename='dark_test'
    )


def plot_averaging_summary(a_hat, bits, freq_err, avg_samples, save_dir=None,
                           fig_number=None, trend_points=1200):
    n_fp = len(a_hat)
    A = {s: np.array([a_hat[g][s] for g in range(n_fp)]) for s in SENSORS}
    cv = np.array([100 * A[s].std() / abs(A[s].mean()) for s in SENSORS])

    # Reduce the raw per-read jitter to a fixed-size trend over the complete run.
    # Each plotted point is the mean of one non-overlapping read bin; the statistics
    # below still use every unsmoothed read.
    trend_points = max(1, int(trend_points))
    smooth_reads = max(1, int(np.ceil(n_fp / trend_points)))
    starts = np.arange(0, n_fp, smooth_reads)
    stops = np.minimum(starts + smooth_reads, n_fp)
    trend_x = 0.5 * (starts + stops - 1) + 1
    trend = {s: np.array([A[s][a:b].mean() for a, b in zip(starts, stops)])
             for s in SENSORS}
    smooth_seconds = smooth_reads * avg_samples / Fs
    smooth_desc = (f'{smooth_reads:,}-read block mean, {smooth_seconds:.3g} s/bin'
                   if smooth_reads > 1 else '1-read block mean')

    panels = [
        # ppm from each channel's own full-run mean: absolute mV differs substantially
        # between channels, so normalization makes the slow trends comparable.
        dict(kind='lines', span=2, height=2.35, x=trend_x,
             series=[
                 (s.split('_')[-1], 1e6 * (trend[s] / A[s].mean() - 1))
                 for s in SENSORS
             ],
             xlabel='Fingerprint Read',
             ylabel='Deviation From Full-Run Mean (ppm)',
             title=f'Smoothed Per-Read Amplitude Trend'),
        dict(kind='bars', span=2,
             labels=[s.split('_')[-1] for s in SENSORS],
             values=cv,
             ylabel='Coefficient of Variation (%)',
             fmt='{:.3f}',
             title='Amplitude Drift per Sensor'),
    ]

    return plot_summary(
        'Fingerprint Repeatability Under Averaging',
        panels,
        subtitle=f'{n_fp} reads of {avg_samples} averaged samples '
                 f'({avg_samples/Fs:.3f} s each), $F_s$ = {Fs:.0f} Hz. '
                 f'Amplitude traces show {smooth_desc}; summary statistics use all reads.',
        save_dir=save_dir,
        fig_number=fig_number,
        filename='averaging'
    )


# --- added when splitting the notebook ---
# In the notebook the setup cell ran after this cell and set VIZ_DPI = 160, which overrode
# the 150 above. This keeps that same final value.
from qrng_common import VIZ_DPI
