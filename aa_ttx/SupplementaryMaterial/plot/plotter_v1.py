import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import LogLocator, NullFormatter, MaxNLocator
import sys
import random

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------
SORT_MODE = 'global'  # 'global' or 'local'
GLOBAL_SORT_ANCHOR_SQRT_S = 500.0  # anchor slice (None to integrate across range)
SQRT_S_MIN_DISPLAY = 0.1
SQRT_S_MAX_DISPLAY = None
LOWER_INSET_WIDTH_SCALE = 1.0
TOP_PERCENTILES_TO_PLOT = [0.1, 0.2, 0.5, 0.8]
X_MIN, X_MAX = 1.0, 9001.0

# Styling
MAGNIFY = 6
FONT_SIZE_TITLE = 14 + MAGNIFY
FONT_SIZE_SUBTITLE = 10 + MAGNIFY
FONT_SIZE_LEGEND = 10 + MAGNIFY
FONT_SIZE_AXIS = 12 + MAGNIFY
FONT_SIZE_TICKS = 8 + MAGNIFY

# Heatmap options
# None -> label all; list of names -> label those; float threshold -> label >threshold; int N -> top N by strength
HEATMAP_GLOBAL_SORT_Y_LABELS = ['GL262','GL88','GL82','GL46','GL370','GL162','GL402','GL152','GL286','GL160', 'GL288']+['GL274','GL258','GL268','GL94','GL222','GL372', 'GL402', 'GL218']
HEATMAP_ALL_GRAPHS_FONTSIZE = 10
HEATMAP_COLOR_SCALE = 'pow'  # 'linear','log','pow'
HEATMAP_LOG_LINTHRESH = 1.0e-2
HEATMAP_COLOR_POWER_EXPONENT = 0.5
MIN_ROW_WIDTH, MAX_ROW_WIDTH = 1.0, 5.0

# Panel y-limits
# nf=1
#YLIMIT_LINE = (0.1, 300.01)
#YLIMIT_PERCENT = (8.0e-4, 20.01)
#YLIMIT_TOTAL = (1.0e-4, 1.01)

# nf=5
YLIMIT_LINE = (0.1, 5000.01)
YLIMIT_PERCENT = (8.0e-7, 1000.01)
YLIMIT_TOTAL = (1.0e-5, 1.01)

# -----------------------------------------------------------------------------

def load_data(filename):
    raw = eval(open(filename).read())
    entries = [e for e in raw if SQRT_S_MIN_DISPLAY <= e[0] <= (SQRT_S_MAX_DISPLAY or e[0])]
    s_vals = np.array([e[0] for e in entries])
    names = [g[0] for g in sorted(entries[0][1], key=lambda x: x[0])]
    N, M = len(names), len(s_vals)
    data = np.zeros((N, M)); errors = np.zeros((N, M))
    for j, (_, graphs) in enumerate(entries):
        for i, (nm, val_err) in enumerate(sorted(graphs, key=lambda x: x[0])):
            val, err = (val_err if isinstance(val_err, (list, tuple)) else (val_err, 0.0))
            data[i, j], errors[i, j] = val, err
    return s_vals, names, data, errors

def normalize_for_heatmap(data, order):
    subset = data[order, :]
    absmax = np.max(np.abs(subset), axis=0)
    absmax[absmax == 0] = 1.0
    return subset / absmax

def compute_strength(disp):
    if HEATMAP_COLOR_SCALE == 'log':
        norm = mcolors.SymLogNorm(HEATMAP_LOG_LINTHRESH, vmin=-1, vmax=1)
    else:
        norm = mcolors.Normalize(-1, 1)
    mapped = norm(disp)
    mid = norm(0.0)
    return np.mean(np.abs(mapped - mid), axis=1)

def compute_row_edges(strength):
    smin, smax = strength.min(), strength.max()
    norm_s = (strength - smin) / (smax - smin) if smax > smin else np.ones_like(strength)
    widths = MIN_ROW_WIDTH + norm_s * (MAX_ROW_WIDTH - MIN_ROW_WIDTH)
    edges = np.concatenate(([0.0], np.cumsum(widths)))
    mids = (edges[:-1] + edges[1:]) / 2.0
    return edges, mids

def apply_global_sort(data, strength, s_vals):
    if GLOBAL_SORT_ANCHOR_SQRT_S is not None:
        idx0 = np.argmin(np.abs(s_vals - GLOBAL_SORT_ANCHOR_SQRT_S))
        order = np.argsort(np.abs(data[:, idx0]))[::-1]
    else:
        order = np.argsort(strength)[::-1]
    return order

def main(raw_file):
    # Load
    s_vals, names, data, errors = load_data(raw_file)
    M = data.shape[0]

    # Initial order and normalized disp for strength
    default_order = np.arange(M)
    disp0 = normalize_for_heatmap(data, default_order)
    strength0 = compute_strength(disp0)

    # Apply global sort
    order = apply_global_sort(data, strength0, s_vals)
    names_ord = [names[i] for i in order]

    # Normalize for plot and compute strength
    disp_norm = normalize_for_heatmap(data, order)
    strength = strength0[order]

    # Color-scale transform
    if HEATMAP_COLOR_SCALE == 'pow':
        disp_plot = np.sign(disp_norm) * (np.abs(disp_norm) ** HEATMAP_COLOR_POWER_EXPONENT)
        norm = mcolors.Normalize(-1, 1)
    elif HEATMAP_COLOR_SCALE == 'log':
        disp_plot = disp_norm
        norm = mcolors.SymLogNorm(HEATMAP_LOG_LINTHRESH, vmin=-1, vmax=1)
    else:
        disp_plot = disp_norm
        norm = mcolors.Normalize(-1, 1)

    # Compute x-edges on log-scale
    log_s = np.log(s_vals)
    if len(s_vals) > 1:
        mids = (log_s[:-1] + log_s[1:]) / 2.0
        first = log_s[0] - (mids[0] - log_s[0])
        last  = log_s[-1] + (log_s[-1] - mids[-1])
        x_edges = np.exp(np.concatenate(([first], mids, [last])))
    else:
        x_edges = np.exp([log_s[0] - 0.1, log_s[0] + 0.1])

    # Compute y-edges and midpoints
    y_edges, y_mid = compute_row_edges(strength)

    # Panels 2–4 data (unchanged)
    sum_pos = np.sum(np.maximum(0, data), axis=0)
    sum_neg = np.sum(np.minimum(0, data), axis=0)
    total_xs = sum_pos + sum_neg
    rel_pos = np.divide(sum_pos, total_xs, out=np.zeros_like(sum_pos), where=np.abs(total_xs)>1e-15)
    rel_neg = np.divide(sum_neg, total_xs, out=np.zeros_like(sum_neg), where=np.abs(total_xs)>1e-15)
    rel_agg = np.abs(rel_pos) + np.abs(rel_neg) - 1.0

    sum_contrib_top = {p: np.zeros_like(s_vals) for p in TOP_PERCENTILES_TO_PLOT}
    rel_contrib_top = {p: np.zeros_like(s_vals) for p in TOP_PERCENTILES_TO_PLOT}
    for j in range(len(s_vals)):
        col = data[:, j]; tot = total_xs[j]
        idxs = np.argsort(np.abs(col))[::-1]
        for p in TOP_PERCENTILES_TO_PLOT:
            n_top = max(1, int(p * M))
            s_top = np.sum(col[idxs[:n_top]])
            sum_contrib_top[p][j] = s_top
            rel_contrib_top[p][j] = abs(1.0 - s_top / tot) if np.abs(tot)>1e-15 else 0.0

    err_tot = np.sqrt(np.sum(errors**2, axis=0))
    err_pos = np.where(total_xs > 1e-9, err_tot, 0.0)
    err_neg = np.where(total_xs < -1e-9, err_tot, 0.0)

    # Plot
    fig, (ax_hm, ax_line, ax_pct, ax_tot) = plt.subplots(
        4, 1, figsize=(15, 22), sharex=True,
        gridspec_kw={'height_ratios': [3, 1.5, 1.5, 1.5]}
    )

    # Panel 1: Heatmap
    heatmap_title = r'Individual forward scattering graph contributions at $\mu_r = m_Z$ as a function of $\sqrt{s}$'
    im = ax_hm.pcolormesh(x_edges, y_edges, disp_plot,
                          cmap='coolwarm_r', norm=norm, shading='flat')
    ax_hm.set_title(heatmap_title, fontsize=FONT_SIZE_TITLE)
    ax_hm.set_xscale('log')
    ax_hm.set_ylim(y_edges[0], y_edges[-1])
    ax_hm.invert_yaxis()

    if GLOBAL_SORT_ANCHOR_SQRT_S is not None:
        i0 = np.argmin(np.abs(s_vals - GLOBAL_SORT_ANCHOR_SQRT_S))
        ax_hm.axvline(x_edges[i0],   color='k', lw=1)
        ax_hm.axvline(x_edges[i0+1], color='k', lw=1)

    cbar = fig.colorbar(im, ax=ax_hm, orientation='vertical',
                        fraction=0.046, pad=0.02)
    cbar.set_label(
        r'Graph contributions normalized to the max. abs. value per $\sqrt{s}$ slice',
        fontsize=FONT_SIZE_SUBTITLE
    )

    # y-ticks for heatmap
    if HEATMAP_GLOBAL_SORT_Y_LABELS is None:
        ticks, labels = y_mid, names_ord
    elif isinstance(HEATMAP_GLOBAL_SORT_Y_LABELS, int):
        topN = HEATMAP_GLOBAL_SORT_Y_LABELS
        idxs = np.argsort(strength)[::-1][:topN]
        ticks = y_mid[idxs]
        labels = [names_ord[i] for i in idxs]
    elif isinstance(HEATMAP_GLOBAL_SORT_Y_LABELS, float):
        thr = HEATMAP_GLOBAL_SORT_Y_LABELS
        sel = [i for i, s in enumerate(strength) if s >= thr]
        ticks = y_mid[sel]
        labels = [names_ord[i] for i in sel]
    else:
        sel = [names_ord.index(nm) for nm in HEATMAP_GLOBAL_SORT_Y_LABELS if nm in names_ord]
        ticks = y_mid[sel]
        labels = [names_ord[i] for i in sel]

    ax_hm.set_yticks(ticks)
    ax_hm.set_yticklabels(labels, fontsize=HEATMAP_ALL_GRAPHS_FONTSIZE)

    # Panel 2: Relative aggregate contributions
    ax_line.set_yscale('log')
    ax_line.plot(
        s_vals, rel_agg,
        label='1 - sum of relative absolute value contributions',
        color='darkgreen', lw=2
    )
    ax_line.set_ylabel('Relative aggregate contributions', fontsize=FONT_SIZE_AXIS)
    ax_line.set_ylim(YLIMIT_LINE)
    ax_line.grid(True, which='both', linestyle=':', linewidth=0.7)
    ax_line.legend(loc='best', fontsize=FONT_SIZE_TICKS)
    ax_line.set_title('Relative aggregate contributions & cancellations',
                      fontsize=FONT_SIZE_SUBTITLE)
    ax_line.set_xlabel('')

    # Panel 3: Relative contribution from top percentiles
    ax_pct.set_yscale('log')
    colors = ['darkgreen', 'darkorchid', 'sienna', 'teal', 'gold']
    for i, p in enumerate(TOP_PERCENTILES_TO_PLOT):
        ax_pct.plot(
            s_vals, rel_contrib_top[p],
            label=f'abs(1 - sum of top {int(p*100)}% relative contributions)',
            color=colors[i % len(colors)], lw=2
        )
    ax_pct.set_ylabel('Relative aggr. contribution [-]', fontsize=FONT_SIZE_AXIS)
    ax_pct.set_ylim(YLIMIT_PERCENT)
    ax_pct.grid(True, which='both', linestyle=':', linewidth=0.7)
    ax_pct.legend(loc='best', fontsize=FONT_SIZE_TICKS)
    ax_pct.set_title(
        r'Abs. value of relative contribution from percentiles of forward scattering graphs (ranked per $\sqrt{s}$ by magnitude)',
        fontsize=FONT_SIZE_SUBTITLE
    )
    ax_pct.set_xlabel('')

    # Panel 4: Total cross-section (log-log)
    ax_tot.set_xscale('log')
    ax_tot.set_yscale('log')


#    tp = np.ma.masked_where(total_xs <= 1e-9, total_xs)
#    tn = np.ma.masked_where(total_xs >= -1e-9, np.abs(total_xs))
#    if JOIN_LEFT:
#        mi = np.flatnonzero(tp.mask)
#        ui = np.flatnonzero(~tn.mask)
#        if mi.size and ui.size:
#            tp[mi[0]] = tn[ui[0]]
#    else:
#        mi = np.flatnonzero(tn.mask)
#        ui = np.flatnonzero(~tp.mask)
#        if mi.size and ui.size:
#            tn[mi[-1]] = tp[ui[-1]]
#
#    if np.any(tp > 1e-9):
#        ax_tot.plot(s_vals, tp, label='+Δσ', color='darkgreen', lw=2)
#        ax_tot.errorbar(s_vals, tp, yerr=err_pos,
#                        fmt='none', ecolor='gray',
#                        elinewidth=1, capsize=2)
#    if np.any(tn > 1e-9):
#        ax_tot.plot(s_vals, tn, label='-Δσ',
#                    color='darkred', lw=2, ls='--')
#        ax_tot.errorbar(s_vals, tn, yerr=err_neg,
#                        fmt='none', ecolor='gray',
#                        elinewidth=1, capsize=2)



    # mask the two branches
    eps = 1e-12
    tp = np.ma.masked_where(total_xs <= eps, total_xs).copy()
    tn = np.ma.masked_where(total_xs >= -eps, np.abs(total_xs)).copy()

    # 2) bridge every sign‐flip
    for j in range(len(total_xs) - 1):
        a, b = total_xs[j], total_xs[j+1]
        if a * b < 0:  # sign change
            glue = min(abs(a), abs(b))
            if a > 0:
                # pos→neg flip: extend pos-branch at j+1, neg-branch at j
                tp[j+1] = glue
            else:
                # neg→pos flip: extend neg-branch at j+1, pos-branch at j
                tn[j+1] = glue

    # 3) drop any 2-point‐only runs by zeroing them
    def drop_two_point_runs(branch):
        valid = ~branch.mask
        idx = np.where(valid)[0]
        # split into contiguous runs
        groups = np.split(idx, np.where(np.diff(idx) != 1)[0] + 1)
        for run in groups:
            if len(run) == 2:
                branch[run] = 0.0

    drop_two_point_runs(tp)
    drop_two_point_runs(tn)

    # 4) plot both branches with error bars
    if np.any(tp > eps):
        ax_tot.plot(s_vals, tp,   label='+Δσ', color='darkgreen', lw=2)
        ax_tot.errorbar(s_vals, tp, yerr=err_pos,
                        fmt='none', ecolor='gray',
                        elinewidth=1, capsize=2)
    if np.any(tn > eps):
        ax_tot.plot(s_vals, tn,   label='-Δσ', color='darkred',
                    lw=2, ls='--')
        ax_tot.errorbar(s_vals, tn, yerr=err_neg,
                        fmt='none', ecolor='gray',
                        elinewidth=1, capsize=2)


    ax_tot.set_ylabel(r'$\Delta\sigma^{(NNLO~QCD)}$ [pb]',
                      fontsize=FONT_SIZE_AXIS)
    ax_tot.set_xlabel(r'$\sqrt{s}-2 m_t$ [GeV]', fontsize=FONT_SIZE_AXIS)
    ax_tot.set_ylim(YLIMIT_TOTAL)
    ax_tot.grid(True, which='both', linestyle=':', linewidth=0.7)

    # enforce display x-limits
    if X_MIN is not None and X_MAX is not None:
        ax_hm.set_xlim(X_MIN, X_MAX)
        for ax in (ax_line, ax_pct, ax_tot):
            ax.set_xlim(X_MIN, X_MAX)


    # enforce display x-limits
    if X_MIN is not None and X_MAX is not None:
        ax_hm.set_xlim(X_MIN, X_MAX)
        for ax in (ax_line, ax_pct, ax_tot):
            ax.set_xlim(X_MIN, X_MAX)

    # 1) let tight_layout do its thing on margins
    fig.tight_layout(rect=[0, 0.01, 1, 0.96])

    # 2) now manually slide the lower panels under the heatmap data area
    fig.canvas.draw()   # ensure positions are finalized
    heatmap_bbox = ax_hm.get_position()
    h_left  = heatmap_bbox.x0
    h_width = heatmap_bbox.width
    lower_w = h_width * LOWER_INSET_WIDTH_SCALE
    offset  = (h_width * (1.0 - LOWER_INSET_WIDTH_SCALE)) / 2.0
    left    = h_left + offset

    for ax_lower in (ax_line, ax_pct, ax_tot):
        bb = ax_lower.get_position()
        ax_lower.set_position([left, bb.y0, lower_w, bb.height])

    # finally save & show
    fig.savefig('supergraph_hierarchy.pdf')
    plt.show()

if __name__ == '__main__':
    main(sys.argv[1])

