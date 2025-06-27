import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import LogLocator, NullFormatter, MaxNLocator

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------
SORT_MODE = 'global'  # 'global' or 'local'
GLOBAL_SORT_ANCHOR_SQRT_S = 500.0  # anchor slice (None to integrate across range)
SQRT_S_MIN_DISPLAY = 0.1
SQRT_S_MAX_DISPLAY = None
LOWER_INSET_WIDTH_SCALE = 1.0
TOP_PERCENTILES_TO_PLOT = [0.05, 0.1, 0.2, 0.3]
X_MIN, X_MAX = 0.1, 9001.0
THRESHOLD = 1e-9  # for masking zero/near-zero in total-xs
USE_ABS_CONTRIBS_FOR_PERCENTILE_PLOT = True

# Styling
MAGNIFY = 6
FONT_SIZE_TITLE = 14 + MAGNIFY
FONT_SIZE_SUBTITLE = 10 + MAGNIFY
FONT_SIZE_LEGEND = 10 + MAGNIFY
FONT_SIZE_AXIS = 12 + MAGNIFY
FONT_SIZE_TICKS = 8 + MAGNIFY

MARKER_LINE = 'o'
MARKER_PERCENT = 'o'
MARKER_TOTAL = 'o'
MARKERSIZE_LINE = 3
MARKERSIZE_PERCENT = 3
MARKERSIZE_TOTAL = 3

# layout margins
LEFT_MARGIN      = 0.06   # fraction of figure width (0…1)
RIGHT_MARGIN     = None
VERTICAL_SPACING = 0.1    # relative height between subplots
BOTTOM_MARGIN    = 0.03

# Heatmap options
# None -> label all; list -> those names; float -> threshold; int -> top N
HEATMAP_GLOBAL_SORT_Y_LABELS = [
    'GL262','GL88','GL82','GL46','GL370','GL162','GL402',
    'GL152','GL286','GL160','GL288','GL274','GL258','GL268',
    'GL94','GL222','GL372','GL402','GL218'
]
HEATMAP_ALL_GRAPHS_FONTSIZE = 10
HEATMAP_COLOR_SCALE = 'pow'  # 'linear','log','pow'
HEATMAP_LOG_LINTHRESH = 1.0e-2
HEATMAP_COLOR_POWER_EXPONENT = 0.5
MIN_ROW_WIDTH, MAX_ROW_WIDTH = 1.0, 5.0

# Panel y-limits
YLIMIT_LINE    = (0.9e-1, 2000.01)
if USE_ABS_CONTRIBS_FOR_PERCENTILE_PLOT:
#    YLIMIT_PERCENT = (8.0e-7, 10.01)
    YLIMIT_PERCENT = (1.99, 21.01)
else:
    YLIMIT_PERCENT = (8.0e-7, 2000.01)
YLIMIT_TOTAL   = (1.0e-5, 2.01)
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
    res = subset / absmax[None, :]
    return res 

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
        return np.argsort(np.abs(data[:, idx0]))[::-1]
    else:
        return np.argsort(strength)[::-1]

def drop_two_point_runs(branch):
    valid = ~branch.mask
    try:
        idx = np.where(valid)[0]
    except ValueError:
        return
    groups = np.split(idx, np.where(np.diff(idx) != 1)[0] + 1)
    for run in groups:
        if len(run) == 2:
            #branch[run] = 0.0
            branch.mask[run] = True

#def glue_sign_flips(xs):
#    tp = np.ma.masked_where(xs <= THRESHOLD, xs).copy()
#    tn = np.ma.masked_where(xs >= -THRESHOLD, np.abs(xs)).copy()
#    for j in range(len(xs) - 1):
#        a, b = xs[j], xs[j+1]
#        if a * b < 0:
#            glue = min(abs(a), abs(b))
#            if a > 0:
#                tp[j+1] = glue
#            else:
#                tn[j+1] = glue
#    drop_two_point_runs(tp)
#    drop_two_point_runs(tn)
#    return tp, tn

def glue_sign_flips(xs, errs=None):
    """
    Glue small sign flips in xs and propagate each point’s original error,
    ensuring that whenever we insert a “glue” value we un-mask its error.

    Parameters
    ----------
    xs : array-like, shape (N,)
        Central values.
    errs : array-like or None, shape (N,)
        Corresponding errors for each xs entry.

    Returns
    -------
    tp : MaskedArray shape (N,)
        Positive-branch contributions, with flips glued and two-point runs masked.
    tn : MaskedArray shape (N,)
        Negative-branch contributions (abs), with flips glued and runs masked.
    err_tp : MaskedArray or None
        Errors for tp, masked in the same spots (None if errs is None).
    err_tn : MaskedArray or None
        Errors for tn, masked in the same spots (None if errs is None).
    """
    # initial masks: anything below threshold goes away
    tp = np.ma.masked_where(xs <= THRESHOLD, xs).copy()
    tn = np.ma.masked_where(xs >= -THRESHOLD, np.abs(xs)).copy()

    # if errors were provided, mask them the same way initially
    if errs is not None:
        errs = np.asarray(errs)
        err_tp = np.ma.masked_where(xs <= THRESHOLD, errs).copy()
        err_tn = np.ma.masked_where(xs >= -THRESHOLD, errs).copy()
    else:
        err_tp = err_tn = None

    # glue any genuine sign flips, and un-mask the error at the glued point
    for j in range(len(xs) - 1):
        a, b = xs[j], xs[j+1]
        if a * b < 0:
            glue = min(abs(a), abs(b))
            if a > 0:
                tp[j+1] = glue
                if err_tp is not None:
                    # un-mask error at this point so the bar is drawn
                    err_tp.mask[j+1] = False
            else:
                tn[j+1] = glue
                if err_tn is not None:
                    err_tn.mask[j+1] = False

    # drop isolated two-point runs in both branches (this may re-mask points)
    drop_two_point_runs(tp)
    drop_two_point_runs(tn)

    # make sure the final masks are applied to the error arrays as well
    if err_tp is not None:
        err_tp.mask = tp.mask.copy()
        err_tn.mask = tn.mask.copy()

    return tp, tn, err_tp, err_tn

def main(raw_file):
    # --- NNLO data ---
    s_vals, names, data, errors = load_data(raw_file)
    M = data.shape[0]

    # --- Optional NLO total x-sec ---
    s_nlo = xs_nlo = None
    if os.path.exists('NLO_xsec.py'):
        raw_nlo = sorted(eval(open('NLO_xsec.py').read()), key=lambda x: x[0])
        s_nlo   = np.array([e[0] for e in raw_nlo])
        xs_nlo  = np.array([e[1] for e in raw_nlo])

    # --- Optional LO total x-sec ---
    s_lo = xs_lo = None
    if os.path.exists('LO_xsec.py'):
        raw_lo = sorted(eval(open('LO_xsec.py').read()), key=lambda x: x[0])
        s_lo   = np.array([e[0] for e in raw_lo])
        xs_lo  = np.array([e[1] for e in raw_lo])

    
    # --- Heatmap ordering & normalization ---
    disp0     = normalize_for_heatmap(data, np.arange(M))
    #print([max(abs(float(i)) for i in disp0[:,j]) for j in range(len(disp0[0,:]))])
    strength0 = compute_strength(disp0)
    order     = apply_global_sort(data, strength0, s_vals)
    names_ord = [names[i] for i in order]
    disp_norm = normalize_for_heatmap(data, order)
    #print([max(abs(float(i)) for i in disp_norm[:,j]) for j in range(len(disp_norm[0,:]))])
    strength  = strength0[order]

    # --- Color-scale transform ---
    if HEATMAP_COLOR_SCALE == 'pow':
        disp_plot = np.sign(disp_norm) * (np.abs(disp_norm)**HEATMAP_COLOR_POWER_EXPONENT)
        norm      = mcolors.Normalize(-1, 1)
    elif HEATMAP_COLOR_SCALE == 'log':
        disp_plot = disp_norm
        norm      = mcolors.SymLogNorm(HEATMAP_LOG_LINTHRESH, vmin=-1, vmax=1)
    else:
        disp_plot = disp_norm
        norm      = mcolors.Normalize(-1, 1)

    # --- X edges for heatmap ---
    log_s = np.log(s_vals)
    if len(s_vals) > 1:
        mids = (log_s[:-1] + log_s[1:]) / 2.0
        first= log_s[0] - (mids[0] - log_s[0])
        last = log_s[-1] + (log_s[-1] - mids[-1])
        x_edges = np.exp(np.concatenate(([first], mids, [last])))
    else:
        x_edges = np.exp([log_s[0] - 0.1, log_s[0] + 0.1])

    # --- Y edges for heatmap ---
    y_edges, y_mid = compute_row_edges(strength)

    # --- Aggregate data for lower panels ---
    sum_pos  = np.sum(np.maximum(0, data), axis=0)
    sum_neg  = np.sum(np.minimum(0, data), axis=0)
    total_xs = sum_pos + sum_neg
    rel_pos  = np.divide(sum_pos, total_xs, out=np.zeros_like(sum_pos), where=np.abs(total_xs)>1e-15)
    rel_neg  = np.divide(sum_neg, total_xs, out=np.zeros_like(sum_neg), where=np.abs(total_xs)>1e-15)
    rel_agg  = np.abs(rel_pos) + np.abs(rel_neg) - 1.0

    sum_contrib_top = {p: np.zeros_like(s_vals) for p in TOP_PERCENTILES_TO_PLOT}
    rel_contrib_top = {p: np.zeros_like(s_vals) for p in TOP_PERCENTILES_TO_PLOT}
    for j in range(len(s_vals)):
        col = data[:, j]; 
        if USE_ABS_CONTRIBS_FOR_PERCENTILE_PLOT:
            tot = np.sum(np.abs(col))
        else:
            tot = total_xs[j]
        idxs = np.argsort(np.abs(col))[::-1]
        for p in TOP_PERCENTILES_TO_PLOT:
            n_top = max(1, int(p * M))
            if USE_ABS_CONTRIBS_FOR_PERCENTILE_PLOT:
                s_top = np.sum(np.abs(col[idxs[:n_top]]))
            else:
                s_top = np.sum(col[idxs[:n_top]])
            sum_contrib_top[p][j] = s_top
            #rel_contrib_top[p][j] = abs(1.0 - s_top / (float(n_top/float(M))*tot)) if abs(tot)>1e-15 else 0.0
            rel_contrib_top[p][j] = s_top / (float(n_top/float(M))*tot)

    err_tot = np.sqrt(np.sum(errors**2, axis=0))
    err_pos = np.where(total_xs > THRESHOLD, err_tot, 0.0)
    err_neg = np.where(total_xs < -THRESHOLD, err_tot, 0.0)

    # --- Create figure & axes ---
    fig, (ax_hm, ax_line, ax_pct, ax_tot) = plt.subplots(
        4, 1, figsize=(15, 22), sharex=True,
        gridspec_kw={'height_ratios': [3, 1.5, 1.5, 1.5]}
    )

    # Panel 1: Heatmap
    heatmap_title = (
        r'Individual forward scattering graph contributions at $\mu_r = m_Z$ '
        r'as a function of $\sqrt{s}$'
    )
    #print([max(abs(float(i)) for i in disp_plot[:,j]) for j in range(len(disp_plot[0,:]))])
    #stop
    im = ax_hm.pcolormesh(
        x_edges, y_edges, disp_plot,
        cmap='coolwarm_r', norm=norm, shading='flat'
    )
    ax_hm.set_title(heatmap_title, fontsize=FONT_SIZE_TITLE)
    ax_hm.set_xscale('log')
    ax_hm.set_ylim(y_edges[0], y_edges[-1])
    ax_hm.invert_yaxis()
    if GLOBAL_SORT_ANCHOR_SQRT_S is not None:
        i0 = np.argmin(np.abs(s_vals - GLOBAL_SORT_ANCHOR_SQRT_S))
        ax_hm.axvline(x_edges[i0],   color='k', lw=1)
        ax_hm.axvline(x_edges[i0+1], color='k', lw=1)
    cbar = fig.colorbar(im, ax=ax_hm, orientation='vertical', fraction=0.046, pad=0.02)
    cbar.set_label(
        r'Graph contributions normalized to the max. abs. value per $\sqrt{s}$ slice',
        fontsize=FONT_SIZE_SUBTITLE
    )

    # Align lower insets

    # first let tight_layout handle things...
    fig.tight_layout()

    # then bump out the left margin and add vertical spacing
    fig.subplots_adjust(
        left=LEFT_MARGIN,
        right=RIGHT_MARGIN,
        hspace=VERTICAL_SPACING,
        bottom=BOTTOM_MARGIN
    )
    #fig.tight_layout(rect=[0, 0.01, 1, 0.96])
    
    fig.canvas.draw()
    hb = ax_hm.get_position(); left, width = hb.x0, hb.width
    target_w = width * LOWER_INSET_WIDTH_SCALE
    offset   = (width * (1.0 - LOWER_INSET_WIDTH_SCALE)) / 2.0
    new_left = left + offset
    for ax in (ax_line, ax_pct, ax_tot):
        bb = ax.get_position()
        ax.set_position([new_left, bb.y0, target_w, bb.height])

    # Heatmap y-ticks
    if HEATMAP_GLOBAL_SORT_Y_LABELS is None:
        ticks, labels = y_mid, names_ord
    elif isinstance(HEATMAP_GLOBAL_SORT_Y_LABELS, int):
        topN = HEATMAP_GLOBAL_SORT_Y_LABELS
        idxs = np.argsort(strength)[::-1][:topN]
        ticks = y_mid[idxs]; labels = [names_ord[i] for i in idxs]
    elif isinstance(HEATMAP_GLOBAL_SORT_Y_LABELS, float):
        thr = HEATMAP_GLOBAL_SORT_Y_LABELS
        sel = [i for i, s in enumerate(strength) if s >= thr]
        ticks = y_mid[sel]; labels = [names_ord[i] for i in sel]
    else:
        sel = [names_ord.index(nm) for nm in HEATMAP_GLOBAL_SORT_Y_LABELS if nm in names_ord]
        ticks = y_mid[sel]; labels = [names_ord[i] for i in sel]
    ax_hm.set_yticks(ticks)
    ax_hm.set_yticklabels(labels, fontsize=HEATMAP_ALL_GRAPHS_FONTSIZE)

    # Panel 2: Relative aggregate contributions
    ax_line.set_yscale('log')
    ax_line.plot(
        s_vals, rel_agg,
        label='sum of relative absolute value contributions - 1',
        color='darkgreen', lw=2, marker=MARKER_LINE, markersize=MARKERSIZE_LINE
    )
    ax_line.set_ylabel('Relative aggregate contributions', fontsize=FONT_SIZE_AXIS)
    ax_line.set_ylim(YLIMIT_LINE)
    ax_line.grid(True, which='both', linestyle=':', linewidth=0.7)
    ax_line.legend(loc='best', fontsize=FONT_SIZE_TICKS)
    ax_line.set_title('Relative aggregate contributions & cancellations', fontsize=FONT_SIZE_SUBTITLE)
    ax_line.set_xlabel('')

    # Panel 3: Top-percentile relative contributions
    #ax_pct.set_yscale('log')
    palette = ['darkgreen','darkorchid','sienna','teal','gold']
    for i, p in enumerate(TOP_PERCENTILES_TO_PLOT):
        if not USE_ABS_CONTRIBS_FOR_PERCENTILE_PLOT:    
            ax_pct.plot(
                s_vals, rel_contrib_top[p],
                label=f'abs(1 - sum of top {int(p*100)}% relative contributions)',
                color=palette[i % len(palette)], lw=2, marker=MARKER_PERCENT, markersize=MARKERSIZE_PERCENT
            )
        else:
            ax_pct.plot(
                s_vals, rel_contrib_top[p],
                label=f'(sum of top {int(p*100)}% relative abs. contributions)) / {p}',
                color=palette[i % len(palette)], lw=2, marker=MARKER_PERCENT, markersize=MARKERSIZE_PERCENT
            )

    ax_pct.set_ylabel('Relative aggr. contribution [-]', fontsize=FONT_SIZE_AXIS)
    ax_pct.set_ylim(YLIMIT_PERCENT)
    ax_pct.grid(True, which='both', linestyle=':', linewidth=0.7)
    ax_pct.legend(loc='best', fontsize=FONT_SIZE_TICKS)
    if USE_ABS_CONTRIBS_FOR_PERCENTILE_PLOT:
        ax_pct.set_title(
            r'Cum. relative |contribution| from percentiles of forward scattering graphs '
            r'(ranked per $\sqrt{s}$ by magnitude)',
            fontsize=FONT_SIZE_SUBTITLE
        )
    else:
        ax_pct.set_title(
            r'Abs. value of relative contribution from percentiles of forward scattering graphs '
            r'(ranked per $\sqrt{s}$ by magnitude)',
            fontsize=FONT_SIZE_SUBTITLE
        )
    ax_pct.set_xlabel('')

    # Panel 4: Total cross-section
    ax_tot.set_xscale('log'); ax_tot.set_yscale('log')
    ax_tot.set_title(
        r'Inclusive cross-section for $\gamma \gamma \rightarrow t \bar{t}$ $(\mu_r = m_Z)$, with $n_f=5$ and $n_h(m_t)=1$',
        fontsize=FONT_SIZE_SUBTITLE
    )
    # NNLO
    #tp_nnlo, tn_nnlo = glue_sign_flips(total_xs)
    tp_nnlo, tn_nnlo, err_pos, err_neg = glue_sign_flips(total_xs, err_tot)
    #from pprint import pprint
    #pprint(list(zip(s_vals,tp_nnlo,err_pos)))
    #pprint(list(zip(s_vals,tn_nnlo,err_neg)))
    has_pos_nnlo = np.any(tp_nnlo > THRESHOLD)
    has_neg_nnlo = np.any(tn_nnlo > THRESHOLD)
    if has_pos_nnlo:
        lbl = (r'$\Delta\sigma^{(\mathrm{NNLO~QCD})}$' if not has_neg_nnlo
               else r'$+\Delta\sigma^{(\mathrm{NNLO~QCD})}$')
        ax_tot.plot(s_vals, tp_nnlo, label=lbl, color='darkgreen', lw=2, marker=MARKER_TOTAL, markersize=MARKERSIZE_TOTAL)
        ax_tot.errorbar(s_vals, tp_nnlo, yerr=err_pos,
                        fmt='none', ecolor='gray', elinewidth=1, capsize=2)
    if has_neg_nnlo:
        ax_tot.plot(s_vals, tn_nnlo,
                    label=r'$-\Delta\sigma^{(\mathrm{NNLO~QCD})}$',
                    color='darkred', lw=2, ls='--', marker=MARKER_TOTAL, markersize=MARKERSIZE_TOTAL)
        ax_tot.errorbar(s_vals, tn_nnlo, yerr=err_neg,
                        fmt='none', ecolor='gray', elinewidth=1, capsize=2)

    # NLO
    if s_nlo is not None:
        tp_nlo, tn_nlo, _, _ = glue_sign_flips(xs_nlo)
        mask_nlo   = (s_nlo >= X_MIN) & (s_nlo <= X_MAX)
        has_pos_nlo = np.any(tp_nlo[mask_nlo] > THRESHOLD)
        has_neg_nlo = np.any(tn_nlo[mask_nlo] > THRESHOLD)
        if has_pos_nlo:
            lbl = (r'$\Delta\sigma^{(\mathrm{NLO~QCD})}$' if not has_neg_nlo
                   else r'$+\Delta\sigma^{(\mathrm{NLO~QCD})}$')
            ax_tot.plot(s_nlo, tp_nlo, label=lbl, color='blue', lw=2, ls='-.')
        if has_neg_nlo:
            ax_tot.plot(s_nlo, tn_nlo,
                        label=r'$-\Delta\sigma^{(\mathrm{NLO~QCD})}$',
                        color='cyan', lw=2, ls='-.')

    # LO
    if s_lo is not None:
        tp_lo, tn_lo, _, _ = glue_sign_flips(xs_lo)
        mask_lo   = (s_lo >= X_MIN) & (s_lo <= X_MAX)
        has_pos_lo = np.any(tp_lo[mask_lo] > THRESHOLD)
        has_neg_lo = np.any(tn_lo[mask_lo] > THRESHOLD)
        if has_pos_lo:
            lbl = (r'$\Delta\sigma^{(\mathrm{LO~QCD})}$' if not has_neg_lo
                   else r'$+\Delta\sigma^{(\mathrm{LO~QCD})}$')
            ax_tot.plot(s_lo, tp_lo, label=lbl, color='black', lw=2, ls=':')
        if has_neg_lo:
            ax_tot.plot(s_lo, tn_lo,
                        label=r'$-\Delta\sigma^{(\mathrm{LO~QCD})}$',
                        color='gray', lw=2, ls=':')

    # Total = NLO + LO
    if (s_nlo is not None) and (s_lo is not None):
        if np.array_equal(s_nlo, s_lo):
            s_tot  = s_nlo
            xs_tot = xs_nlo + xs_lo
        else:
            s_tot   = s_nlo
            xs_lo_i = np.interp(s_nlo, s_lo, xs_lo)
            xs_tot  = xs_nlo + xs_lo_i
        tp_tot, tn_tot, _, _ = glue_sign_flips(xs_tot)
        mask_tot       = (s_tot >= X_MIN) & (s_tot <= X_MAX)
        has_pos_tot    = np.any(tp_tot[mask_tot] > THRESHOLD)
        has_neg_tot    = np.any(tn_tot[mask_tot] > THRESHOLD)
        if has_pos_tot:
            ax_tot.plot(s_tot, tp_tot, label='Total', color='k', lw=2)
        if has_neg_tot:
            ax_tot.plot(s_tot, tn_tot, label='-Total', color='k', lw=2, ls='--')

    ax_tot.set_ylabel(r'$\Delta\sigma^{(\mathrm{NNLO~QCD})}$ [pb]', fontsize=FONT_SIZE_AXIS)
    ax_tot.set_xlabel(r'$\sqrt{s}-2 m_t$ [GeV]', fontsize=FONT_SIZE_AXIS)
    ax_tot.set_ylim(YLIMIT_TOTAL)
    ax_tot.grid(True, which='both', linestyle=':', linewidth=0.7)
    ax_tot.legend(fontsize=FONT_SIZE_TICKS, loc='best')

    # Enforce display x-limits
    if X_MIN is not None and X_MAX is not None:
        ax_hm.set_xlim(X_MIN, X_MAX)
        for ax in (ax_line, ax_pct, ax_tot):
            ax.set_xlim(X_MIN, X_MAX)

    fig.savefig('supergraph_hierarchy.pdf')
    plt.show()

if __name__ == '__main__':
    main(sys.argv[1])

