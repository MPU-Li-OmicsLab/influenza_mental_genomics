# -*- coding: utf-8 -*-

"""
Created on Wed Aug  5 15:11:01 2026
This scripts writen by Abao Xing
@author        : Abao Xing
@Email         : albertxn7@gmail.com

   ┏┓　　┏┓
  ┏┛┻━━━━┛┻┓
  ┃　　　　  ┃
  ┃　━　　━　 ┃
  ┃　┳┛　┗┳　 ┃
  ┃　　　　　 ┃
  ┃　　　┻　　┃
  ┃　　　　　 ┃
  ┗━━┓　　　┏━┛
  　　┃　　 ┃ 神兽保佑
  　　┃　　 ┃ 代码无BUG！！！
  　　┃　　 ┗━━━━━┓
  　　┃　　　　　　  ┣┓
 　　┃　　　　　　  ┏┛┃
 　　┗┓┓┏━━━━━┳┓┏━━┛
  　　┃┫┫　   ┃┫┫
  　　┗┻┛　   ┗┻┛

#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =====================================================================
# MiXeR_grid_figures.py
# ---------------------------------------------------------------------
# Scan MiXeR results (Primary / Replication distinguished via
# GWAS.infos.xlsx), sort to match LDSC forest layout, and draw the four
# standard MiXeR panels. All parameters are set below __main__.
# =====================================================================
"""

import os
import re
import glob
import json
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
from matplotlib_venn import venn2, venn2_circles
from scipy.interpolate import interp1d

import warnings
warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------
# GWAS info loader
# ---------------------------------------------------------------------
def load_gwas_info(excel_path="GWAS.infos.xlsx"):
    import pandas as pd
    df = pd.read_excel(excel_path)
    df_valid = df[df['group'].notna()].copy()
    lookup = {}
    for _, row in df_valid.iterrows():
        gid = str(row['gwas_id'])
        lookup[gid] = {'group': row['group']}
        if pd.notna(row['data']):
            stem = str(row['data']).replace('.tsv.gz', '').replace('.tsv', '')
            lookup[stem] = {'group': row['group']}
    return lookup


def extract_gwas_key(dirname):
    m = re.search(r'\(([A-Za-z0-9_\-]+)\)$', dirname)
    return m.group(1) if m else None


def infer_group(dirname, gwas_info):
    if gwas_info is None:
        return 'Unknown'
    key = extract_gwas_key(os.path.basename(os.path.normpath(dirname)))
    if key and key in gwas_info:
        return gwas_info[key]['group']
    key_lower = key.lower() if key else ''
    for k, v in gwas_info.items():
        if key_lower in k.lower() or k.lower() in key_lower:
            return v['group']
    return 'Unknown'


# ---------------------------------------------------------------------
# Sort order matching LDSC forest plot
# ---------------------------------------------------------------------
LDSC_ORDER = [
    'Opiate dependence',
    'Tobacco use disorder',
    'Alcohol dependence',
    'Schizophrenia',
    'Bipolar disorder',
    'PTSD',
    'Anxiety disorder',
    'Major depressive disorder',
    'Sleep disorder',
    'Tourette syndrome',
    'Anorexia nervosa',
    'Obsessive-compulsive disorder',
    'ADHD',
    'Autism spectrum disorder',
    'Parkinson disease',
    'Alzheimers disease',
    'Suicide attempt',
]

TRAIT_ABBR = {
    'Post traumatic stress disorder': 'PTSD',
    'Attention deficit hyperactivity disorder': 'ADHD',
}


def _abbr(name):
    return TRAIT_ABBR.get(name, name)


def sort_results(results):
    order = {name: i for i, name in enumerate(LDSC_ORDER)}
    for full, short in TRAIT_ABBR.items():
        if short in order:
            order[full] = order[short]
    def key_fn(r):
        return order.get(r['trait2'], 999)
    return sorted(results, key=key_fn)


# ---------------------------------------------------------------------
# Panel drawing routines
# ---------------------------------------------------------------------
def make_venn_plot(data, traits, colors = (0, 1), max_size = None, statistic = ('mean', 'std'), font_size = 11):
    cm = plt.cm.get_cmap('tab10')
    scale_factor = 1000
    with_std = (len(statistic) > 1)
    n1 = data['ci']['nc1@p9'][statistic[0]] / scale_factor
    n1_se = data['ci']['nc1@p9'][statistic[1]] / scale_factor if with_std else None
    n2 = data['ci']['nc2@p9'][statistic[0]] / scale_factor
    n2_se = data['ci']['nc2@p9'][statistic[1]] / scale_factor if with_std else None
    n12 = data['ci']['nc12@p9'][statistic[0]] / scale_factor
    n12_se = data['ci']['nc12@p9'][statistic[1]] / scale_factor if with_std else None
    rg = data['ci']['rg'][statistic[0]]

    if max_size is None: max_size = n1 + n2 + n12
    f = lambda x: x if x < 7 else x + 1

    v = venn2(subsets = (n1, n2, n12), normalize_to = (n1 + n2 + n12) / max_size, set_labels = ("", ""))
    v.get_patch_by_id('100').set_color(cm.colors[f(colors[0])])
    v.get_patch_by_id('010').set_color(cm.colors[f(colors[1])])
    v.get_patch_by_id('110').set_color(cm.colors[7])

    if (n1_se is not None) and (n2_se is not None) and (n12_se is not None):
        formatter1 = '{:.2f}\n({:.2f})' if ((n1 + n12 + n2) < 1) else '{:.1f}\n({:.1f})'
    else:
        formatter1 = '{:.2f}' if ((n1 + n12 + n2) < 1) else '{:.1f}'
    formatter = [formatter1, formatter1, formatter1]
    v.get_label_by_id('100').set_text(formatter[0].format(n1, n1_se))
    v.get_label_by_id('010').set_text(formatter[1].format(n2, n2_se))
    v.get_label_by_id('110').set_text(formatter[2].format(n12, n12_se))

    plt.xlim([-0.75, 0.75]), plt.ylim([-0.7, 0.6])

    clr = plt.cm.get_cmap('seismic')((rg + 1) / 2)
    plt.gca().add_patch(patches.Rectangle(((-abs(0.7 * rg) if (rg < 0) else 0), -0.7), abs(0.7 * rg), 0.15, fill = True, clip_on = False, color = clr))
    plt.gca().add_patch(patches.Rectangle((-0.70, -0.7), 1.4, 0.15, fill = False, clip_on = False))
    plt.gca().add_patch(patches.Rectangle((0, -0.7), 0, 0.15, fill = False, clip_on = False, linewidth = 3))
    plt.gca().text(-0.35 if (rg > 0) else 0.35, -0.7 + 0.13 / 2, '$r_g$={:.2f}'.format(rg), fontsize = font_size, horizontalalignment = 'center', verticalalignment = 'center')


def make_strat_qq_plots(data, flip = False, tick_fontsize = 10, label_fontsize = 10):
    cm = plt.cm.get_cmap('tab10')
    for i in np.array(range(0, 4)) + (4 if flip else 0):
        plt.plot(data['qqplot'][i]['data_logpvec'], data['qqplot'][i]['hv_logp'], color = cm.colors[i % 4], linestyle = 'solid')
    plt.plot(data['qqplot'][i]['hv_logp'], data['qqplot'][i]['hv_logp'], 'k--')
    for i in np.array(range(0, 4)) + (4 if flip else 0):
        plt.plot(data['qqplot'][i]['model_logpvec'], data['qqplot'][i]['hv_logp'], color = cm.colors[i % 4], linestyle = 'dashed')
    plt.ylim(0, 7.3); plt.xlim(0, 7.3)
    plt.xticks(np.arange(0, 7.3, step = 1), fontsize = tick_fontsize)
    plt.yticks(np.arange(0, 7.3, step = 1), fontsize = tick_fontsize)
    plt.xlabel('Expected -$log_{10}(p)$', fontsize = label_fontsize)
    plt.ylabel('Observed -$log_{10}(p)$', fontsize = label_fontsize)


def plot_likelihood(data, label_fontsize = 12):
    if 'likelihood' not in data:
        print('json does not contain likelihood curves, skip likelihood plot')
        return
    cm = plt.cm.get_cmap('tab10')
    like_x = np.array([like_x for like_x, like_y in data['likelihood']])
    like_y = np.array([like_y for like_x, like_y in data['likelihood']])
    like_x_min = np.min(like_x)
    like_x_max = np.max(like_x)
    plot_x = np.arange(like_x_min, like_x_max, (like_x_max - like_x_min) / 100)
    plot_y = np.zeros((like_y.shape[0], len(plot_x)))
    for i in range(like_x.shape[0]):
        plot_y[i, :] = interp1d(like_x[i, :], like_y[i, :], bounds_error = False)(plot_x)
        plot_y[i, :] = plot_y[i, :] - np.nanmin(plot_y[i, :])
    plot_y_def = np.sum(np.isfinite(plot_y), 0)
    plot_y = np.nanmean(plot_y, 0)
    plot_y[plot_y_def < 3] = np.nan
    plt.plot(plot_x / 1000, plot_y)
    for like_x, like_y in data['likelihood']:
        plt.plot(np.array(like_x) / 1000, like_y - np.min(like_y), linestyle = 'dotted', color = cm.colors[0], alpha = 0.3)
    plt.xlabel('Shared variant number [k]', fontsize = label_fontsize)
    plt.ylabel('-log(L) + const', fontsize = label_fontsize)


# ---------------------------------------------------------------------
# Result discovery
# ---------------------------------------------------------------------
def parse_traits(dirname):
    parts = dirname.split('_vs_')
    if len(parts) != 2:
        return dirname, ''
    clean = lambda s: re.sub(r'\s*\([^)]*\)', '', s).replace('_', ' ').strip()
    return clean(parts[0]), clean(parts[1])


def find_mixer_results(root, gwas_info = None, group_filter = None):
    results = []
    for d in sorted(glob.glob(os.path.join(root, '*/'))):
        fit  = glob.glob(os.path.join(d, '*.fit.json'))
        test = glob.glob(os.path.join(d, '*.test.json'))
        if len(fit) != 1 or len(test) != 1:
            continue
        group = infer_group(d, gwas_info)
        if group_filter and group != group_filter:
            continue
        trait1, trait2 = parse_traits(os.path.basename(os.path.normpath(d)))
        try:
            data_fit  = json.loads(open(fit[0]).read())
            data_test = json.loads(open(test[0]).read())
        except Exception:
            print('error reading json in {}, skip'.format(d))
            continue
        results.append({
            'dir': d, 'trait1': trait1, 'trait2': trait2, 'group': group,
            'fit': data_fit, 'test': data_test,
        })
    return results


def shared_venn_max(results, statistic):
    totals = []
    for res in results:
        ci = res['fit']['ci']
        total = sum(ci[k][statistic] for k in ['nc1@p9', 'nc2@p9', 'nc12@p9']) / 1000.0
        totals.append(total)
    return max(totals)


def draw_panel(panel, res, venn_max = None, venn_statistic = ('mean', 'std'),
               venn_fontsize = 11, qq_tick_fontsize = 10, qq_label_fontsize = 10,
               likelihood_label_fontsize = 12):
    t1, t2 = res['trait1'], res['trait2']
    if panel == 'venn':
        make_venn_plot(res['fit'], traits = [t1, t2], max_size = venn_max, statistic = list(venn_statistic), font_size = venn_fontsize)
    elif panel == 'qq_t1_given_t2':
        if 'qqplot' not in res['test']:
            plt.axis('off')
            return
        make_strat_qq_plots(res['test'], flip = False, tick_fontsize = qq_tick_fontsize, label_fontsize = qq_label_fontsize)
    elif panel == 'qq_t2_given_t1':
        if 'qqplot' not in res['test']:
            plt.axis('off')
            return
        make_strat_qq_plots(res['test'], flip = True, tick_fontsize = qq_tick_fontsize, label_fontsize = qq_label_fontsize)
    elif panel == 'likelihood':
        plot_likelihood(res['fit'], label_fontsize = likelihood_label_fontsize)
    else:
        raise ValueError('unknown panel: {}'.format(panel))


def panel_title(panel, res):
    t1, t2 = _abbr(res['trait1']), _abbr(res['trait2'])
    if panel == 'qq_t1_given_t2':
        return '{} | {}'.format(t1, t2)
    if panel == 'qq_t2_given_t1':
        return '{} | {}'.format(t2, t1)
    return '{} & {}'.format(t1, t2)


def qq_legend_handles():
    cm = plt.cm.get_cmap('tab10')
    handles = [Line2D([0], [0], color = cm.colors[i], linewidth = 1.5) for i in range(4)]
    handles.append(Line2D([0], [0], color = 'grey', linewidth = 1.5, linestyle = 'dashed'))
    labels = ['All SNPs', '$P\\leq0.1$', '$P\\leq0.01$', '$P\\leq0.001$', 'MiXeR model (dashed)']
    return handles, labels


def add_qq_legend(fig, legend_fontsize = 10):
    handles, labels = qq_legend_handles()
    fig.legend(handles, labels, loc = 'lower right', bbox_to_anchor = (0.996, 0.065), fontsize = legend_fontsize, frameon = False)


def add_qq_legend_axes(ax, legend_fontsize = 10):
    handles, labels = qq_legend_handles()
    ax.legend(handles, labels, loc = 'lower right', fontsize = legend_fontsize, frameon = False)


def save_figure(fig, out_stem, exts, png_dpi):
    for ext in exts:
        plt.savefig(out_stem + '.' + ext, bbox_inches = 'tight', dpi = (png_dpi if ext == 'png' else None))
        print('Generated ' + out_stem + '.' + ext)
    plt.close(fig)


# ---------------------------------------------------------------------
# Grid figure
# ---------------------------------------------------------------------
def make_grid_figure(results, panel, m, n, out_dir, out_prefix, exts,
                     cell_size, aspect_ratio = 1, venn_max = None,
                     venn_statistic = ('mean', 'std'), png_dpi = 300,
                     title_fontsize = 14, cell_title_fontsize = 10,
                     venn_fontsize = 11, qq_tick_fontsize = 10, qq_label_fontsize = 10,
                     likelihood_label_fontsize = 12, legend_fontsize = 10):
    fig = plt.figure(figsize = (cell_size * aspect_ratio * n, cell_size * m))
    same_t1 = len(set(r['trait1'] for r in results)) == 1

    for i, res in enumerate(results):
        plt.subplot(m, n, i + 1)
        try:
            draw_panel(panel, res, venn_max = venn_max, venn_statistic = venn_statistic,
                       venn_fontsize = venn_fontsize, qq_tick_fontsize = qq_tick_fontsize,
                       qq_label_fontsize = qq_label_fontsize,
                       likelihood_label_fontsize = likelihood_label_fontsize)
        except Exception as e:
            print('WARNING: panel {} failed for {}: {}'.format(panel, res['trait2'], e))
        cell_title = _abbr(res['trait2']) if same_t1 else '{} & {}'.format(_abbr(res['trait1']), _abbr(res['trait2']))
        plt.title(cell_title, fontsize = cell_title_fontsize, fontweight = 'bold')

    for j in range(len(results), m * n):
        ax = plt.subplot(m, n, j + 1)
        ax.axis('off')

    t1 = results[0]['trait1']
    t2 = 'trait2' if same_t1 else ''
    if PANEL_SUPTITLE[panel] != None: 
        fig.suptitle(PANEL_SUPTITLE[panel].format(t1 = t1, t2 = t2), fontsize = title_fontsize)
    plt.tight_layout(rect = [0, 0, 1, 0.98], pad = 1.2)
    if panel.startswith('qq_'):
        add_qq_legend(fig, legend_fontsize = legend_fontsize)

    os.makedirs(out_dir, exist_ok = True)
    save_figure(fig, os.path.join(out_dir, '{}_{}'.format(out_prefix, panel)), exts, png_dpi)


# ---------------------------------------------------------------------
# Split figure
# ---------------------------------------------------------------------
def make_single_figures(results, out_dir, exts,
                        venn_statistic = ('mean', 'std'), png_dpi = 300,
                        single_title_fontsize = 11, venn_fontsize = 11,
                        qq_tick_fontsize = 10, qq_label_fontsize = 10,
                        likelihood_label_fontsize = 12, legend_fontsize = 10):
    os.makedirs(out_dir, exist_ok = True)
    same_t1 = len(set(r['trait1'] for r in results)) == 1

    for res in results:
        stem = re.sub(r'\s+', '_', res['trait2']) if same_t1 else re.sub(r'\s+', '_', '{}_vs_{}'.format(res['trait1'], res['trait2']))
        for panel, (cell_size, aspect_ratio) in PANELS.items():
            plt.figure(figsize = (cell_size * aspect_ratio, cell_size))
            try:
                draw_panel(panel, res, venn_statistic = venn_statistic,
                           venn_fontsize = venn_fontsize,
                           qq_tick_fontsize = qq_tick_fontsize,
                           qq_label_fontsize = qq_label_fontsize,
                           likelihood_label_fontsize = likelihood_label_fontsize)
            except Exception as e:
                print('WARNING: panel {} failed for {}: {}'.format(panel, res['trait2'], e))
                plt.close()
                continue
            plt.title(panel_title(panel, res), fontsize = single_title_fontsize)
            if panel.startswith('qq_'):
                add_qq_legend(plt.gcf(), legend_fontsize = legend_fontsize)
            save_figure(plt.gcf(), os.path.join(out_dir, '{}_{}'.format(stem, panel)),
                        exts, png_dpi)


PANELS = {
    'venn': [3.2, 1],
    'qq_t1_given_t2': [3.2, 0.9],
    'qq_t2_given_t1': [3.2, 0.9],
    'likelihood': [3.2, 1],
}

PANEL_SUPTITLE = {
    'venn'           : None,
    'qq_t1_given_t2' : None,
    'qq_t2_given_t1' : None,
    'likelihood'     : None,
}


if __name__ == '__main__':

    # ---- Preset parameters ----------------------------------------------
    mixer_dir   = 'Results/MiXeR'
    out_dir     = 'Figures'
    gwas_info   = 'GWAS.infos.xlsx'
    group       = 'Primary'           # Primary | Replication | None (all)
    panels      = ['venn', 'qq_t1_given_t2', 'qq_t2_given_t1', 'likelihood']
    split       = False               # True: one 1x1 figure per result per panel
    m, n        = 3, 6                # grid layout (m*n >= number of results)
    exts        = ['pdf', 'png']
    png_dpi     = 300
    # out_prefix auto-assembled from group below
    shared_venn_scale = False         # True: all Venn diagrams share one size scale
    venn_statistic = ('mean', 'std')

    # ---- Font settings --------------------------------------------------
    font_family = 'Times New Roman'
    title_fontsize = 14
    cell_title_fontsize = 12
    single_title_fontsize = 12
    label_fontsize = 12
    tick_fontsize = 12
    venn_fontsize = 12
    legend_fontsize = 8
    # ----------------------------------------------------------------------

    plt.rcParams['font.family'] = font_family

    gwas_map = load_gwas_info(gwas_info) if os.path.exists(gwas_info) else None
    results = find_mixer_results(mixer_dir, gwas_map, group)
    if not results:
        raise SystemExit('no MiXeR results found under: ' + mixer_dir)

    results = sort_results(results)
    out_prefix = group if group else MiXeR
    print('found {} MiXeR results under {}'.format(len(results), mixer_dir))
    print('order: ' + ' -> '.join(r['trait2'] for r in results))

    if split:
        make_single_figures(
            results, out_dir, exts,
            venn_statistic = venn_statistic, png_dpi = png_dpi,
            single_title_fontsize = single_title_fontsize,
            venn_fontsize = venn_fontsize,
            qq_tick_fontsize = tick_fontsize,
            qq_label_fontsize = label_fontsize,
            likelihood_label_fontsize = label_fontsize,
            legend_fontsize = legend_fontsize
        )
    else:
        if len(results) == 1:
            m, n = 1, 1
        if m * n < len(results):
            raise SystemExit('grid {}x{} is too small for {} results'.format(m, n, len(results)))
        venn_max = shared_venn_max(results, venn_statistic[0]) if shared_venn_scale else None
        for panel, (cell_size, aspect_ratio) in PANELS.items():
            if panel not in panels:
                continue
            make_grid_figure(
                results, panel, m, n, out_dir, out_prefix, exts,
                cell_size, aspect_ratio,
                venn_max = venn_max, venn_statistic = venn_statistic,
                png_dpi = png_dpi, title_fontsize = title_fontsize,
                cell_title_fontsize = cell_title_fontsize,
                venn_fontsize = venn_fontsize,
                qq_tick_fontsize = tick_fontsize,
                qq_label_fontsize = label_fontsize,
                likelihood_label_fontsize = label_fontsize,
                legend_fontsize = legend_fontsize
            )

    print('Done.')
