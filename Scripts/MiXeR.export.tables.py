# -*- coding: utf-8 -*-

import os
import re
import glob
import json
import pandas as pd
from pathlib import Path

ROOT = Path(r'E:\BaiduSyncdisk\005.Bioinformatics\SCI\008\influenza_mental_genomics')
TABLES = ROOT / 'Tables'

os.chdir(ROOT)

# ---------------------------------------------------------------------------
# Part 1: export MiXeR results
# ---------------------------------------------------------------------------
def load_gwas_info(excel_path = "GWAS.infos.xlsx"):
    df = pd.read_excel(excel_path)
    df_valid = df[df['group'].notna()].copy()
    lookup = {}
    for _, row in df_valid.iterrows():
        gid = str(row['gwas_id']).strip()
        lookup[gid] = {'group': str(row['group']).strip()}
    return lookup

def find_mixer_results(root, gwas_info = None, group_filter = None):
    results = []
    for d in sorted(glob.glob(os.path.join(root, '*/'))):
        fit = glob.glob(os.path.join(d, '*.fit.json'))
        test = glob.glob(os.path.join(d, '*.test.json'))
        if len(fit) != 1 or len(test) != 1: continue
        key = re.search(r'\(([A-Za-z0-9_\-]+)\)$', os.path.basename(os.path.normpath(d)))
        grp = gwas_info.get(key.group(1), {}).get('group', 'Unknown') if gwas_info else 'Unknown'
        if group_filter and grp != group_filter: continue
        parts = os.path.basename(os.path.normpath(d)).split('_vs_')
        clean = lambda s: re.sub(r'\s*\([^)]*\)', '', s).replace('_', ' ').strip()
        try:
            data_fit = json.loads(open(fit[0], encoding='utf-8').read())
        except Exception:
            print('error reading json in {}, skip'.format(d))
            continue
        results.append({'dir': d, 'trait1': clean(parts[0]), 'trait2': clean(parts[1]), 'group': grp, 'fit': data_fit})
    return results

def export_mixer_table():
    results = find_mixer_results('results/MiXeR', load_gwas_info(), None)
    print('pairs:', len(results))

    rows = []
    for res in results:
        ci = res['fit']['ci']
        ms = res['fit'].get('modelselection', {})
        key = re.search(r'\(([A-Za-z0-9_\-]+)\)$', os.path.basename(os.path.normpath(res['dir']))).group(1)
        row = {'gwas_id': key, 'trait': res['trait2'], 'group': res['group']}
        for k in ['h2_T1', 'h2_T2', 'rg', 'nc1@p9', 'nc2@p9', 'nc1u@p9', 'nc2u@p9', 'nc12@p9', 'dice']:
            row[k + '_mean'] = ci[k]['mean']
            row[k + '_std'] = ci[k]['std']
        row['dAIC_min'] = ms.get('best_vs_min_AIC')
        row['dBIC_min'] = ms.get('best_vs_min_BIC')
        rows.append(row)

    df = pd.DataFrame(rows).sort_values(['group', 'trait']).reset_index(drop = True)
    df = df.loc[:, ['gwas_id', 'trait', 'group', 'h2_T1_mean', 'h2_T1_std', 'h2_T2_mean', 'h2_T2_std', 'rg_mean', 'rg_std',
       'nc1@p9_mean', 'nc1@p9_std', 'nc2@p9_mean', 'nc2@p9_std', 'nc12@p9_mean', 'nc12@p9_std',
       'dice_mean', 'dice_std', 'dAIC_min', 'dBIC_min']]

    
    colnoms = ['GWAS accession', 'Trait', 'Panel', 'Influenza h2 (mean)',
                  'Influenza h2 (SD)', 'Trait h2 (mean)', 'Trait h2 (SD)',
                  'rg (mean)', 'rg (SD)',
                  'Influenza-specific SNPs (mean)', 'Influenza-specific SNPs (SD)',
                  'Trait-specific SNPs (mean)', 'Trait-specific SNPs (SD)',
                  'Shared SNPs (mean)', 'Shared SNPs (SD)',
                  'Dice (mean)', 'Dice (SD)',
                  'Delta AIC', 'Delta BIC']
    
    TABLES.mkdir(parents = True, exist_ok = True)
    out = TABLES / 'Table.S04.mixer_all.xlsx'
    df.to_excel(out, index = False, header = colnoms)
    print(df.shape, '->', out)
    print(df[['gwas_id', 'group', 'rg_mean', 'dice_mean', 'dAIC_min']].head(6).to_string())
    
    
    
    return df


if __name__ == '__main__':
    
    mx = export_mixer_table()