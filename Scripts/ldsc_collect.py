# -*- coding: utf-8 -*-

import os
import re
import pandas as pd

def get_phenotype_name(filepath):
    return os.path.splitext(os.path.basename(filepath))[0]

# ============================================================================
# rg
# ============================================================================
def parse_rg_log(filepath):
    result = {}
    try:
        with open(filepath, 'r', encoding = 'utf-8', errors = 'ignore') as f:
            content = f.read()
    except Exception as e:
        result['parse_error'] = str(e)
        return result

    lines = content.splitlines()

    # --- Phenotype 1 h2 (influenza) ---
    idx1 = content.find('Heritability of phenotype 1')
    idx2 = content.find('Heritability of phenotype 2')
    block1 = content[idx1:idx2] if (idx1 != -1 and idx2 != -1) else ''

    m = re.search(r'Total Observed scale h2:\s*([-\d\.eE]+)\s*\(([\d\.eE]+)\)', block1)
    if m:
        result['h2_obs_p1'] = float(m.group(1))
        result['h2_obs_se_p1'] = float(m.group(2))
    m = re.search(r'Intercept:\s*([-\d\.eE]+)\s*\(([\d\.eE]+)\)', block1)
    if m:
        result['h2_int_p1'] = float(m.group(1))
        result['h2_int_se_p1'] = float(m.group(2))

    # --- Phenotype 2 h2 (target trait) ---
    idx3 = content.find('Heritability of phenotype 2/2')
    idx4 = content.find('Genetic Covariance')
    block2 = content[idx3:idx4] if (idx3 != -1 and idx4 != -1) else ''
    m = re.search(r'Total Observed scale h2:\s*([-\d\.eE]+)\s*\(([\d\.eE]+)\)', block2)
    if m:
        result['h2_obs_p2_block'] = float(m.group(1))
        result['h2_obs_se_p2_block'] = float(m.group(2))

    # --- Summary table ---
    for i, line in enumerate(lines):
        if not line.startswith('Summary of Genetic Correlation Results'):
            continue

        if i + 2 >= len(lines):
            break
        data_line = lines[i + 2].strip()
        if not data_line or data_line.startswith('Analysis'):
            break

        tokens = data_line.split()
        if len(tokens) < 12:
            break

        # Last 10 tokens are numeric
        numeric_part = tokens[-10:]

        keys = ['rg', 'se', 'z', 'p', 'h2_obs', 'h2_obs_se',
                'h2_int', 'h2_int_se', 'gcov_int', 'gcov_int_se']
        for k, v in zip(keys, numeric_part):
            try:
                result[k] = float(v)
            except ValueError:
                result[k] = None if v in ('NA', 'nan', 'NaN') else v
        break

    return result


def collect_rg_results(base_dir):
    records = []
    for root, dirs, files in os.walk(base_dir):
        if os.path.basename(root) == 'h2':
            continue
        disease = os.path.basename(root)
        for fname in files:
            if not fname.endswith('.log'):
                continue
            fpath = os.path.join(root, fname)
            phenotype = get_phenotype_name(fpath)
            parsed = parse_rg_log(fpath)

            record = {'disease': disease, 'gwas_id': phenotype}
            record.update(parsed)
            records.append(record)
    return pd.DataFrame(records)


if __name__ == '__main__':

    BASE_DIR = 'influenza_mental_genomics/results/ldsc'
    OUT_DIR = 'influenza_mental_genomics/Tables'

    eids = ['J10_INFLUENZA', 'GCST90432107']
    
    for eid in eids:
        
        rg_df = collect_rg_results(f'{BASE_DIR}/{eid}')
        
        if not rg_df.empty:
            rg_df = rg_df.sort_values(['disease', 'gwas_id']).reset_index(drop = True)
    
        rg_out = os.path.join(OUT_DIR, f'Table.SS.rg_summary_{eid}.xlsx')
    
        rg_df.to_excel(rg_out, index = False)
    
        print(f'Total rg results: {len(rg_df)}, saved to: {rg_out}')















