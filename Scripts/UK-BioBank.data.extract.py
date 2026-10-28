# -*- coding: utf-8 -*-

import time
from pathlib import Path

import pandas as pd

PROJDIR = '/xbio/articles/influenza_mental_genomics'
DATADIR = '/xbio/data/UK-BioBank'
SRC = Path(DATADIR) / 'uk'
OUTDIR = Path(PROJDIR) / 'results.temp' / 'UKB_NHANES'
RESDIR = Path(PROJDIR) / 'results' / 'UKB_NHANES'
COHORT = 'ukb_cohort'
PREFIX = 'UK-BioBank.data.extract'
EXPECTED_ROWS = 502131

# requested fields: field_id -> (instance, block); instance None takes the single instance-less column
REQUESTED = {
    31:    (None, 'demographics'),
    53:    ('0',  'demographics'),
    54:    ('0',  'demographics'),
    189:   ('0',  'demographics'),
    21000: ('0',  'demographics'),
    21001: ('0',  'demographics'),
    21003: ('0',  'demographics'),
    1558:  ('0',  'demographics'),
    20116: ('0',  'demographics'),
    20161: ('0',  'demographics'),
    30710: ('0',  'demographics'),
    2090:  ('0',  'demographics'),
    2100:  ('0',  'demographics'),
    41235: (None, 'demographics'),
    40000: ('0',  'death'),
    40001: ('0',  'death'),
}

# ICD10 first-report blocks: code -> (date field, source field)
INFECTION_BLOCKS = {
    'J09': (131438, 131439), 'J10': (131440, 131441), 'J11': (131442, 131443),
    'J12': (131444, 131445), 'J13': (131446, 131447), 'J14': (131448, 131449),
    'J15': (131450, 131451), 'J16': (131452, 131453), 'J17': (131454, 131455),
    'J18': (131456, 131457),
    'J00': (131424, 131425), 'J01': (131426, 131427), 'J02': (131428, 131429),
    'J03': (131430, 131431), 'J04': (131432, 131433), 'J05': (131434, 131435),
    'J06': (131436, 131437),
}

# psychiatric ICD10 first-report date fields: code -> field_id
PSYCH_BLOCKS = {
    'F17': 130868, 'F20': 130874, 'F25': 130884, 'F30': 130890, 'F31': 130892,
    'F32': 130894, 'F33': 130896, 'F34': 130898, 'F38': 130900, 'F39': 130902,
    'F40': 130904, 'F41': 130906, 'F42': 130908, 'F43': 130910,
}

CHUNK_FOLDERS = ['data_11', 'data_21', 'data_22', 'data_31', 'data_41', 'data_51', 'data_61']

def resolve_path(path_in_csv):
    '''field_id.csv paths lack the _De withdrawal-removal suffix of the real files.'''
    p = Path(path_in_csv.strip().removeprefix('./'))
    return SRC / p.with_name(p.stem + '_De' + p.suffix)

def load_field_index():
    '''Stack the per-chunk field_id.csv files into one field index.'''
    parts = []
    for folder in CHUNK_FOLDERS:
        f = SRC / folder / 'field_id.csv'
        if not f.exists():
            print(f'[WARN] missing {f}, skipped')
            continue
        tab = pd.read_csv(f, dtype = str)
        tab['chunk'] = tab['path'].map(resolve_path)
        parts.append(tab)
    return pd.concat(parts, ignore_index = True)

def build_column_plan(index):
    '''Resolve the requested fields to (field_name, chunk, block) plan entries.'''
    plan, missing = [], []

    def add(fid, inst, block):
        sub = index[index['field_id'] == str(fid)]
        if inst is not None:
            sub = sub[sub['field_name'] == f'p{fid}_i{inst}']
        if sub.empty:
            missing.append((fid, inst, block))
            return
        for _, r in sub.iterrows():
            plan.append((r['field_name'], r['chunk'], block))

    for fid, (inst, block) in REQUESTED.items():
        add(fid, inst, block)
    for _, (d_fid, s_fid) in INFECTION_BLOCKS.items():
        add(d_fid, None, 'infection')
        add(s_fid, None, 'infection')
    for _, d_fid in PSYCH_BLOCKS.items():
        add(d_fid, None, 'psychiatry')
    return plan, missing

def build_cohort(plan):
    '''Read each chunk once and join on eid; every chunk holds the same eid set.'''
    by_chunk = {}
    for fname, chunk, _ in plan:
        by_chunk.setdefault(chunk, []).append(fname)
    cohort = None
    for i, (chunk, names) in enumerate(sorted(by_chunk.items()), 1):
        print(f'[STEP] reading {chunk.relative_to(SRC)} ({len(names)} cols) [{i}/{len(by_chunk)}]')
        tab = pd.read_csv(chunk, usecols = ['eid'] + names, low_memory = False).set_index('eid')
        cohort = tab if cohort is None else cohort.join(tab, how = 'left')
        print(f'       rows so far: {len(cohort):,}')
    cohort = cohort.reset_index()
    if len(cohort) != EXPECTED_ROWS:
        raise SystemExit(f'[ERROR] row count {len(cohort):,} != expected {EXPECTED_ROWS:,}')
    key = lambda c: (int(''.join(ch for ch in c.split('_')[0] if ch.isdigit())), c)
    cols = ['eid'] + sorted([c for c in cohort.columns if c != 'eid'], key = key)
    return cohort[cols]

def write_cohort(cohort):
    '''Write the cohort table as parquet, or tsv.gz when pyarrow is missing.'''
    OUTDIR.mkdir(parents = True, exist_ok = True)
    try:
        import pyarrow  # noqa: F401
        out = OUTDIR / f'{COHORT}.parquet'
        cohort.to_parquet(out, index = False)
    except ImportError:
        out = OUTDIR / f'{COHORT}.tsv.gz'
        cohort.to_csv(out, sep = '\t', index = False, compression = 'gzip')
    print(f'[OK] {cohort.shape} -> {out}')
    return out

def write_summary(cohort, plan):
    '''Write the per-column completeness table and the long plot data.'''
    RESDIR.mkdir(parents = True, exist_ok = True)
    block_of = {fname: block for fname, _, block in plan}
    total = len(cohort)
    rows = [{'column': col, 'block': block_of.get(col, 'key'),
             'non_missing': int(cohort[col].notna().sum()), 'missing': int(cohort[col].isna().sum()),
             'dtype': str(cohort[col].dtype)} for col in cohort.columns]
    summ = pd.DataFrame(rows)
    res = RESDIR / f'{PREFIX}.results.tsv'
    summ.to_csv(res, sep = '\t', index = False)
    print(f'[OK] {summ.shape} -> {res}')
    keep = summ[summ['column'] != 'eid'].reset_index(drop = True)
    plot = pd.concat([
        keep[['column', 'block']].assign(metric = 'non_missing', value = keep['non_missing']),
        keep[['column', 'block']].assign(metric = 'missing', value = keep['missing']),
        keep[['column', 'block']].assign(metric = 'pct_non_missing',
                                         value = keep['non_missing'] / total * 100),
    ], ignore_index = True)
    pdata = RESDIR / f'{PREFIX}.plot_data.tsv'
    plot.to_csv(pdata, sep = '\t', index = False)
    print(f'[OK] {plot.shape} -> {pdata}')

def load_cohort():
    '''Read back a previously written cohort table.'''
    for out in (OUTDIR / f'{COHORT}.parquet', OUTDIR / f'{COHORT}.tsv.gz'):
        if out.exists():
            if out.suffix == '.parquet':
                return pd.read_parquet(out)
            return pd.read_csv(out, sep = '\t')
    raise FileNotFoundError(f'no {COHORT} table under {OUTDIR}')

if __name__ == '__main__':

    # ==== preset parameters ====
    FORCE = True
    STEPS = ['cohort_table', 'summary']

    t0 = time.time()

    print('=' * 70)
    print('[STEP 1/3] load field index and build the column plan')
    print('=' * 70)
    index = load_field_index()
    plan, missing = build_column_plan(index)
    for fid, inst, block in missing:
        print(f'[WARN] field {fid} (instance {inst}, {block}) not in any field_id.csv -> column skipped')
    print(f'[OK] column plan: {len(plan)} columns from {len(set(c[1] for c in plan))} chunks')

    cohort = None
    if 'cohort_table' in STEPS:
        print('=' * 70)
        print('[STEP 2/3] merge field chunks -> ukb_cohort.parquet')
        print('=' * 70)
        if not FORCE and (OUTDIR / f'{COHORT}.parquet').exists():
            print(f'[SKIP] {COHORT}.parquet already exists (FORCE = False)')
        else:
            cohort = build_cohort(plan)
            write_cohort(cohort)

    if 'summary' in STEPS:
        print('=' * 70)
        print('[STEP 3/3] per-column completeness summary')
        print('=' * 70)
        if cohort is None:
            cohort = load_cohort()
        write_summary(cohort, plan)

    print('=' * 70)
    print(f'[DONE] extract finished in {(time.time() - t0) / 60:.1f} min')
    print('=' * 70)
