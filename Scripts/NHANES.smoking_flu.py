# -*- coding: utf-8 -*-

import os
import time
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import ConvergenceWarning

# sparse per-cycle event counts: some fits do not converge
warnings.filterwarnings('ignore', category = ConvergenceWarning)
warnings.filterwarnings('ignore', category = RuntimeWarning)

PROJDIR = '/xbio/articles/influenza_mental_genomics'
DATADIR = '/xbio/data/NHANES'
COHORT = 'NHANES_2005_2018'
OUTDIR = os.path.join(PROJDIR, 'results', 'UKB_NHANES')
RESULTS_TSV = os.path.join(OUTDIR, 'NHANES.smoking_flu.results.tsv')
PLOT_TSV = os.path.join(OUTDIR, 'NHANES.smoking_flu.plot_data.tsv')

OUTCOME = 'flu'
MODELS = ['M1', 'M2']
SMOKE_LEVELS = ['never', 'former', 'current']
DOSE_LEVELS = ['never', 'cpd_1_9', 'cpd_10_19', 'cpd_20p']
ADJUST_LABELS = {
    'M1': 'age+sex+race(+cycle)',
    'M2': 'M1+poverty_ratio+education',
}

# cycle -> public file stems
CYCLES = {
    '2005-2006': dict(demo = 'demo_d', hsq = 'hsq_d', smq = 'smq_d'),
    '2007-2008': dict(demo = 'demo_e', hsq = 'hsq_e', smq = 'smq_e'),
    '2009-2010': dict(demo = 'demo_f', hsq = 'hsq_f', smq = 'smq_f'),
    '2011-2012': dict(demo = 'demo_g', hsq = 'hsq_g', smq = 'smq_g'),
    '2013-2014': dict(demo = 'demo_h', hsq = 'hsq_h', smq = 'smq_h'),
    '2015-2016': dict(demo = 'demo_i', hsq = 'hsq_i', smq = 'smq_i'),
    '2017-2018': dict(demo = 'demo_j', hsq = 'hsq_j', smq = 'smq_j'),
}
N_CYC = len(CYCLES)
CYCLE_ORDER = list(CYCLES.keys()) + ['Pooled_2005_2018']


def read_tsv(path):
    '''Read one NHANES TSV with lower-cased column names.'''
    df = pd.read_csv(path, sep = '\t', low_memory = False)
    df.columns = [c.strip().lower() for c in df.columns]
    return df


def load_cycle(cyc, spec):
    '''Merge demographics, flu-like illness and smoking of one cycle.'''
    root = os.path.join(DATADIR, cyc)
    demo = read_tsv(os.path.join(root, 'Demographics', spec['demo'] + '.tsv'))
    keep = ['seqn', 'riagendr', 'ridageyr', 'indfmpir', 'dmdeduc2', 'wtmec2yr']
    racevar = 'ridreth3' if 'ridreth3' in demo.columns else 'ridreth1'
    demo = demo[keep + [racevar]].rename(columns = {racevar: 'race'})
    for c in keep[1:]:
        demo[c] = pd.to_numeric(demo[c], errors = 'coerce')

    hsq = read_tsv(os.path.join(root, 'Questionnaire', spec['hsq'] + '.tsv'))
    hsq = hsq[['seqn', 'hsq520']]
    hsq['hsq520'] = pd.to_numeric(hsq['hsq520'], errors = 'coerce')

    smq = read_tsv(os.path.join(root, 'Questionnaire', spec['smq'] + '.tsv'))
    smq = smq[['seqn', 'smq020', 'smq040', 'smd650']]
    for c in ['smq020', 'smq040', 'smd650']:
        smq[c] = pd.to_numeric(smq[c], errors = 'coerce')

    df = demo.merge(hsq, on = 'seqn', how = 'left').merge(smq, on = 'seqn', how = 'left')
    df['cycle'] = cyc
    return df


def clean(df):
    '''Adults >= 18 with coded smoking status, outcome and model covariates.'''
    d = df[df['ridageyr'] >= 18].copy()
    # outcome: 1 = flu-like illness in past 30 d, 0 = none (7/9 -> missing)
    d['flu'] = np.where(d['hsq520'].isin([1, 2]), (d['hsq520'] == 1).astype(float), np.nan)
    d['male'] = (d['riagendr'] == 1).astype(float)
    d.loc[~d['riagendr'].isin([1, 2]), 'male'] = np.nan
    # exposure: never / former / current from SMQ020 + SMQ040
    d['smoke_c'] = np.where(d['smq020'] == 2, 'never',
                    np.where((d['smq020'] == 1) & (d['smq040'].isin([1, 2])), 'current',
                    np.where((d['smq020'] == 1) & (d['smq040'] == 3), 'former', 'missing')))
    # cigarettes per day (SMD650); 777/999 = refused / don't know -> missing
    d['cpd'] = np.where(d['smd650'].isin([777, 999]), np.nan, d['smd650'])
    d['dose_c'] = np.where(d['smoke_c'] == 'never', 'never',
                   np.where((d['smoke_c'] == 'current') & d['cpd'].notna(),
                            pd.cut(d['cpd'], bins = [0, 9, 19, 200],
                                   labels = DOSE_LEVELS[1:]).astype(object), 'missing'))
    # education categorical with its own missing level;
    # poverty ratio median-imputed with a missingness indicator
    d['educ_c'] = pd.cut(d['dmdeduc2'], bins = [0, 2, 3, 4, 5],
                         labels = ['lt_hs', 'hs', 'some_col', 'college']).astype(object)
    d.loc[~d['dmdeduc2'].isin([1, 2, 3, 4, 5]), 'educ_c'] = 'missing'
    d['pir'] = d['indfmpir']
    d['pir_miss'] = d['pir'].isna().astype(float)
    d['pir'] = d['pir'].fillna(d['pir'].median())
    d['wt'] = d['wtmec2yr']
    return d


def design_matrix(d, model, pooled, exposure):
    '''Build the design matrix of one weighted logistic model.'''
    levels = SMOKE_LEVELS if exposure == 'smoke_c' else DOSE_LEVELS
    lv = d[exposure].astype(object)
    expo = pd.concat([(lv == l).astype(float).rename('expo_%s' % l)
                      for l in levels], axis = 1).iloc[:, 1:]
    expo = expo.loc[:, expo.sum(axis = 0) > 0]
    cols = [expo, d['ridageyr'].rename('age'), d['male'].rename('male')]
    race = pd.get_dummies(d['race'].astype('Int64'), prefix = 'race', drop_first = False)
    race = race.astype(float)
    race.columns = ['race_%s' % c.split('_')[-1] for c in race.columns]
    X = pd.concat(cols, axis = 1)
    if pooled:
        cyc = pd.get_dummies(d['cycle'], prefix = 'cyc', drop_first = True).astype(float)
        X = pd.concat([X, race.iloc[:, 1:], cyc], axis = 1)
    else:
        X = pd.concat([X, race.iloc[:, 1:]], axis = 1)
    if model == 'M2':
        edu = pd.get_dummies(d['educ_c'], prefix = 'edu').astype(float).iloc[:, 1:]
        X = pd.concat([X, d['pir'].rename('pir'), d['pir_miss'].rename('pir_miss'),
                       edu], axis = 1)
    X.insert(0, 'const', 1.0)
    return X.astype(float)


def fit_logit(d, model, exposure, pooled, label):
    '''Fit one weighted logistic model; return one row per contrast as dicts.'''
    need = ['flu', 'ridageyr', 'male', 'race', exposure]
    dd = d[d[exposure].astype(object).isin(
        SMOKE_LEVELS if exposure == 'smoke_c' else DOSE_LEVELS)].copy()
    dd = dd.dropna(subset = need).copy()
    if pooled:
        dd['wt'] = dd['wt'] / N_CYC
    X = design_matrix(dd, model, pooled, exposure)
    y = dd['flu'].astype(float)
    w = dd['wt'].astype(float)
    n_flu = int((dd['flu'] == 1).sum())
    res = sm.GLM(y, X, family = sm.families.Binomial(), var_weights = w) \
             .fit(cov_type = 'HC1')
    ref = SMOKE_LEVELS[0] if exposure == 'smoke_c' else DOSE_LEVELS[0]
    levels = SMOKE_LEVELS if exposure == 'smoke_c' else DOSE_LEVELS
    rows = []
    for lv in levels[1:]:
        key = 'expo_%s' % lv
        if key not in res.params.index:
            continue
        b = res.params[key]; se = res.bse[key]; p = res.pvalues[key]
        n_lv = int((dd[exposure] == lv).sum())
        n_lv_flu = int(((dd[exposure] == lv) & (dd['flu'] == 1)).sum())
        rows.append(dict(analysis = 'smoking_status' if exposure == 'smoke_c'
                         else 'smoking_dose', cycle = label, model = model,
                         outcome = OUTCOME, exposure = lv,
                         contrast = '%s_vs_%s' % (lv, ref),
                         n = int(len(dd)), n_exposed = n_lv, n_exposed_flu = n_lv_flu,
                         n_events = n_flu, estimate = np.exp(b),
                         ci_lo = np.exp(b - 1.96 * se), ci_hi = np.exp(b + 1.96 * se),
                         p_value = p, est_type = 'OR'))
    return rows


def run_all(all_df):
    '''Per-cycle and pooled weighted logistic models for status and dose.'''
    rows = []
    for cyc, d in all_df.items():
        d = clean(d)
        for exposure in ['smoke_c', 'dose_c']:
            for model in MODELS:
                try:
                    rows += fit_logit(d, model, exposure, False, cyc)
                except Exception as e:
                    print('[SKIP] %s %s %s: %s' % (cyc, exposure, model, e))
    pooled = pd.concat([clean(d) for d in all_df.values()], ignore_index = True)
    for exposure in ['smoke_c', 'dose_c']:
        for model in MODELS:
            try:
                rows += fit_logit(pooled, model, exposure, True, 'Pooled_2005_2018')
            except Exception as e:
                print('[SKIP] pooled %s %s: %s' % (exposure, model, e))
    return pd.DataFrame(rows)


def build_results(res):
    '''Add adjustment labels to the long plain-column results table.'''
    d = res.copy()
    d['adjustment'] = d['model'].map(ADJUST_LABELS)
    cols = ['analysis', 'cycle', 'model', 'outcome', 'exposure', 'contrast',
            'adjustment', 'n', 'n_exposed', 'n_exposed_flu', 'n_events',
            'estimate', 'ci_lo', 'ci_hi', 'p_value', 'est_type']
    return d[cols]


def build_plot_data(results):
    '''Add panel and cycle-order columns so R can draw the forest figure.'''
    d = results.copy()
    d['panel'] = np.where(d['analysis'] == 'smoking_dose', 'C_dose_response_or',
                  np.where(d['exposure'] == 'current', 'A_current_vs_never_or',
                           'B_former_vs_never_or'))
    d['cycle_order'] = d['cycle'].map({c: i + 1 for i, c in enumerate(CYCLE_ORDER)})
    d['cycle_is_pooled'] = (d['cycle'] == 'Pooled_2005_2018').astype(int)
    cols = ['panel', 'cycle', 'cycle_order', 'cycle_is_pooled', 'analysis', 'model',
            'outcome', 'exposure', 'contrast', 'adjustment', 'estimate', 'ci_lo',
            'ci_hi', 'p_value', 'est_type', 'n', 'n_exposed', 'n_exposed_flu',
            'n_events']
    d = d[cols].sort_values(['panel', 'cycle_order', 'model'], kind = 'stable')
    return d.reset_index(drop = True)


if __name__ == '__main__':

    # ==== preset parameters ====
    FORCE = True                                        # overwrite existing tables
    STEPS = ['load', 'logit', 'write']

    t0 = time.time()
    os.makedirs(OUTDIR, exist_ok = True)

    all_df = {}
    if 'load' in STEPS:
        print('=' * 70)
        print('[STEP 1/3] load NHANES cycles 2005-2018')
        print('=' * 70)
        for cyc, spec in CYCLES.items():
            raw = load_cycle(cyc, spec)
            all_df[cyc] = raw
            print('[OK] %s: raw n=%d, adults>=18 n=%d' %
                  (cyc, len(raw), int((raw['ridageyr'] >= 18).sum())))

    res = pd.DataFrame()
    if 'logit' in STEPS:
        print('=' * 70)
        print('[STEP 2/3] analysis: smoking status / dose -> flu-like illness')
        print('=' * 70)
        res = run_all(all_df)
        pooled_desc = pd.concat([clean(d) for d in all_df.values()], ignore_index = True)
        a = pooled_desc[pooled_desc['smoke_c'].isin(SMOKE_LEVELS)].dropna(subset = ['flu', 'wt'])
        print('[OK] pooled adults with smoking status+outcome+weight: n=%d' % len(a))
        print('[OK] HSQ520=1 unweighted total: %d | weighted %.2f%%' % (
            int((a['flu'] == 1).sum()),
            100 * a.loc[a['flu'] == 1, 'wt'].sum() / a['wt'].sum()))
        print('[OK] smoking distribution (unweighted):')
        print(a['smoke_c'].value_counts().to_string())
        print('[OK] flu prevalence by group (unweighted %%):')
        print((a.groupby('smoke_c')['flu'].mean() * 100).round(2).to_string())

    if 'write' in STEPS:
        print('=' * 70)
        print('[STEP 3/3] write results and plot-data tables')
        print('=' * 70)
        results = build_results(res)
        plot_data = build_plot_data(results)
        if FORCE or not os.path.exists(RESULTS_TSV):
            results.to_csv(RESULTS_TSV, sep = '\t', index = False)
            print('[OK] wrote %s (%d rows)' % (RESULTS_TSV, len(results)))
        else:
            print('[SKIP] %s already exists' % RESULTS_TSV)
        if FORCE or not os.path.exists(PLOT_TSV):
            plot_data.to_csv(PLOT_TSV, sep = '\t', index = False)
            print('[OK] wrote %s (%d rows)' % (PLOT_TSV, len(plot_data)))
        else:
            print('[SKIP] %s already exists' % PLOT_TSV)

    print('=' * 70)
    print('[DONE] %s finished in %.1f min' % (COHORT, (time.time() - t0) / 60))
    print('=' * 70)
