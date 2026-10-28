# -*- coding: utf-8 -*-

import time
import os
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

warnings.filterwarnings('ignore', category = FutureWarning)

PROJDIR = '/xbio/articles/influenza_mental_genomics'
EXTRACT = PROJDIR + '/data/KNHANES/KNHANES_2016_2019.extract.tsv'
OUTDIR = PROJDIR + '/results/KNHANES'
RESULTS_TSV = OUTDIR + '/KNHANES.depression_vaccination.results.tsv'
PLOT_TSV = OUTDIR + '/KNHANES.depression_vaccination.plot_data.tsv'

ADULT_AGE = 19
CYCLES = ['2016', '2017', '2018', '2019']
N_CYC = len(CYCLES)
MODELS = ['M1', 'M2']

# M1 = demographics + survey year; M2 = M1 + socioeconomic position
ADJUST_LABELS = {
    'M1': 'age+sex(+cycle)',
    'M2': 'M1+income+education+marital',
}


def load_frame():
    '''Read the KNHANES extract; build vaccination outcome and exposures.'''
    df = pd.read_csv(EXTRACT, sep = '\t', low_memory = False)
    df = df[df['year'].astype(str).isin(CYCLES)].copy()
    df['age'] = pd.to_numeric(df['age'], errors = 'coerce')
    df = df[df['age'] >= ADULT_AGE].copy()
    df['wt'] = pd.to_numeric(df['wt_itvex'], errors = 'coerce') / N_CYC
    df['vax'] = pd.to_numeric(df['BH9_11'], errors = 'coerce')
    df['vax'] = df['vax'].where(df['vax'].isin([1, 2]))
    df['vaccinated'] = (df['vax'] == 1).astype(float).where(df['vax'].notna())
    dep = pd.to_numeric(df['DF2_dg'], errors = 'coerce')
    df['dep_dx'] = dep.where(dep.isin([0, 1]))
    stress = pd.to_numeric(df['mh_stress'], errors = 'coerce')
    df['stress'] = stress.where(stress.isin([0, 1]))
    bs3 = pd.to_numeric(df['BS3_1'], errors = 'coerce')
    df['smoke'] = bs3.map({1: 'current', 2: 'current', 3: 'former', 8: 'never'})
    df['sex'] = pd.to_numeric(df['sex'], errors = 'coerce').map({1: 'M', 2: 'F'})
    df['incm'] = pd.to_numeric(df['incm'], errors = 'coerce').where(lambda s: s.between(1, 4))
    df['edu'] = pd.to_numeric(df['edu'], errors = 'coerce').where(lambda s: s.between(1, 4))
    df['marri_1'] = pd.to_numeric(df['marri_1'], errors = 'coerce').where(lambda s: s.isin([1, 2]))
    df['cycle'] = df['year'].astype(str)
    print('[OK] adults', len(df), '| vax events', int((df['vax'] == 1).sum()),
          '| dep_dx', int(df['dep_dx'].sum()), '| smoke current', int((df['smoke'] == 'current').sum()))
    return df


def fit_weighted_logit(d, formula, term, raw_col):
    '''Weighted logistic regression with HC1 SEs; returns OR row for one term.'''
    d = d.dropna(subset = ['vaccinated', 'wt', raw_col])
    n = len(d)
    events = int(d['vaccinated'].sum())
    exposed = int((d[raw_col].isin([1, 'current', 'former'])).sum())
    if events < 20 or d[raw_col].nunique() < 2:
        return dict(estimate = np.nan, ci_lo = np.nan, ci_hi = np.nan, p = np.nan,
                    n = n, events = events, exposed = exposed, note = 'skipped_low_events')
    try:
        m = smf.glm(formula, data = d, family = sm.families.Binomial(),
                    var_weights = d['wt']).fit(cov_type = 'HC1')
        return dict(estimate = float(np.exp(m.params[term])), ci_lo = float(np.exp(m.conf_int().loc[term, 0])),
                    ci_hi = float(np.exp(m.conf_int().loc[term, 1])), p = float(m.pvalues[term]),
                    n = n, events = events, exposed = exposed, note = '')
    except Exception as e:
        return dict(estimate = np.nan, ci_lo = np.nan, ci_hi = np.nan, p = np.nan,
                    n = n, events = events, exposed = exposed, note = 'fit_failed:%s' % e)


def build_formulas(exposure):
    '''Model formulas for a binary exposure (reference level last).'''
    if exposure == 'smoke':
        return {
            'M1': 'vaccinated ~ C(smoke, Treatment(reference="never")) + age + C(sex) + C(cycle)',
            'M2': ('vaccinated ~ C(smoke, Treatment(reference="never")) + age + C(sex) + C(cycle)'
                   ' + C(incm) + C(edu) + C(marri_1)'),
        }, [('smoke', 'current'), ('smoke', 'former')]
    if exposure == 'dep_dx':
        return {
            'M1': 'vaccinated ~ dep_dx + age + C(sex) + C(cycle)',
            'M2': 'vaccinated ~ dep_dx + age + C(sex) + C(cycle) + C(incm) + C(edu) + C(marri_1)',
        }, [('dep_dx', 1)]
    return {
        'M1': 'vaccinated ~ stress + age + C(sex) + C(cycle)',
        'M2': 'vaccinated ~ stress + age + C(sex) + C(cycle) + C(incm) + C(edu) + C(marri_1)',
    }, [('stress', 1)]


def term_name(exposure, level):
    if exposure == 'smoke':
        return 'C(smoke, Treatment(reference="never"))[T.%s]' % level
    return exposure


def run_all(df):
    '''Fit every exposure x model x cycle combination.'''
    rows = []
    for exposure in ['dep_dx', 'stress', 'smoke']:
        formulas, contrasts = build_formulas(exposure)
        for model in MODELS:
            for cycle in CYCLES + ['Pooled']:
                d = df if cycle == 'Pooled' else df[df['cycle'] == cycle]
                for _, level in contrasts:
                    term = term_name(exposure, level)
                    r = fit_weighted_logit(d, formulas[model], term, exposure)
                    rows.append(dict(exposure = exposure, level = level, model = model, cycle = cycle,
                                     adjustment = ADJUST_LABELS[model], **r))
    return pd.DataFrame(rows)


def plot_frame(res):
    '''Plot table: pooled M2 primary, then per-cycle M2, ordered for the figure.'''
    order = ['dep_dx', 'stress', 'smoke']
    lvl = {(1,): 1}
    rows = []
    o = 0
    for exposure in order:
        sub = res[(res['exposure'] == exposure) & (res['model'] == 'M2')].copy()
        sub['lv'] = sub['level'].map(lambda x: 0 if x in (1, 'current') else (1 if x == 'former' else 0))
        for _, r in sub.sort_values(['cycle', 'lv'], ascending = [True, True]).iterrows():
            o += 1
            rows.append(dict(row_order = o, exposure = r['exposure'], level = r['level'],
                             model = r['model'], cycle = r['cycle'], estimate = r['estimate'],
                             ci_lo = r['ci_lo'], ci_hi = r['ci_hi'], p = r['p'], n = r['n'],
                             events = r['events'], exposed = r['exposed'], note = r['note']))
    return pd.DataFrame(rows)


if __name__ == '__main__':

    t0 = time.time()
    os.makedirs(OUTDIR, exist_ok = True)
    df = load_frame()
    res = run_all(df)
    res = res.sort_values(['exposure', 'model', 'cycle']).reset_index(drop = True)
    res.to_csv(RESULTS_TSV, sep = '\t', index = False)
    print('[OK] written', RESULTS_TSV)
    plot_frame(res).to_csv(PLOT_TSV, sep = '\t', index = False)
    print('[OK] written', PLOT_TSV)

    cols = ['exposure', 'level', 'model', 'cycle', 'estimate', 'ci_lo', 'ci_hi', 'p', 'n', 'events']
    print(res[res['cycle'] == 'Pooled'][cols].to_string(index = False))
    print('[DONE] %.1f s' % (time.time() - t0))
