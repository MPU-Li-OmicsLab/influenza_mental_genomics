# -*- coding: utf-8 -*-

import time
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

PROJDIR = '/xbio/articles/influenza_mental_genomics'
COHORT = Path(PROJDIR) / '.xbio' / 'UKB_NHANES' / 'temp' / 'ukb_cohort.parquet'
OUTDIR = Path(PROJDIR) / 'results' / 'UKB_NHANES'
RESULTS_TSV = OUTDIR / 'UKB.smoking_flu.results.tsv'
PLOT_TSV = OUTDIR / 'UKB.smoking_flu.plot_data.tsv'

FLU_COLS = ['p131438', 'p131440', 'p131442']                  # first-reported J09 / J10 / J11
DEATH_COL = 'p40000_i0'                                       # date of death
SMOKING_MAP = {0: 'never', 1: 'previous', 2: 'current'}       # p20116_i0 coding
COVARIATES = 'p21003_i0 + p31 + p21001_i0 + C(eth_grp, Treatment("white"))'
PACKYEAR_UNIT = 10.0                                          # OR per 10 pack-years
BASELINE_TERMS = ['p21003_i0', 'p31', 'p21001_i0', 'eth_grp']
SMK_TERM = 'C(smoking, Treatment("never"))[T.current]'
SMK_PREV_TERM = 'C(smoking, Treatment("never"))[T.previous]'

FOREST_ROWS = [                                               # forest display order
    ('A1_main_logit', 'current', 'Main: current vs never (logit)'),
    ('A2a_3level_logit', 'previous', '3-level: previous vs never'),
    ('A2a_3level_logit', 'current', '3-level: current vs never'),
    ('A2b_packyears_logit', 'py10', 'Per 10 pack-years (ever smokers)'),
    ('A2c_cox_ph', 'current', 'Cox PH: current vs never'),
    ('A2d_crp_logit', 'current', 'Main + log2(CRP)'),
]


def eth_grp(x):
    '''Map UKB ethnicity code p21000_i0 to 6 broad groups; missing/negative -> missing.'''
    if pd.isna(x) or x < 0:
        return 'missing'
    if x == 5:
        return 'chinese'
    if x == 6:
        return 'other'
    if 1000 <= x < 2000:
        return 'white'
    if 2000 <= x < 3000:
        return 'mixed'
    if 3000 <= x < 4000:
        return 'asian'
    if 4000 <= x < 5000:
        return 'black'
    return 'other'


def load_cohort(path):
    '''Read the cohort parquet; derive baseline/flu/death dates, outcome and covariates.'''
    df = pd.read_parquet(path)
    print('[OK] cohort loaded', df.shape, 'from', path)
    df['baseline'] = pd.to_datetime(df['p53_i0'], errors = 'coerce')
    flu = pd.concat([pd.to_datetime(df[c], errors = 'coerce') for c in FLU_COLS], axis = 1)
    df['flu_date'] = flu.min(axis = 1)
    df['incident_flu'] = (df['flu_date'].notna() & (df['flu_date'] > df['baseline'])).astype(int)
    df['prevalent_flu'] = df['flu_date'].notna() & (df['flu_date'] <= df['baseline'])
    df['death'] = pd.to_datetime(df[DEATH_COL], errors = 'coerce')
    admin_censor = pd.concat([flu[c] for c in flu.columns] + [df['death']]).max()
    df['smoking'] = df['p20116_i0'].map(SMOKING_MAP)
    df['eth_grp'] = df['p21000_i0'].apply(eth_grp)
    df['log2_crp'] = np.log2(df['p30710_i0'])
    df['py10'] = df['p20161_i0'] / PACKYEAR_UNIT
    print('[OK] baseline', str(df['baseline'].min().date()), '->',
          str(df['baseline'].max().date()), '| admin censor', str(admin_censor.date()))
    print('[OK] incident flu', int(df['incident_flu'].sum()),
          '| prevalent (excluded)', int(df['prevalent_flu'].sum()))
    return df, admin_censor


def analysis_subset(df, groups, extra_terms = None):
    '''Non-prevalent rows of the given smoking groups with complete covariates.'''
    terms = BASELINE_TERMS + list(extra_terms or [])
    return df[df['smoking'].isin(groups) & ~df['prevalent_flu']].dropna(subset = terms).copy()


def glm_row(model, term, analysis, contrast, d, estimate_type):
    '''One long-format result row (exponentiated estimate with 95% CI) from a GLM fit.'''
    i = list(model.params.index).index(term)
    ci = np.exp(model.conf_int().iloc[i])
    return {'section': 'model', 'analysis': analysis, 'term': term, 'contrast': contrast,
            'estimate_type': estimate_type, 'estimate': round(float(np.exp(model.params.iloc[i])), 4),
            'ci_lo': round(float(ci.iloc[0]), 4), 'ci_hi': round(float(ci.iloc[1]), 4),
            'p': float(model.pvalues.iloc[i]), 'n_total': int(model.nobs),
            'n_events': int(d['incident_flu'].sum())}


def cox_row(d1, admin_censor):
    '''Cox PH current vs never; follow-up ends at flu, death or the admin censoring date.'''
    d = d1.copy()
    d['end'] = admin_censor
    d.loc[d['death'].notna(), 'end'] = d['death']
    is_inc = d['incident_flu'] == 1
    d.loc[is_inc, 'end'] = d.loc[is_inc, 'flu_date']
    earlier_death = is_inc & d['death'].notna() & (d['death'] < d['flu_date'])
    d.loc[earlier_death, 'end'] = d.loc[earlier_death, 'death']
    d['t'] = (d['end'] - d['baseline']).dt.days.astype(float)
    d['event'] = d['incident_flu']
    d = d[d['t'] > 0]
    fit = sm.PHReg.from_formula('t ~ C(smoking, Treatment("never")) + ' + COVARIATES,
                                data = d, status = d['event'].values).fit()
    i = list(fit.model.exog_names).index(SMK_TERM)
    return {'section': 'model', 'analysis': 'A2c_cox_ph', 'term': SMK_TERM,
            'contrast': 'current vs never', 'estimate_type': 'HR',
            'estimate': round(float(np.exp(fit.params[i])), 4),
            'ci_lo': round(float(np.exp(fit.params[i] - 1.96 * fit.bse[i])), 4),
            'ci_hi': round(float(np.exp(fit.params[i] + 1.96 * fit.bse[i])), 4),
            'p': float(fit.pvalues[i]), 'n_total': int(len(d)),
            'n_events': int(d['event'].sum())}


def crude_rate_rows(df):
    '''Crude incident flu counts and rates by smoking group (prevalent cases excluded).'''
    rows = []
    for g in ['never', 'previous', 'current']:
        s = df[(df['smoking'] == g) & (~df['prevalent_flu'])]
        rows.append({'section': 'crude_rate', 'analysis': 'crude_rate', 'term': g,
                     'contrast': np.nan, 'estimate_type': 'rate_pct',
                     'estimate': round(float(100 * s['incident_flu'].mean()), 2),
                     'ci_lo': np.nan, 'ci_hi': np.nan, 'p': np.nan,
                     'n_total': int(len(s)), 'n_events': int(s['incident_flu'].sum())})
    return rows


def fit_all(df, admin_censor):
    '''Fit every prespecified model; return the combined long result table.'''
    rows = []
    d1 = analysis_subset(df, ['never', 'current'])
    print('[OK] A1 main logistic | N=%d events=%d' % (len(d1), d1['incident_flu'].sum()))
    m1 = smf.logit('incident_flu ~ C(smoking, Treatment("never")) + ' + COVARIATES,
                   data = d1).fit(disp = False)
    rows.append(glm_row(m1, SMK_TERM, 'A1_main_logit', 'current vs never', d1, 'OR'))

    d2 = analysis_subset(df, ['never', 'previous', 'current'])
    print('[OK] A2a 3-level logistic | N=%d events=%d' % (len(d2), d2['incident_flu'].sum()))
    m2 = smf.logit('incident_flu ~ C(smoking, Treatment("never")) + ' + COVARIATES,
                   data = d2).fit(disp = False)
    rows.append(glm_row(m2, SMK_PREV_TERM, 'A2a_3level_logit', 'previous vs never', d2, 'OR'))
    rows.append(glm_row(m2, SMK_TERM, 'A2a_3level_logit', 'current vs never', d2, 'OR'))

    d3 = analysis_subset(df, ['previous', 'current'], extra_terms = ['py10'])
    print('[OK] A2b pack-years logistic | N=%d events=%d' % (len(d3), d3['incident_flu'].sum()))
    m3 = smf.logit('incident_flu ~ py10 + ' + COVARIATES, data = d3).fit(disp = False)
    rows.append(glm_row(m3, 'py10', 'A2b_packyears_logit', 'per 10 pack-years', d3, 'OR'))

    r4 = cox_row(d1, admin_censor)
    print('[OK] A2c Cox PH | N=%d events=%d' % (r4['n_total'], r4['n_events']))
    rows.append(r4)

    d5 = d1.dropna(subset = ['log2_crp']).copy()
    print('[OK] A2d main + log2(CRP) | N=%d events=%d' % (len(d5), d5['incident_flu'].sum()))
    m5 = smf.logit('incident_flu ~ C(smoking, Treatment("never")) + log2_crp + ' + COVARIATES,
                   data = d5).fit(disp = False)
    rows.append(glm_row(m5, SMK_TERM, 'A2d_crp_logit', 'current vs never', d5, 'OR'))
    return pd.DataFrame(crude_rate_rows(df) + rows)


def plot_frame(res):
    '''Assemble the long plot table: forest rows in display order plus crude rates.'''
    models = res[res['section'] == 'model']
    rows = []
    for analysis, key, label in FOREST_ROWS:
        r = models[(models['analysis'] == analysis)
                   & models['term'].str.contains(key, regex = False)].iloc[0]
        rows.append({'panel': 'forest', 'label': label, 'analysis': r['analysis'],
                     'term': r['term'], 'contrast': r['contrast'],
                     'estimate_type': r['estimate_type'], 'estimate': r['estimate'],
                     'ci_lo': r['ci_lo'], 'ci_hi': r['ci_hi'], 'p': r['p'],
                     'n_total': r['n_total'], 'n_events': r['n_events']})
    for _, r in res[res['section'] == 'crude_rate'].iterrows():
        rows.append({'panel': 'crude_rate', 'label': 'Crude rate: %s' % r['term'],
                     'analysis': r['analysis'], 'term': r['term'], 'contrast': r['contrast'],
                     'estimate_type': r['estimate_type'], 'estimate': r['estimate'],
                     'ci_lo': r['ci_lo'], 'ci_hi': r['ci_hi'], 'p': r['p'],
                     'n_total': r['n_total'], 'n_events': r['n_events']})
    out = pd.DataFrame(rows)
    out.insert(1, 'row_order', range(1, len(out) + 1))
    return out


def write_tables(res):
    '''Write the results table and the plot table as TSV.'''
    res.to_csv(RESULTS_TSV, sep = '\t', index = False)
    print('[OK] written', RESULTS_TSV)
    plot = plot_frame(res)
    plot.to_csv(PLOT_TSV, sep = '\t', index = False)
    print('[OK] written', PLOT_TSV)


if __name__ == '__main__':

    # ==== preset parameters ====
    FORCE = False      # refit models and overwrite existing output tables
    STEPS = ['models', 'tables']

    t0 = time.time()
    OUTDIR.mkdir(parents = True, exist_ok = True)

    if RESULTS_TSV.exists() and PLOT_TSV.exists() and not FORCE:
        print('[SKIP] both output tables exist (set FORCE = True to refit)')
    else:
        if 'models' in STEPS:
            print('=' * 70)
            print('[STEP 1/2] fit prespecified smoking -> flu models')
            print('=' * 70)
            df, admin_censor = load_cohort(COHORT)
            res = fit_all(df, admin_censor)
            print(res[res['section'] == 'model'].to_string(index = False))
        if 'tables' in STEPS:
            print('=' * 70)
            print('[STEP 2/2] write results and plot tables')
            print('=' * 70)
            write_tables(res)

    print('=' * 70)
    print('[DONE] analysis finished in %.1f min' % ((time.time() - t0) / 60))
    print('=' * 70)
