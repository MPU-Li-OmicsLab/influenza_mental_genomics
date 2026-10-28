# -*- coding: utf-8 -*-

import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats as st

PROJDIR = '/xbio/articles/influenza_mental_genomics'
COHORT = Path(PROJDIR) / '.xbio' / 'UKB_NHANES' / 'temp' / 'ukb_cohort.parquet'
OLINK = '/xbio/data/UK-BioBank/Protein/olink_instance_0_all.csv'
OUTDIR = Path(PROJDIR) / 'results' / 'UKB_NHANES'
RESULTS_TSV = OUTDIR / 'UKB.protein_flu.results.tsv'
PLOT_TSV = OUTDIR / 'UKB.protein_flu.plot_data.tsv'

FLU_J09J11 = ['p131438', 'p131440', 'p131442']                  # first-reported J09 / J10 / J11
FLU_J09J12 = FLU_J09J11 + ['p131444']                           # + J12 viral pneumonia
FLU_SEASON_MONTHS = [10, 11, 12, 1, 2, 3]                      # Oct-Mar flu season
EXCL_DAYS = 730                                                 # 2-year event exclusion

# (analysis, spec, outcome col, outcome label, CRP, season, subset, exclusion)
SPECS = [
    ('R0_replication',   'Replication of the original MAIN model',
     'incident_j09j11', 'J09-J11', True,  False, None, None),
    ('A1_composite_crp', 'Composite outcome (adds J12 viral pneumonia), +CRP',
     'incident_j09j12', 'J09-J12', True,  False, None, None),
    ('P_primary',        'PRIMARY: composite, -CRP, +season',
     'incident_j09j12', 'J09-J12', False, True,  None, None),
    ('A2_primary_crp',   'Primary + CRP',
     'incident_j09j12', 'J09-J12', True,  True,  None, None),
    ('A4_fluseason',     'Primary, flu-season blood draws only (Oct-Mar)',
     'incident_j09j12', 'J09-J12', False, True,  'fluseason', None),
    ('A5_excl_730d',     'Primary, exclude events <=730 d post-draw',
     'incident_j09j12', 'J09-J12', False, True,  None, EXCL_DAYS),
    ('A6_excl_730d_crp', 'A5 + CRP',
     'incident_j09j12', 'J09-J12', True,  True,  None, EXCL_DAYS),
]

# decomposition diagnostics (caveat-clarifying, not data-mined)
DIAGNOSTICS = [
    ('D1_j12_only',    'J12-only (viral pneumonia)',
     'incident_j12',    'J12',     False, True),
    ('D2_j09j10_only', 'J09/J10-only (identified influenza)',
     'incident_j09j10', 'J09-J10', False, True),
]

# forest display order: primary first, then +/- CRP, then restriction variants
FOREST_ROWS = [                                                 # (analysis, label)
    ('P_primary',        'P  Primary: J09-J12 composite, -CRP, +season'),
    ('A2_primary_crp',   'A2  Primary + CRP'),
    ('A1_composite_crp', 'A1  J09-J12 composite, +CRP (no season)'),
    ('R0_replication',   'R0  Original main replication: J09-J11, +CRP'),
    ('A4_fluseason',     'A4  Primary, flu-season blood draws only'),
    ('A5_excl_730d',     'A5  Primary, events >730 d post-draw only'),
    ('A6_excl_730d_crp', 'A6  A5 + CRP'),
]


def load_data():
    '''Read cohort + ST6GAL1 NPX; derive incident outcomes, season and covariates.'''
    need = (['eid', 'p31', 'p53_i0', 'p54_i0', 'p21000_i0', 'p20116_i0', 'p21003_i0',
             'p21001_i0', 'p30710_i0'] + FLU_J09J12)
    need = list(dict.fromkeys(need))
    df = pd.read_parquet(COHORT, columns = need)
    print('[OK] cohort loaded', df.shape, 'from', COHORT)
    for c in FLU_J09J12 + ['p53_i0']:
        df[c] = pd.to_datetime(df[c], errors = 'coerce')

    # incident = first-reported date strictly after the baseline blood draw
    df['flu_first_j09j11'] = df[FLU_J09J11].min(axis = 1)
    df['flu_first_j09j12'] = df[FLU_J09J12].min(axis = 1)
    df['incident_j09j11'] = (df['flu_first_j09j11'].notna() & (df['flu_first_j09j11'] > df['p53_i0'])).astype(int)
    df['incident_j09j12'] = (df['flu_first_j09j12'].notna() & (df['flu_first_j09j12'] > df['p53_i0'])).astype(int)
    df['incident_j12'] = (df['p131444'].notna() & (df['p131444'] > df['p53_i0'])).astype(int)
    df['incident_j09j10'] = (df[['p131438', 'p131440']].min(axis = 1).notna() & (df[['p131438', 'p131440']].min(axis = 1) > df['p53_i0'])).astype(int)
    df['event_date'] = df['flu_first_j09j12']

    # season of blood draw (index date): Oct-Mar vs Apr-Sep
    df['season_flu'] = df['p53_i0'].dt.month.isin(FLU_SEASON_MONTHS).astype(float)

    print('[STEP] reading ST6GAL1 NPX from', OLINK)
    prot = pd.read_csv(OLINK, usecols = ['eid', 'pst6gal1_i0'])
    n_raw = len(prot)
    prot = prot.drop_duplicates(subset = 'eid')
    prot['st6gal1'] = pd.to_numeric(prot['pst6gal1_i0'], errors = 'coerce')
    if len(prot) < n_raw:
        print('[WARN] dropped %d duplicated eid rows in protein file' % (n_raw - len(prot)))

    d = df.merge(prot[['eid', 'st6gal1']], on = 'eid', how = 'inner')
    print('[OK] merged N=%d | ST6GAL1 NPX n=%d' % (len(d), d['st6gal1'].notna().sum()))
    print('[OK] incident events J09-J11=%d | J09-J12=%d (J12-only col=%d) | J09/J10=%d | flu-season draws=%d (%.1f%%)'
          % (d['incident_j09j11'].sum(), d['incident_j09j12'].sum(),
             (d['p131444'].notna() & (d['p131444'] > d['p53_i0'])).sum(),
             d['incident_j09j10'].sum(), d['season_flu'].sum(),
             100 * d['season_flu'].mean()))

    # covariates: age / sex / BMI / smoking (3-level) / ethnicity flags / centre
    d['sex'] = d['p31'].map({0: 0, 1: 1}).astype(float)
    d['age'] = d['p21003_i0'].astype(float)
    d['bmi'] = d['p21001_i0'].astype(float)
    d.loc[~d['p20116_i0'].isin([0, 1, 2]), 'p20116_i0'] = np.nan
    d['smoking'] = d['p20116_i0'].astype(float)
    eth = d['p21000_i0']
    d['eth_missing'] = np.where(eth.isin([-1, -3]) | eth.isna(), 1.0, 0.0)
    d['eth_nonwhite'] = np.where(eth.notna() & ~eth.isin([-1, -3]), (eth >= 2000).astype(float), 0.0)
    d['logcrp'] = np.log1p(d['p30710_i0'].astype(float))
    d['st6gal1_z'] = (d['st6gal1'] - d['st6gal1'].mean()) / d['st6gal1'].std()
    return d


def covariate_block(dd, use_crp, use_season):
    '''Covariates shared by every spec (optionally +CRP and +season), exposure excluded.'''
    X = pd.DataFrame({'age': dd['age'], 'sex': dd['sex'], 'bmi': dd['bmi'],
                      'smoking': dd['smoking'], 'eth_nonwhite': dd['eth_nonwhite'],
                      'eth_missing': dd['eth_missing']})
    if use_crp:
        X['logcrp'] = dd['logcrp']
    if use_season:
        X['season_flu'] = dd['season_flu']
    return pd.concat([X, pd.get_dummies(dd['p54_i0'], prefix = 'center', drop_first = True).astype(float)], axis = 1)


def base_design(dd, use_crp, use_season):
    '''Covariate block shared by every spec plus the per-SD exposure term.'''
    X = covariate_block(dd, use_crp, use_season)
    X.insert(0, 'st6gal1_z', dd['st6gal1_z'])
    return X


def fit_logit(y, X):
    '''Complete-case logistic fit; bfgs with a Newton polish, Wald SEs from the Hessian.'''
    mod = sm.Logit(y, sm.add_constant(X))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        r1 = mod.fit(disp = 0, maxiter = 1000, method = 'bfgs')
        try:
            res = mod.fit(disp = 0, maxiter = 100, start_params = r1.params, method = 'newton')
            if not res.mle_retvals.get('converged', False):
                raise RuntimeError('newton polish did not converge')
        except Exception:
            res = r1
    return mod, res


def or_per_sd(mod, res):
    '''Exponentiated per-SD estimate with 95% CI and Wald p from the analytic Hessian.'''
    b = float(res.params['st6gal1_z'])
    cov = np.linalg.inv(-mod.hessian(res.params))
    se = float(np.sqrt(cov[list(res.params.index).index('st6gal1_z'), list(res.params.index).index('st6gal1_z')]))
    return b, se, float(np.exp(b)), float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se)), \
        float(2 * st.norm.sf(abs(b / se)))


def fit_spec(d, analysis, spec, outcome, outcome_label, use_crp, use_season, subset, excl_days):
    '''Fit one specification-tree row (OR per +1 SD ST6GAL1 NPX).'''
    dd = d.copy()
    if subset == 'fluseason':
        dd = dd[dd['season_flu'] == 1.0]
    if excl_days is not None:
        dd = dd[~((dd[outcome] == 1) & (dd['event_date'] <= dd['p53_i0'] + pd.Timedelta(days = excl_days)))]
    X = base_design(dd, use_crp, use_season)
    m = pd.concat([dd[[outcome]].astype(float), X], axis = 1).dropna()
    y = m.pop(outcome)
    mod, res = fit_logit(y, m)
    b, se, OR, lo, hi, p = or_per_sd(mod, res)
    restriction = 'none'
    if subset == 'fluseason':
        restriction = 'flu-season draws only'
    elif excl_days is not None:
        restriction = 'exclude events <=%d d' % excl_days
    print('[OK] %-18s N=%6d ev=%4d OR=%.3f (%.3f-%.3f) p=%.4g'
          % (analysis, len(m), y.sum(), OR, lo, hi, p))
    return {'section': 'model', 'analysis': analysis, 'spec': spec, 'term': 'st6gal1_z',
            'contrast': 'per +1 SD NPX', 'estimate_type': 'OR',
            'estimate': round(OR, 4), 'ci_lo': round(lo, 4), 'ci_hi': round(hi, 4), 'p': p,
            'n_total': int(len(m)), 'n_events': int(y.sum()),
            'outcome': outcome_label, 'crp': '+CRP' if use_crp else '-CRP',
            'season': '+season' if use_season else '-', 'restriction': restriction,
            'p_trend': np.nan, 'note': 'beta=%.6f se=%.6f' % (b, se)}


def fit_quintiles(d):
    '''Quintile dose-response on the PRIMARY spec (Q2-Q5 vs Q1) + p-for-trend.'''
    dd = d.copy()
    dd['st6gal1_q'] = pd.qcut(dd['st6gal1'], 5, labels = False)
    X = covariate_block(dd, use_crp = False, use_season = True)
    X = pd.concat([X, pd.get_dummies(dd['st6gal1_q'], prefix = 'q',
                                     drop_first = True).astype(float)], axis = 1)
    m = pd.concat([dd[['incident_j09j12']].astype(float), X], axis = 1).dropna()
    y = m.pop('incident_j09j12')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res_q = sm.Logit(y, sm.add_constant(m)).fit(disp = 0, maxiter = 1000, method = 'bfgs')
    qcols = [c for c in m.columns if c.startswith('q_')]
    cov = np.linalg.inv(-sm.Logit(y, sm.add_constant(m)).hessian(res_q.params))
    names = list(res_q.params.index)

    # p-trend: quintile medians as a single linear term
    med = dd.groupby(pd.qcut(dd['st6gal1'], 5, labels = False))['st6gal1_z'].median()
    Xt = covariate_block(dd, use_crp = False, use_season = True)
    Xt['q_med'] = pd.qcut(dd['st6gal1'], 5, labels = False).map(med)
    mt = pd.concat([dd[['incident_j09j12']].astype(float), Xt], axis = 1).dropna()
    yt = mt.pop('incident_j09j12')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res_t = sm.Logit(yt, sm.add_constant(mt)).fit(disp = 0, maxiter = 1000, method = 'bfgs')
    p_trend = float(2 * st.norm.sf(abs(res_t.params['q_med'] / res_t.bse['q_med'])))

    rows = []
    for qc in qcols:
        b = float(res_q.params[qc])
        se = float(np.sqrt(cov[names.index(qc), names.index(qc)]))
        lvl = int(float(qc.split('_')[1])) + 1
        OR, lo, hi = float(np.exp(b)), float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se))
        p = float(2 * st.norm.sf(abs(b / se)))
        in_q = dd['st6gal1_q'] == lvl - 1
        n_in_q, ev_in_q = int(in_q.sum()), int(dd.loc[in_q, 'incident_j09j12'].sum())
        print('[OK] Q%d (n=%d, events=%d) OR=%.3f (%.3f-%.3f) p=%.4g'
              % (lvl, n_in_q, ev_in_q, OR, lo, hi, p))
        rows.append({'section': 'quintile', 'analysis': 'quintiles_primary',
                     'spec': 'Quintiles of ST6GAL1 NPX, primary spec', 'term': 'Q%d' % lvl,
                     'contrast': 'vs Q1', 'estimate_type': 'OR',
                     'estimate': round(OR, 4), 'ci_lo': round(lo, 4), 'ci_hi': round(hi, 4),
                     'p': p, 'n_total': n_in_q, 'n_events': ev_in_q,
                     'outcome': 'J09-J12', 'crp': '-CRP',
                     'season': '+season', 'restriction': 'none',
                     'p_trend': p_trend, 'note': 'quintile median z=%.4f' % med.loc[lvl - 1]})
    print('[OK] p-trend=%.4g' % p_trend)
    return rows


def fit_rcs(d):
    '''Restricted cubic spline (3 knots) nonlinearity test on the PRIMARY spec.'''
    dd = d.copy()
    knots = np.nanquantile(dd['st6gal1_z'], [0.1, 0.5, 0.9])

    def rcs_basis(x, kn):
        '''Standard 3-knot restricted cubic spline basis (one nonlinear term).'''
        t1, t2, t3 = kn
        x = np.asarray(x, dtype = float)
        return (np.maximum(x - t1, 0) ** 3
                - (t3 - t1) / (t3 - t2) * np.maximum(x - t2, 0) ** 3
                + (t2 - t1) / (t3 - t2) * np.maximum(x - t3, 0) ** 3)

    Xs = covariate_block(dd, use_crp = False, use_season = True)
    Xs['st6_linear'] = dd['st6gal1_z']
    Xs['st6_rcs'] = rcs_basis(dd['st6gal1_z'].values, knots)
    ms = pd.concat([dd[['incident_j09j12']].astype(float), Xs], axis = 1).dropna()
    ys = ms.pop('incident_j09j12')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res_s = sm.Logit(ys, sm.add_constant(ms)).fit(disp = 0, maxiter = 1000, method = 'bfgs')
        res_l = sm.Logit(ys, sm.add_constant(ms.drop(columns = ['st6_rcs']))).fit(
            disp = 0, maxiter = 1000, method = 'bfgs')
    lr_stat = 2 * (res_s.llf - res_l.llf)
    p_nonlin = float(st.chi2.sf(lr_stat, 1))
    b_l, se_l = float(res_l.params['st6_linear']), float(res_l.bse['st6_linear'])
    b_r, se_r = float(res_s.params['st6_rcs']), float(res_s.bse['st6_rcs'])
    note = 'knots z p10/p50/p90 = %.3f / %.3f / %.3f' % (knots[0], knots[1], knots[2])
    print('[OK] RCS knots z = %.3f / %.3f / %.3f | linear-only OR=%.3f p=%.4g | LR p-nonlinearity=%.4g'
          % (knots[0], knots[1], knots[2], np.exp(b_l), 2 * st.norm.sf(abs(b_l / se_l)), p_nonlin))
    return [{'section': 'rcs', 'analysis': 'rcs_primary',
             'spec': 'Restricted cubic spline (3 knots), primary spec', 'term': 'st6gal1_z_linear',
             'contrast': 'per +1 SD NPX', 'estimate_type': 'OR',
             'estimate': round(float(np.exp(b_l)), 4),
             'ci_lo': round(float(np.exp(b_l - 1.96 * se_l)), 4),
             'ci_hi': round(float(np.exp(b_l + 1.96 * se_l)), 4),
             'p': float(2 * st.norm.sf(abs(b_l / se_l))), 'n_total': int(len(ms)),
             'n_events': int(ys.sum()), 'outcome': 'J09-J12', 'crp': '-CRP', 'season': '+season',
             'restriction': 'none', 'p_trend': np.nan, 'note': note + ' (linear-only model)'},
            {'section': 'rcs', 'analysis': 'rcs_primary',
             'spec': 'Restricted cubic spline (3 knots), primary spec', 'term': 'st6_rcs',
             'contrast': 'nonlinear term', 'estimate_type': 'beta',
             'estimate': round(b_r, 6), 'ci_lo': round(b_r - 1.96 * se_r, 6),
             'ci_hi': round(b_r + 1.96 * se_r, 6),
             'p': float(2 * st.norm.sf(abs(b_r / se_r))), 'n_total': int(len(ms)),
             'n_events': int(ys.sum()), 'outcome': 'J09-J12', 'crp': '-CRP', 'season': '+season',
             'restriction': 'none', 'p_trend': np.nan,
             'note': note + ' | spline linear piece beta=%.6f se=%.6f'
             % (float(res_s.params['st6_linear']), float(res_s.bse['st6_linear']))},
            {'section': 'rcs', 'analysis': 'rcs_primary',
             'spec': 'Restricted cubic spline (3 knots), primary spec', 'term': 'nonlinearity_LR',
             'contrast': 'spline vs linear', 'estimate_type': 'LR_chi2',
             'estimate': round(lr_stat, 6), 'ci_lo': np.nan, 'ci_hi': np.nan, 'p': p_nonlin,
             'n_total': int(len(ms)), 'n_events': int(ys.sum()), 'outcome': 'J09-J12',
             'crp': '-CRP', 'season': '+season', 'restriction': 'none', 'p_trend': np.nan,
             'note': note + ' | df=1'}]


def fit_all(d):
    '''Fit every prespecified model; return the combined long result table.'''
    rows = []
    print('[STEP] specification tree (OR per +1 SD ST6GAL1 NPX)')
    for analysis, spec, outcome, outcome_label, use_crp, use_season, subset, excl_days in SPECS:
        rows.append(fit_spec(d, analysis, spec, outcome, outcome_label,
                             use_crp, use_season, subset, excl_days))
    print('[STEP] decomposition diagnostics')
    for analysis, spec, outcome, outcome_label, use_crp, use_season in DIAGNOSTICS:
        rows.append(fit_spec(d, analysis, spec, outcome, outcome_label,
                             use_crp, use_season, None, None))
    print('[STEP] quintile dose-response (primary spec)')
    rows += fit_quintiles(d)
    print('[STEP] restricted cubic spline nonlinearity (primary spec)')
    rows += fit_rcs(d)
    return pd.DataFrame(rows)


def plot_frame(res):
    '''Assemble the long plot table: spec-tree forest rows plus quintile panel rows.'''
    models = res[res['section'] == 'model']
    rows = []
    for analysis, label in FOREST_ROWS:
        r = models[models['analysis'] == analysis].iloc[0]
        rows.append({'panel': 'forest', 'label': label, 'analysis': r['analysis'],
                     'term': r['term'], 'contrast': r['contrast'],
                     'estimate_type': r['estimate_type'], 'estimate': r['estimate'],
                     'ci_lo': r['ci_lo'], 'ci_hi': r['ci_hi'], 'p': r['p'],
                     'n_total': r['n_total'], 'n_events': r['n_events']})
    for _, r in res[res['section'] == 'quintile'].iterrows():
        rows.append({'panel': 'quintile', 'label': r['term'], 'analysis': r['analysis'],
                     'term': r['term'], 'contrast': r['contrast'],
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
            print('[STEP 1/2] fit prespecified ST6GAL1 -> flu models')
            print('=' * 70)
            d = load_data()
            res = fit_all(d)
        if 'tables' in STEPS:
            print('=' * 70)
            print('[STEP 2/2] write results and plot tables')
            print('=' * 70)
            write_tables(res)

    print('=' * 70)
    print('[DONE] analysis finished in %.1f min' % ((time.time() - t0) / 60))
    print('=' * 70)
