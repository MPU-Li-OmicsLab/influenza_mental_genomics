# -*- coding: utf-8 -*-

import os
import time
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.tools.sm_exceptions import ConvergenceWarning

# sparse per-cycle event counts: some fits do not converge, and the
# partial-likelihood score/hessian can divide by empty risk sets
warnings.filterwarnings('ignore', category = ConvergenceWarning)
warnings.filterwarnings('ignore', category = RuntimeWarning)

PROJDIR = '/xbio/articles/influenza_mental_genomics'
DATADIR = '/xbio/data/NHANES'
COHORT = 'NHANES_2005_2018'
OUTDIR = os.path.join(PROJDIR, 'results', 'UKB_NHANES')
RESULTS_TSV = os.path.join(OUTDIR, 'NHANES.flu_depression.results.tsv')
PLOT_TSV = os.path.join(OUTDIR, 'NHANES.flu_depression.plot_data.tsv')

OUTCOMES = ['phq9', 'phq9_ge10']
MODELS = ['M1', 'M2', 'M3']
EXPOSURES_COX = ['phq9_per5', 'phq9_ge10']

# Cox solver settings; the observed information of the pooled fit is
# near-singular, so the covariance is inverted by eigen decomposition
COX_MAXITER = 500
COX_TOL = 1e-8
COX_RIDGE = 1e-10
COX_EIG_RTOL = 1e-10
ADJUST_LABELS = {
    'M1': 'age+sex+race(+cycle)',
    'M2': 'M1+poverty_ratio+education+smoking',
    'M3': 'M2+log(CRP)',
}

# cycle -> public file stems; crp = (lab stem, crp variable, multiplier to
# mg/L); None = CRP not measured in that cycle
CYCLES = {
    '2005-2006': dict(demo = 'demo_d', hsq = 'hsq_d', dpq = 'dpq_d', smq = 'smq_d', crp = ('crp_d', 'lbxcrp', 10.0)),
    '2007-2008': dict(demo = 'demo_e', hsq = 'hsq_e', dpq = 'dpq_e', smq = 'smq_e', crp = ('crp_e', 'lbxcrp', 10.0)),
    '2009-2010': dict(demo = 'demo_f', hsq = 'hsq_f', dpq = 'dpq_f', smq = 'smq_f', crp = ('crp_f', 'lbxcrp', 10.0)),
    '2011-2012': dict(demo = 'demo_g', hsq = 'hsq_g', dpq = 'dpq_g', smq = 'smq_g', crp = None),
    '2013-2014': dict(demo = 'demo_h', hsq = 'hsq_h', dpq = 'dpq_h', smq = 'smq_h', crp = None),
    '2015-2016': dict(demo = 'demo_i', hsq = 'hsq_i', dpq = 'dpq_i', smq = 'smq_i', crp = ('hscrp_i', 'lbxhscrp', 1.0)),
    '2017-2018': dict(demo = 'demo_j', hsq = 'hsq_j', dpq = 'dpq_j', smq = 'smq_j', crp = ('hscrp_j', 'lbxhscrp', 1.0)),
}
N_CYC = len(CYCLES)
PHQ_ITEMS = ['dpq%03d' % i for i in range(10, 100, 10)]
CYCLE_ORDER = list(CYCLES.keys()) + ['Pooled_2005_2018']


def read_tsv(path):
    '''Read one NHANES TSV with lower-cased column names.'''
    df = pd.read_csv(path, sep = '\t', low_memory = False)
    df.columns = [c.strip().lower() for c in df.columns]
    return df


def load_cycle(cyc, spec):
    '''Merge demographics, flu-like illness, PHQ-9, smoking and CRP of one cycle.'''
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

    dpq = read_tsv(os.path.join(root, 'Questionnaire', spec['dpq'] + '.tsv'))
    items = [i for i in PHQ_ITEMS if i in dpq.columns]
    dpq = dpq[['seqn'] + items]
    for i in items:
        dpq[i] = pd.to_numeric(dpq[i], errors = 'coerce')

    smq = read_tsv(os.path.join(root, 'Questionnaire', spec['smq'] + '.tsv'))
    smq = smq[['seqn', 'smq020', 'smq040']]
    for c in ['smq020', 'smq040']:
        smq[c] = pd.to_numeric(smq[c], errors = 'coerce')

    df = demo.merge(hsq, on = 'seqn', how = 'left').merge(dpq, on = 'seqn', how = 'left') \
             .merge(smq, on = 'seqn', how = 'left')

    if spec['crp'] is not None:
        stem, var, mult = spec['crp']
        crp = read_tsv(os.path.join(root, 'Laboratory', stem + '.tsv'))
        crp = crp[['seqn', var]].rename(columns = {var: 'crp_raw'})
        crp['crp_raw'] = pd.to_numeric(crp['crp_raw'], errors = 'coerce')
        crp['crp_mgl'] = crp['crp_raw'] * mult
    else:
        crp = pd.DataFrame({'seqn': df['seqn'], 'crp_mgl': np.nan})
    df = df.merge(crp[['seqn', 'crp_mgl']], on = 'seqn', how = 'left')
    df['cycle'] = cyc
    return df


def clean(df):
    '''Adults >= 18 with coded exposure, complete PHQ-9 and model covariates.'''
    d = df[df['ridageyr'] >= 18].copy()
    # exposure: 1 = flu-like illness in past 30 d, 0 = none (7/9 -> missing)
    d['flu'] = np.where(d['hsq520'].isin([1, 2]), (d['hsq520'] == 1).astype(float), np.nan)
    # PHQ-9: item scores 0-3 only, all nine items required
    items = [i for i in PHQ_ITEMS if i in d.columns]
    valid = d[items].apply(lambda s: s.isin([0, 1, 2, 3]))
    complete = valid.all(axis = 1)
    d['phq9'] = np.where(complete, d[items].sum(axis = 1), np.nan)
    d['phq9_ge10'] = np.where(complete, (d['phq9'] >= 10).astype(float), np.nan)
    d['male'] = (d['riagendr'] == 1).astype(float)
    d.loc[~d['riagendr'].isin([1, 2]), 'male'] = np.nan
    d['smoker'] = np.where((d['smq020'] == 1) & (d['smq040'].isin([1, 2])), 1.0,
                   np.where(d['smq020'] == 2, 0.0, np.nan))
    # education and smoking as categorical with their own missing level;
    # poverty ratio median-imputed with a missingness indicator
    d['educ_c'] = pd.cut(d['dmdeduc2'], bins = [0, 2, 3, 4, 5],
                         labels = ['lt_hs', 'hs', 'some_col', 'college']).astype(object)
    d.loc[~d['dmdeduc2'].isin([1, 2, 3, 4, 5]), 'educ_c'] = 'missing'
    d['smoke_c'] = np.where(d['smoker'] == 1, 'current',
                    np.where(d['smoker'] == 0, 'non_current', 'missing'))
    d['pir'] = d['indfmpir']
    d['pir_miss'] = d['pir'].isna().astype(float)
    d['pir'] = d['pir'].fillna(d['pir'].median())
    d['logcrp'] = np.where(d['crp_mgl'] > 0, np.log(d['crp_mgl']), np.nan)
    d['wt'] = d['wtmec2yr']
    return d


def design_matrix(d, model, pooled):
    '''Build the design matrix of one cross-sectional model specification.'''
    cols = [d['flu'].rename('flu'), d['ridageyr'].rename('age'), d['male'].rename('male')]
    race = pd.get_dummies(d['race'].astype('Int64'), prefix = 'race', drop_first = False)
    race = race.astype(float)
    race.columns = ['race_%s' % c.split('_')[-1] for c in race.columns]
    X = pd.concat(cols, axis = 1)
    if pooled:
        cyc = pd.get_dummies(d['cycle'], prefix = 'cyc', drop_first = True).astype(float)
        X = pd.concat([X, race.iloc[:, 1:], cyc], axis = 1)
    else:
        X = pd.concat([X, race.iloc[:, 1:]], axis = 1)
    if model in ('M2', 'M3'):
        edu = pd.get_dummies(d['educ_c'], prefix = 'edu').astype(float).iloc[:, 1:]
        smo = pd.get_dummies(d['smoke_c'], prefix = 'smo').astype(float).iloc[:, 1:]
        X = pd.concat([X, d['pir'].rename('pir'), d['pir_miss'].rename('pir_miss'),
                       edu, smo], axis = 1)
    if model == 'M3':
        X = pd.concat([X, d['logcrp'].rename('logcrp')], axis = 1)
    X.insert(0, 'const', 1.0)
    return X.astype(float)


def fit_cross_sectional(d, model, outcome, pooled, label):
    '''Fit one weighted linear/logistic model; return the exposure row as a dict.'''
    need = ['flu', 'ridageyr', 'male', 'race', 'wt']
    if model == 'M3':
        need += ['logcrp']
    dd = d.dropna(subset = need + [outcome]).copy()
    if pooled:
        dd['wt'] = dd['wt'] / N_CYC
    X = design_matrix(dd, model, pooled)
    y = dd[outcome].astype(float)
    w = dd['wt'].astype(float)
    n_flu = int((dd['flu'] == 1).sum())
    if outcome == 'phq9':
        res = sm.WLS(y, X, weights = w).fit(cov_type = 'HC1')
    else:
        res = sm.GLM(y, X, family = sm.families.Binomial(), var_weights = w) \
                 .fit(cov_type = 'HC1')
    b = res.params['flu']; se = res.bse['flu']; p = res.pvalues['flu']
    out = dict(analysis = 'cross_sectional', outcome = outcome, cycle = label,
               model = model, n = int(len(dd)), n_flu = n_flu,
               log_est = b, ci_lo = b - 1.96 * se, ci_hi = b + 1.96 * se, p = p)
    if outcome == 'phq9_ge10':
        out['est'] = np.exp(b)
        out['ci_lo'] = np.exp(b - 1.96 * se)
        out['ci_hi'] = np.exp(b + 1.96 * se)
        out['est_type'] = 'OR'
    else:
        out['est'] = b
        out['est_type'] = 'beta'
    return out


def run_cross_sectional(all_df):
    '''Per-cycle and pooled weighted models for both PHQ-9 outcomes.'''
    rows = []
    for cyc, d in all_df.items():
        d = clean(d)
        for outcome in OUTCOMES:
            for model in MODELS:
                if model == 'M3' and not d['logcrp'].notna().any():
                    continue
                try:
                    rows.append(fit_cross_sectional(d, model, outcome, False, cyc))
                except Exception as e:
                    print('[SKIP] %s %s %s: %s' % (cyc, outcome, model, e))
    pooled = pd.concat([clean(d) for d in all_df.values()], ignore_index = True)
    for outcome in OUTCOMES:
        for model in MODELS:
            try:
                rows.append(fit_cross_sectional(pooled, model, outcome, True,
                                                'Pooled_2005_2018'))
            except Exception as e:
                print('[SKIP] pooled %s %s: %s' % (outcome, model, e))
    return pd.DataFrame(rows)


def load_mortality(cyc):
    '''Read the NDI public mortality follow-up file of one cycle.'''
    f = 'nhanes_%s_mort_2019_public.tsv' % cyc.replace('-', '_')
    m = read_tsv(os.path.join(DATADIR, 'mort', f))
    for c in ['eligstat', 'mortstat', 'ucod_leading', 'permth_int']:
        m[c] = pd.to_numeric(m[c], errors = 'coerce')
    return m


def fit_cox(dd, exposure, label, subset):
    '''Fit one cycle-stratified Cox model; return the exposure row as a dict.'''
    need = ['time_y', 'event', exposure, 'ridageyr', 'male', 'race']
    d = dd.dropna(subset = need).copy()
    X = pd.concat([d[exposure].rename('expo'), d['ridageyr'].rename('age'),
                   d['male'].rename('male'), d['pir'].rename('pir'),
                   d['pir_miss'].rename('pir_miss')], axis = 1)
    race = pd.get_dummies(d['race'].astype('Int64'), prefix = 'race').astype(float)
    race.columns = ['race_%s' % c.split('_')[-1] for c in race.columns]
    edu = pd.get_dummies(d['educ_c'], prefix = 'edu').astype(float).iloc[:, 1:]
    smo = pd.get_dummies(d['smoke_c'], prefix = 'smo').astype(float).iloc[:, 1:]
    X = pd.concat([X, race.iloc[:, 1:], edu, smo], axis = 1).astype(float)
    ev = int(d['event'].sum())
    if ev < 5:
        print('[SKIP] cox %s %s: only %d events' % (label, subset, ev))
        return None
    b, se, p = cox_fit(d, X)
    if b is None or not np.isfinite(b):
        print('[WARN] cox non-estimable: %s %s' % (label, subset))
        return None
    return dict(analysis = 'flu_pneumonia_mortality', exposure = exposure,
                subset = subset, cycle = label, n = int(len(d)), events = ev,
                HR = np.exp(b), ci_lo = np.exp(b - 1.96 * se),
                ci_hi = np.exp(b + 1.96 * se), p = p)


def cox_fit(d, X):
    '''Newton fit of the PHReg partial likelihood with eigen-based covariance.

    Same update rule as statsmodels PHReg.fit(method='newton'); the stock fit
    ends in an LU Hessian inversion that raises LinAlgError on this
    near-singular design, so the observed information is inverted by eigen
    decomposition instead (point estimates are unchanged).
    '''
    ph = PHReg(d['time_y'].astype(float).values, X.values,
               status = d['event'].astype(float).values,
               strata = d['cyc_code'].astype(str).values)
    beta = np.zeros(X.shape[1])
    for _ in range(COX_MAXITER):
        H = -ph.hessian(beta)
        H[np.diag_indices(H.shape[0])] += COX_RIDGE
        step = np.linalg.solve(H, -ph.score(beta))
        new_beta = beta - step
        done = not np.any(np.abs(new_beta - beta) > COX_TOL)
        beta = new_beta
        if done:
            break
    if not np.isfinite(beta[0]):
        return None, None, None
    H = -ph.hessian(beta)
    w, V = np.linalg.eigh(H)
    keep = w > COX_EIG_RTOL * w.max()
    cov = (V[:, keep] / w[keep]) @ V[:, keep].T
    se = np.sqrt(cov[0, 0])
    p = 2 * norm.sf(abs(beta[0] / se))
    return beta[0], se, p


def run_mortality(all_df):
    '''Link PHQ-9 completers to NDI and fit pooled and per-cycle Cox models.'''
    rows = []
    parts = []
    for cyc, d0 in all_df.items():
        d = clean(d0)
        d = d[d['phq9'].notna()]
        m = load_mortality(cyc)
        dd = d.merge(m[['seqn', 'eligstat', 'mortstat', 'ucod_leading', 'permth_int']],
                     on = 'seqn', how = 'inner')
        dd = dd[dd['eligstat'] == 1].copy()
        dd['time_y'] = dd['permth_int'] / 12.0
        dd['event'] = ((dd['mortstat'] == 1) & (dd['ucod_leading'] == 8)).astype(float)
        dd['cyc_code'] = cyc
        parts.append(dd)
        print('  [OK] mort %s: linked n=%d, flu/pneumonia deaths=%d' %
              (cyc, len(dd), int(dd['event'].sum())))
    md = pd.concat(parts, ignore_index = True)

    for subset, dd in [('all', md), ('no_baseline_flu', md[md['flu'] != 1])]:
        for exposure in EXPOSURES_COX:
            dd2 = dd.copy()
            dd2['phq9_per5'] = dd2['phq9'] / 5.0
            r = fit_cox(dd2, exposure, 'Pooled_2005_2018', subset)
            if r:
                rows.append(r)
    for cyc in CYCLES:
        dd = md[md['cyc_code'] == cyc].copy()
        dd['phq9_per5'] = dd['phq9'] / 5.0
        r = fit_cox(dd, 'phq9_per5', cyc, 'all')
        if r:
            rows.append(r)
    return pd.DataFrame(rows)


def build_results(cross, mort):
    '''Combine both analysis blocks into one long plain-column results table.'''
    c = cross.copy()
    c['exposure'] = 'flu_like_illness_past30d'
    c['adjustment'] = c['model'].map(ADJUST_LABELS)
    c['n_events'] = np.nan
    c = c.rename(columns = {'n_flu': 'n_exposed', 'est': 'estimate', 'p': 'p_value'})
    m = mort.copy()
    m['model'] = 'cox_ph_cycle_stratified'
    m['outcome'] = 'influenza_pneumonia_death'
    m['n_exposed'] = np.nan
    m = m.rename(columns = {'subset': 'adjustment', 'events': 'n_events',
                            'HR': 'estimate', 'p': 'p_value'})
    m['est_type'] = 'HR'
    cols = ['analysis', 'cycle', 'model', 'outcome', 'exposure', 'adjustment',
            'n', 'n_exposed', 'n_events', 'estimate', 'ci_lo', 'ci_hi', 'p_value',
            'est_type']
    return pd.concat([c[cols], m[cols]], ignore_index = True)


def build_plot_data(results):
    '''Add panel and cycle-order columns so R can draw the forest figure.'''
    panel_map = {
        ('cross_sectional', 'phq9'): 'A_phq9_score_linear',
        ('cross_sectional', 'phq9_ge10'): 'B_phq9_ge10_or',
        ('flu_pneumonia_mortality', 'influenza_pneumonia_death'): 'C_flu_pneumonia_death_hr',
    }
    d = results.copy()
    d['panel'] = [panel_map[(a, o)] for a, o in zip(d['analysis'], d['outcome'])]
    d['cycle_order'] = d['cycle'].map({c: i + 1 for i, c in enumerate(CYCLE_ORDER)})
    d['cycle_is_pooled'] = (d['cycle'] == 'Pooled_2005_2018').astype(int)
    cols = ['panel', 'cycle', 'cycle_order', 'cycle_is_pooled', 'analysis', 'model',
            'outcome', 'exposure', 'adjustment', 'estimate', 'ci_lo', 'ci_hi',
            'p_value', 'est_type', 'n', 'n_exposed', 'n_events']
    d = d[cols].sort_values(['panel', 'cycle_order', 'model'], kind = 'stable')
    return d.reset_index(drop = True)


if __name__ == '__main__':

    # ==== preset parameters ====
    FORCE = True                                        # overwrite existing tables
    STEPS = ['load', 'cross_sectional', 'mortality', 'write']

    t0 = time.time()
    os.makedirs(OUTDIR, exist_ok = True)

    all_df = {}
    if 'load' in STEPS:
        print('=' * 70)
        print('[STEP 1/4] load NHANES cycles 2005-2018')
        print('=' * 70)
        for cyc, spec in CYCLES.items():
            raw = load_cycle(cyc, spec)
            all_df[cyc] = raw
            print('[OK] %s: raw n=%d, adults>=18 n=%d' %
                  (cyc, len(raw), int((raw['ridageyr'] >= 18).sum())))

    cross = pd.DataFrame()
    if 'cross_sectional' in STEPS:
        print('=' * 70)
        print('[STEP 2/4] analysis A: flu-like illness x PHQ-9 (cross-sectional)')
        print('=' * 70)
        cross = run_cross_sectional(all_df)
        pooled_desc = pd.concat([clean(d) for d in all_df.values()], ignore_index = True)
        a = pooled_desc.dropna(subset = ['flu', 'phq9', 'wt'])
        print('[OK] pooled adults with exposure+PHQ-9+weight: n=%d' % len(a))
        print('[OK] HSQ520=1 weighted %.2f%% | unweighted %.2f%% | PHQ-9>=10 %.2f%%' % (
            100 * a.loc[a['flu'] == 1, 'wt'].sum() / a['wt'].sum(),
            100 * (a['flu'] == 1).mean(), 100 * a['phq9_ge10'].mean()))

    mort = pd.DataFrame()
    if 'mortality' in STEPS:
        print('=' * 70)
        print('[STEP 3/4] analysis B: PHQ-9 -> influenza/pneumonia mortality (NDI)')
        print('=' * 70)
        mort = run_mortality(all_df)

    if 'write' in STEPS:
        print('=' * 70)
        print('[STEP 4/4] write results and plot-data tables')
        print('=' * 70)
        results = build_results(cross, mort)
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
