# -*- coding: utf-8 -*-

import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

warnings.filterwarnings('ignore', category = FutureWarning)

PROJDIR = '/xbio/articles/influenza_mental_genomics'
COHORT = Path(PROJDIR) / '.xbio' / 'UKB_NHANES' / 'temp' / 'ukb_cohort.parquet'
OUTDIR = Path(PROJDIR) / 'results' / 'UKB_NHANES'
RESULTS_TSV = OUTDIR / 'UKB.uri_pne_psych_incidence.results.tsv'
PLOT_TSV = OUTDIR / 'UKB.uri_pne_psych_incidence.plot_data.tsv'

CENSOR = pd.Timestamp('2023-07-01')                   # admin censor = latest event date in the extract
FLU_COLS = ['p131438', 'p131440', 'p131442']          # first-reported J09 / J10 / J11
URI_COLS = ['p131424', 'p131426', 'p131428', 'p131430',
            'p131432', 'p131434', 'p131436']          # first-reported J00-J06
PNE_COLS = ['p131444', 'p131446', 'p131448', 'p131450',
            'p131452', 'p131454', 'p131456']          # first-reported J12-J18
EXPOSURES = {                                         # exposure -> (date columns, contrast label)
    'flu': dict(cols = FLU_COLS, contrast = 'incident flu vs never-flu'),
    'uri': dict(cols = URI_COLS, contrast = 'incident uri vs never-uri'),
    'pne': dict(cols = PNE_COLS, contrast = 'incident pne vs never-pne'),
}
PSYCH = {                                             # first-reported psychiatric ICD-10 blocks
    'F17_tobacco': 'p130868',
    'F20_schizophrenia': 'p130874',
    'F25_schizoaffective': 'p130884',
    'F30_mania': 'p130890',
    'F31_bipolar': 'p130892',
    'F32_dep_single': 'p130894',
    'F33_dep_recur': 'p130896',
    'F34_persistent_mood': 'p130898',
    'F38_other_mood': 'p130900',
    'F39_unspec_mood': 'p130902',
    'F40_phobic': 'p130904',
    'F41_anxiety': 'p130906',
    'F42_ocd': 'p130908',
    'F43_stress_ptsd': 'p130910',
}
KEY_OUTCOMES = ['any_psych', 'F17_tobacco', 'F32_dep_single', 'F33_dep_recur',
                'F41_anxiety', 'F43_stress_ptsd', 'F20_schizophrenia', 'F31_bipolar']
MIN_EVENTS = 5                                        # per-arm event floor for the sensitivity battery
STRATIFIERS = ['gp_nerves', 'psych_spec', 'hosp_ever']

# model formula fragments (rate ratios from Poisson models; S4 is a 2-year window logit);
# the exposure name and the other-infection priors are plugged in per exposure
F_COV = ' + age_index + C(sex) + bmi + C(ethnic)'
F_LIFE = ' + C(smoke) + C(alcohol)'
F_HC = ' + C(gp_nerves) + C(psych_spec)'
F_YEAR = ' + C(index_year)'

RESULT_COLS = ['section', 'exposure', 'analysis', 'outcome', 'term', 'contrast', 'estimate_type',
               'estimate', 'ci_lo', 'ci_hi', 'p', 'n_total', 'n_events', 'n_exposed',
               'events_exposed', 'py_exposed', 'events_unexposed', 'py_unexposed',
               'crude_rr', 'stratifier', 'stratum', 'note']


def eth_grp(x):
    '''Map UKB ethnicity code p21000_i0 to 6 broad groups; missing/negative -> missing.'''
    if pd.isna(x) or x < 0:
        return 'missing'
    if 1000 <= x < 2000:
        return 'white'
    if 2000 <= x < 3000:
        return 'mixed'
    if 3000 <= x < 4000:
        return 'south_asian'
    if 4000 <= x < 5000:
        return 'black'
    if x == 5:
        return 'chinese'
    return 'other'


def exposure_formulas(expo, others):
    '''Model sequence for one exposure; mutual = primary + other infections strictly before index.'''
    f_unadj = 'events ~ %s' % expo
    f_nolife = f_unadj + F_COV
    f_primary = f_nolife + F_LIFE
    f_mutual = f_primary + ''.join(' + %s_prior' % o for o in sorted(others))
    f_hc = f_primary + F_HC
    f_year = f_primary + F_YEAR
    f_logit = f_primary
    return dict(unadj = f_unadj, nolife = f_nolife, primary = f_primary, mutual = f_mutual,
                hc = f_hc, year = f_year, logit = f_logit)


def build_covariates(df, index, expo):
    '''Covariate frame for one exposure; others = 0/1 indicators for the other
    respiratory infections first diagnosed strictly before the index date.'''
    base = pd.to_datetime(df['p53_i0'], errors = 'coerce')
    others = [e for e in EXPOSURES if e != expo]
    cov = pd.DataFrame(index = df.index)
    cov[expo] = df[expo]
    for o in others:
        cov['%s_prior' % o] = (df['%s_first' % o] < index).astype(int)
    cov['age_index'] = pd.to_numeric(df['p21003_i0'], errors = 'coerce') + (index - base).dt.days / 365.25
    cov['sex'] = pd.to_numeric(df['p31'], errors = 'coerce').map({0: 'F', 1: 'M'}).astype('category')
    cov['bmi'] = pd.to_numeric(df['p21001_i0'], errors = 'coerce')
    cov['ethnic'] = pd.to_numeric(df['p21000_i0'], errors = 'coerce').map(eth_grp).astype('category')
    cov['smoke'] = (pd.to_numeric(df['p20116_i0'], errors = 'coerce')
                    .map({0: 'never', 1: 'prev', 2: 'curr'}).fillna('missing').astype('category'))
    cov['alcohol'] = (pd.to_numeric(df['p1558_i0'], errors = 'coerce').where(lambda s: s > 0)
                      .map({1: 'daily', 2: 'wk34', 3: 'wk12', 4: 'monthly', 5: 'special', 6: 'never'})
                      .fillna('missing').astype('category'))
    gp = pd.to_numeric(df['p2090_i0'], errors = 'coerce')
    sp = pd.to_numeric(df['p2100_i0'], errors = 'coerce')
    cov['gp_nerves'] = gp.where(gp.isin([0, 1])).map({0: 'no', 1: 'yes'}).fillna('missing').astype('category')
    cov['psych_spec'] = sp.where(sp.isin([0, 1])).map({0: 'no', 1: 'yes'}).fillna('missing').astype('category')
    cov['index_year'] = index.dt.year.astype('category')
    cov['hosp_ever'] = (pd.to_numeric(df['p41235'], errors = 'coerce').fillna(0) > 0).map(
        {True: 'ever', False: 'never'})
    return cov


def load_cohort():
    '''Read the cohort parquet; derive first-event dates, exposure flags and index dates.'''
    df = pd.read_parquet(COHORT)
    print('[OK] cohort loaded', df.shape, 'from', COHORT)
    df['base'] = pd.to_datetime(df['p53_i0'], errors = 'coerce')
    df['death'] = pd.to_datetime(df['p40000_i0'], errors = 'coerce')
    for expo, spec in EXPOSURES.items():
        df['%s_first' % expo] = df[spec['cols']].apply(pd.to_datetime, errors = 'coerce').min(axis = 1)
        exposed = df['%s_first' % expo].notna() & (df['%s_first' % expo] > df['base'])
        prevalent = df['%s_first' % expo].notna() & (df['%s_first' % expo] <= df['base'])
        df[expo] = exposed.astype(int)
        print('[OK] %-3s incident after baseline (exposed): %-7d | prevalent excluded: %-7d '
              '| never-exposed pool: %d'
              % (expo, int(exposed.sum()), int(prevalent.sum()),
                 int(df['%s_first' % expo].isna().sum())))
    psych = pd.DataFrame({k: pd.to_datetime(df[c], errors = 'coerce') for k, c in PSYCH.items()})
    df['end'] = pd.concat([df['death'], pd.Series(CENSOR, index = df.index)], axis = 1).min(axis = 1)
    return df, psych


def empty_result(n = 0, events = 0, note = ''):
    '''All-NaN result dict for contrasts that are not estimable.'''
    return dict(estimate = np.nan, ci_lo = np.nan, ci_hi = np.nan, p = np.nan,
                n = n, events = events, note = note)


def fit_poisson(d, formula, term):
    '''Poisson GLM with log(person-year) offset and HC0 sandwich SEs; returns the exposure RR.'''
    d = d.dropna(subset = ['py', 'events', term])
    n, events = len(d), int(d['events'].sum())
    if events == 0 or d[term].nunique() < 2:
        return empty_result(n = n, events = events)
    if d.loc[d[term] == 1, 'events'].sum() == 0:
        return empty_result(n = n, events = events, note = 'zero_events_exposed')
    if d.loc[d[term] == 0, 'events'].sum() == 0:
        return empty_result(n = n, events = events, note = 'zero_events_unexposed')
    try:
        m = smf.glm(formula, data = d, family = sm.families.Poisson(),
                    offset = np.log(d['py'])).fit(cov_type = 'HC0')
        return dict(estimate = float(np.exp(m.params[term])),
                    ci_lo = float(np.exp(m.conf_int().loc[term, 0])),
                    ci_hi = float(np.exp(m.conf_int().loc[term, 1])),
                    p = float(m.pvalues[term]), n = n, events = events,
                    note = '' if m.converged else 'not_converged')
    except Exception as e:
        return empty_result(n = n, events = events, note = 'fit_failed:%s' % e)


def fit_logit_window(d, formula, term, years = 2):
    '''Logistic model on a fixed post-index window; returns the exposure odds ratio.'''
    d = d.dropna(subset = ['events', term])
    n, events = len(d), int(d['events'].sum())
    if events == 0 or d[term].nunique() < 2:
        return empty_result(n = n, events = events)
    if d.loc[d[term] == 1, 'events'].sum() == 0:
        return empty_result(n = n, events = events, note = 'zero_events_exposed')
    if d.loc[d[term] == 0, 'events'].sum() == 0:
        return empty_result(n = n, events = events, note = 'zero_events_unexposed')
    try:
        m = smf.logit(formula, data = d).fit(disp = 0)
        return dict(estimate = float(np.exp(m.params[term])),
                    ci_lo = float(np.exp(m.conf_int().loc[term, 0])),
                    ci_hi = float(np.exp(m.conf_int().loc[term, 1])),
                    p = float(m.pvalues[term]), n = n, events = events,
                    note = '' if m.mle_retvals.get('converged', True) else 'not_converged')
    except Exception as e:
        return empty_result(n = n, events = events, note = 'fit_failed:%s' % e)


def model_row(section, analysis, outcome, res, expo, term, contrast, **extra):
    '''Expand one fit result dict into a long-format result row.'''
    r = dict(section = section, exposure = expo, analysis = analysis, outcome = outcome,
             term = term, contrast = contrast, estimate_type = 'RR',
             estimate = res['estimate'], ci_lo = res['ci_lo'], ci_hi = res['ci_hi'],
             p = res['p'], n_total = res['n'], n_events = res['events'], n_exposed = None,
             events_exposed = None, py_exposed = None, events_unexposed = None,
             py_unexposed = None, crude_rr = None, stratifier = None, stratum = None,
             note = res.get('note', ''))
    r.update(extra)
    return r


def fit_outcome(oname, odate, df, cov, index, keep, expo, spec):
    '''Fit every pre-specified model for one outcome and exposure; returns the long rows.'''
    F = exposure_formulas(expo, [e for e in EXPOSURES if e != expo])
    d = cov.loc[keep].copy()
    idx = index[keep]
    end = df['end'][keep]
    o = odate[keep]
    atrisk = o.isna() | (o > idx)                      # incident-risk set for this diagnosis
    d = d.loc[atrisk]
    idx = idx[atrisk]
    end = end[atrisk]
    o = o[atrisk]
    d['events'] = ((o.notna()) & (o <= end)).astype(int)
    d['py'] = (end - idx).dt.days / 365.25
    d = d.loc[d['py'] > 0]
    n_tot = len(d)
    ev_tot = int(d['events'].sum())
    n_exp = int(d[expo].sum())
    ev_exp = int(d.loc[d[expo] == 1, 'events'].sum())
    py_exp = d.loc[d[expo] == 1, 'py'].sum()
    ev_un = int(d.loc[d[expo] == 0, 'events'].sum())
    py_un = d.loc[d[expo] == 0, 'py'].sum()
    crude = np.nan
    if py_exp > 0 and py_un > 0 and ev_un > 0 and ev_exp > 0:
        crude = (ev_exp / py_exp) / (ev_un / py_un)

    # main block: crude rate ratio, unadjusted, no-lifestyle, primary and mutually adjusted RRs
    if ev_exp == 0 or ev_un == 0:
        znote = 'zero_events_exposed' if ev_exp == 0 else 'zero_events_unexposed'
        r0 = r1 = empty_result(n_tot, ev_tot)
        rP = empty_result(n_tot, ev_tot, znote)
        rM = empty_result(n_tot, ev_tot, znote)
    else:
        r0 = fit_poisson(d, F['unadj'], expo)
        r1 = fit_poisson(d, F['nolife'], expo)
        rP = fit_poisson(d, F['primary'], expo)
        rM = fit_poisson(d, F['mutual'], expo)
    desc = dict(n_exposed = n_exp, events_exposed = ev_exp, py_exposed = py_exp,
                events_unexposed = ev_un, py_unexposed = py_un, crude_rr = crude)
    rows = [model_row('main', 'crude', oname, empty_result(n_tot, ev_tot), expo, expo,
                      spec['contrast'], estimate = crude, **desc),
            model_row('main', 'unadjusted', oname, r0, expo, expo, spec['contrast'], **desc),
            model_row('main', 'no_lifestyle', oname, r1, expo, expo, spec['contrast'], **desc),
            model_row('main', 'primary', oname, rP, expo, expo, spec['contrast'], **desc),
            model_row('main', 'mutual', oname, rM, expo, expo, spec['contrast'], **desc)]
    if pd.notna(rP['estimate']):
        print('[OK] %-4s %-18s N=%d exposed=%d events=%d RR=%.3f (%.3f-%.3f) p=%.3g '
              '| mutual RR=%.3f'
              % (expo, oname, n_tot, n_exp, ev_tot, rP['estimate'], rP['ci_lo'], rP['ci_hi'],
                 rP['p'], rM['estimate'] if pd.notna(rM['estimate']) else np.nan))
    else:
        print('[WARN] %-4s %-18s N=%d exposed=%d events=%d not estimable (%s)'
              % (expo, oname, n_tot, n_exp, ev_tot, rP.get('note', '')))

    # sensitivity battery and contact stratification for the key outcomes
    if oname in KEY_OUTCOMES and ev_exp >= MIN_EVENTS and ev_un >= MIN_EVENTS:
        rows.append(model_row('sensitivity', 'S0_unadjusted', oname, r0, expo, expo, spec['contrast']))
        rows.append(model_row('sensitivity', 'S1_nolifestyle', oname, r1, expo, expo, spec['contrast']))
        rows.append(model_row('sensitivity', 'S2_plus_healthcare_contact', oname,
                              fit_poisson(d, F['hc'], expo), expo, expo, spec['contrast']))
        rows.append(model_row('sensitivity', 'S3_plus_index_year', oname,
                              fit_poisson(d, F['year'], expo), expo, expo, spec['contrast']))
        # S4: logistic model on a fixed 2-year post-index window (equal opportunity window)
        idx_full = index[keep][atrisk]
        end_full = df['end'][keep][atrisk]
        full_window = end_full >= (idx_full + pd.DateOffset(years = 2))
        d2 = d.loc[full_window.values].copy()
        o2 = o.loc[full_window.values]
        idx2 = idx_full[full_window.values]
        d2['events'] = ((o2.notna()) & (o2 > idx2) & (o2 <= idx2 + pd.DateOffset(years = 2))).astype(int)
        rows.append(model_row('sensitivity', 'S4_logit_2y_window', oname,
                              fit_logit_window(d2, F['logit'], expo), expo, expo,
                              spec['contrast'], estimate_type = 'OR'))
        # S5: primary model restricted to participants with >= 2 years potential follow-up
        d5 = d.loc[((end_full - idx_full).dt.days >= 730).values]
        rows.append(model_row('sensitivity', 'S5_min2y_followup', oname,
                              fit_poisson(d5, F['primary'], expo), expo, expo, spec['contrast']))
        # healthcare-contact stratification (same adjustment inside strata)
        for strat in STRATIFIERS:
            for lvl in sorted(d[strat].dropna().unique()):
                rs = fit_poisson(d.loc[d[strat] == lvl], F['nolife'], expo)
                rows.append(model_row('strata', 'stratified', oname, rs, expo, expo,
                                      spec['contrast'], stratifier = strat, stratum = lvl))
    return rows


def fit_all(df, psych):
    '''Fit every outcome across the three exposures and return the combined long table.'''
    rows = []
    outcomes = {'any_psych': psych.min(axis = 1)}
    outcomes.update({k: psych[k] for k in PSYCH})
    for expo, spec in EXPOSURES.items():
        first = df['%s_first' % expo]
        index = pd.Series(pd.NaT, index = df.index)
        index[df[expo] == 1] = first[df[expo] == 1]
        index[first.isna()] = df['base'][first.isna()]
        keep = index.notna() & (df['end'] > index)
        cov = build_covariates(df, index, expo)
        print('=' * 70)
        print('[STEP] exposure = %-3s | primary analysis population (index < censor): %d'
              % (expo, int(keep.sum())))
        print('=' * 70)
        for oname, odate in outcomes.items():
            rows.extend(fit_outcome(oname, odate, df, cov, index, keep, expo, spec))
    return pd.DataFrame(rows)[RESULT_COLS]


def plot_row(panel, order, label, r):
    '''One plot-table row (full precision) from a result row.'''
    return dict(panel = panel, row_order = order, label = label, exposure = r['exposure'],
                section = r['section'], analysis = r['analysis'], outcome = r['outcome'],
                term = r['term'], contrast = r['contrast'], estimate_type = r['estimate_type'],
                estimate = r['estimate'], ci_lo = r['ci_lo'], ci_hi = r['ci_hi'], p = r['p'],
                n_total = r['n_total'], n_events = r['n_events'], note = r['note'])


def plot_frame(res):
    '''Long plot table: main forest panels sorted by the flu primary RR (shared order),
    then sensitivity and strata panels.'''
    rows = []
    order = 0
    flu_main = res[(res['section'] == 'main') & (res['analysis'] == 'primary')
                   & (res['exposure'] == 'flu')].dropna(subset = ['estimate']
                                                        ).sort_values('estimate')
    order_map = {o: i for i, o in enumerate(flu_main['outcome'])}
    for expo in EXPOSURES:
        main = res[(res['section'] == 'main') & (res['analysis'] == 'primary')
                   & (res['exposure'] == expo)].copy()
        main['_o'] = main['outcome'].map(order_map)
        for _, r in main.sort_values('_o').iterrows():
            order += 1
            rows.append(plot_row('main_forest', order, r['outcome'], r))
    for _, r in res[res['section'] == 'sensitivity'].iterrows():
        order += 1
        rows.append(plot_row('sensitivity', order,
                             '%s | %s | %s' % (r['exposure'], r['outcome'], r['analysis']), r))
    for _, r in res[res['section'] == 'strata'].iterrows():
        order += 1
        rows.append(plot_row('strata', order,
                             '%s | %s | %s = %s'
                             % (r['exposure'], r['outcome'], r['stratifier'], r['stratum']), r))
    return pd.DataFrame(rows)


def write_tables(res, plot):
    '''Write results.tsv (rounded values) and plot_data.tsv (full precision).'''
    out = res.copy()
    for c in ['estimate', 'ci_lo', 'ci_hi', 'crude_rr']:
        out[c] = pd.to_numeric(out[c], errors = 'coerce').round(4)
    for c in ['py_exposed', 'py_unexposed']:
        out[c] = pd.to_numeric(out[c], errors = 'coerce').round(1)
    for c in ['n_exposed', 'events_exposed', 'events_unexposed']:
        out[c] = out[c].astype('Int64')
    out.to_csv(RESULTS_TSV, sep = '\t', index = False)
    print('[OK] written', RESULTS_TSV)
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
            print('[STEP 1/2] fit incident respiratory infection -> incident psychiatric diagnosis')
            print('=' * 70)
            df, psych = load_cohort()
            res = fit_all(df, psych)
            cols = ['exposure', 'analysis', 'outcome', 'estimate', 'ci_lo', 'ci_hi',
                    'p', 'n_total', 'n_events']
            print(res[res['section'] == 'main'][cols].to_string(index = False))
        if 'tables' in STEPS:
            print('=' * 70)
            print('[STEP 2/2] write results and plot tables')
            print('=' * 70)
            write_tables(res, plot_frame(res))

    print('=' * 70)
    print('[DONE] analysis finished in %.1f min' % ((time.time() - t0) / 60))
    print('=' * 70)
