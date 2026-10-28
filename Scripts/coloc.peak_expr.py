#!/usr/bin/env python3

import argparse
import os
from concurrent.futures import ThreadPoolExecutor

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import statsmodels.api as sm
from pyliftover import LiftOver
from scipy.sparse import csr_matrix
from scipy.stats import spearmanr

BASE = "/home/hello/influenza_mental_genomics"
RNA = os.path.join(BASE, "data/scrna/SEAAD_MTG_RNAseq_final-nuclei.2026-06-22.h5ad")
ATAC = os.path.join(BASE, "data/SEAAD_ATAC/SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad")
GENCODE = os.path.join(BASE, "references/gencode_hg19_genes.txt")
CHAIN = "/home/hello/wkdir/influenza_mental_genomics_ai/data/liftOver/hg19ToHg38.over.chain.gz"

CHR, SNP_POS = "chr3", 186_933_001           # rs10937275
LEAD_POS = 186_977_520                       # ST6GAL1 GWAS lead SNP (negative control)
TARGET_GENE = "ST6GAL1"
WINDOW = 1_000_000                           # gene specificity window

# ATAC subclass -> RNA subclass (SEA-AD naming differs across assays)
CELL_PAIRS = [
    ("L6 IT", "L6 IT"),
    ("Oligodendrocyte", "Oligodendrocyte"),
    ("Astrocyte", "Astrocyte"),
    ("Microglia-PVM", "Immune"),             # RNA aggregates microglia/PVM as "Immune"
]


def load_obs(path):
    """Load obs dataframe (and var_names) from a backed h5ad."""
    a = sc.read_h5ad(path, backed="r")
    obs = a.obs.copy()
    var_names = np.array(a.var_names)
    a.file.close()
    return obs, var_names


def parse_peak(name):
    c, rest = name.split(":")
    s, e = rest.split("-")
    return c, int(s), int(e)


def find_covering_peak(var_names, pos, chrom=CHR):
    """Return (peak name, distance) for the peak covering pos (else nearest)."""
    cands = []
    for v in var_names:
        if not v.startswith(chrom + ":"):
            continue
        _, s, e = parse_peak(v)
        if s <= pos <= e:
            return v, 0
        cands.append((min(abs(pos - s), abs(pos - e)), v))
    cands.sort()
    return cands[0][1], cands[0][0]


class CSRReader:
    """Row gatherer over a CSR matrix stored in an h5ad group.

    The full data/indices arrays are read once (sequential h5py read, fast)
    and cached; row subsets are then sliced in numpy. Indices are cast to
    int32 to halve memory.

    group="X" is the (normalized) matrix; group="layers/UMIs" is the raw
    RNA count layer used for pseudobulk expression.
    """

    def __init__(self, h5ad_path, group="X"):
        self.f = h5py.File(h5ad_path, "r")
        self.x = self.f[group]
        self.indptr = self.x["indptr"][:]
        self.shape = tuple(self.x["shape"][:]) if "shape" in self.x \
            else tuple(self.x.attrs["shape"])
        self._data = None
        self._indices = None

    def _load(self):
        if self._data is None:
            print(f"    loading CSR arrays: nnz={self.x['data'].shape[0]:,}",
                  flush=True)
            self._data = self.x["data"][:]
            self._indices = self.x["indices"][:].astype(np.int32)

    def gather(self, rows, chunk=20000):
        """Return (sub_csr, totals): CSR of selected rows and per-row sums."""
        self._load()
        out_data, out_idx = [], []
        out_ptr = [0]
        totals = np.empty(len(rows), dtype=np.float64)
        for lo in range(0, len(rows), chunk):
            sel = rows[lo:lo + chunk]
            starts, ends = self.indptr[sel], self.indptr[sel + 1]
            lens = ends - starts
            offs = np.zeros(len(lens) + 1, dtype=np.int64)
            np.cumsum(lens, out=offs[1:])
            idx_arr = (np.arange(offs[-1]) - np.repeat(offs[:-1], lens)
                       + np.repeat(starts, lens))
            d = self._data[idx_arr]
            out_data.append(d)
            out_idx.append(self._indices[idx_arr])
            out_ptr.extend((out_ptr[-1] + offs[1:]).tolist())
            totals[lo:lo + len(sel)] = np.add.reduceat(
                d, offs[:-1]) if len(d) else 0.0
        sub = csr_matrix(
            (np.concatenate(out_data), np.concatenate(out_idx),
             np.array(out_ptr, dtype=np.int64)),
            shape=(len(rows), self.shape[1]),
        )
        return sub, totals

    def close(self):
        self._data = None
        self._indices = None
        self.f.close()


def group_sum(sub, totals, codes, n_groups):
    """Per-group sparse column sums and weighted totals."""
    g = csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))),
                   shape=(n_groups, len(codes)))
    return g @ sub, np.bincount(codes, weights=totals, minlength=n_groups)


def fit_test(df, acc_col, expr_col):
    """OLS expr ~ acc + age + sex + ADNC + log depth + log n_cells; plus Spearman."""
    adnc_map = {"Not AD": 0, "Low": 1, "Intermediate": 2, "High": 3}
    d = df.dropna(subset=[acc_col, expr_col]).copy()
    d["adnc"] = d["adnc"].map(adnc_map)
    X = np.column_stack([
        d[acc_col].values.astype(float),
        pd.to_numeric(d["age"], errors="coerce").values,
        (d["sex"].astype(str).str.upper().str.startswith("M")).astype(float).values,
        d["adnc"].values.astype(float),
        np.log10(d["depth"].values + 1.0),
        np.log10(d["n_cells"].values),
    ])
    ok = ~np.isnan(X).any(axis=1) & ~np.isnan(d[expr_col].values)
    X, y = X[ok], d[expr_col].values[ok]
    n = len(y)
    if n < 10 or np.std(X[:, 0]) == 0:
        return dict(n=n, beta=np.nan, se=np.nan, p=np.nan,
                    rho=np.nan, p_spear=np.nan)
    res = sm.OLS(y, sm.add_constant(X)).fit()
    rho, p_sp = spearmanr(X[:, 0], y)
    return dict(n=n, beta=res.params[1], se=res.bse[1], p=res.pvalues[1],
                rho=rho, p_spear=p_sp)


def run_tests(tasks, threads):
    """Run (family, celltype, test, df, acc_col, expr_col) tasks in parallel."""
    def one(t):
        fam, ct, name, df, acc_col, expr_col = t
        return dict(family=fam, celltype=ct, test=name,
                    **fit_test(df, acc_col, expr_col))
    if threads <= 1:
        return [one(t) for t in tasks]
    with ThreadPoolExecutor(max_workers=threads) as ex:
        return list(ex.map(one, tasks))


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    q = np.full(len(p), np.nan)
    ok = ~np.isnan(p)
    if ok.sum() == 0:
        return q
    pv = p[ok]
    order = np.argsort(pv)
    ranks = np.empty(len(pv), dtype=int)
    ranks[order] = np.arange(1, len(pv) + 1)
    qv = np.minimum.accumulate((pv * len(pv) / ranks)[order][::-1])[::-1]
    out = np.empty(len(pv))
    out[order] = np.clip(qv, 0, 1)
    q[ok] = out
    return q


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rna", default = RNA)
    parser.add_argument("--atac", default = ATAC)
    parser.add_argument("--outdir", default = os.path.join(BASE, "results/coloc"))
    parser.add_argument("--tempdir", default = os.path.join(BASE, "results.temp/coloc/peak_expr"))
    parser.add_argument("--min-cells", type=int, default = 20)
    parser.add_argument("--n-random-peaks", type=int, default = 100)
    parser.add_argument("--n-candidate-peaks", type=int, default = 400)
    parser.add_argument("--threads", type=int, default = 1,
                        help = "parallel workers for donor-level association tests")
    parser.add_argument("--seed", type=int, default = 42)
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.tempdir, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    # ---------- 1. obs / var ----------
    print("[1/6] loading obs", flush=True)
    rna_obs, rna_var = load_obs(args.rna)
    atac_obs, atac_var = load_obs(args.atac)
    rna_donors = set(rna_obs["Donor ID"].astype(str))
    atac_donors = set(atac_obs["Donor ID"].astype(str))
    common = sorted(rna_donors & atac_donors)
    print(f"donors: RNA={len(rna_donors)} ATAC={len(atac_donors)} common={len(common)}")
    atac_obs = atac_obs[atac_obs["Donor ID"].astype(str).isin(common)].copy()
    rna_obs = rna_obs[rna_obs["Donor ID"].astype(str).isin(common)].copy()

    # ---------- 2. target peaks ----------
    print("[2/6] locating target peaks", flush=True)
    target_peak, d1 = find_covering_peak(atac_var, SNP_POS)
    lead_peak, d2 = find_covering_peak(atac_var, LEAD_POS)
    _, ts, te = parse_peak(target_peak)
    print(f"rs10937275 ({SNP_POS:,}) -> {target_peak} (dist {d1})")
    print(f"lead SNP   ({LEAD_POS:,}) -> {lead_peak} (dist {d2})")
    pd.DataFrame([
        dict(site="rs10937275", pos=SNP_POS, peak=target_peak, dist_bp=d1),
        dict(site="lead_SNP", pos=LEAD_POS, peak=lead_peak, dist_bp=d2),
    ]).to_csv(os.path.join(args.outdir, "peak_expr_target_peaks.tsv"),
              sep="\t", index=False)

    # ---------- 3. cis gene set (+/-1 Mb, gencode hg19 -> hg38) ----------
    print("[3/6] building cis gene set", flush=True)
    cis_path = os.path.join(args.tempdir, "peak_expr_cis_genes.tsv")
    rna_gene_set = set(rna_var)
    mid = (ts + te) // 2
    if os.path.exists(cis_path):
        gdf = pd.read_csv(cis_path, sep="\t")
        print(f"  reusing lifted cis gene table: {cis_path}")
    else:
        lo = LiftOver(CHAIN)
        genc = pd.read_csv(GENCODE, sep="\t")
        genc = genc[genc["chr"] == CHR]
        recs = []
        for r in genc.itertuples(index=False):
            conv = lo.convert_coordinate(r.chr, int((r.start + r.end) // 2))
            if conv:
                recs.append((r.gene_name, int(conv[0][1])))
        gdf = pd.DataFrame(recs, columns=["gene", "pos_hg38"])
        gdf.to_csv(cis_path, sep="\t", index=False)
    gdf = gdf[gdf["pos_hg38"].between(mid - WINDOW, mid + WINDOW)]
    gdf = gdf[gdf["gene"].isin(rna_gene_set)].drop_duplicates("gene")
    cis_genes = sorted(set(gdf["gene"]) | {TARGET_GENE})
    rng_genes = [g for g in rng.choice(rna_var, size=2000, replace=False)
                 if g not in cis_genes][:2]
    print(f"cis genes within +/-1 Mb: {len(cis_genes)}; random controls: {rng_genes}")

    # ---------- 4. ATAC pseudobulk ----------
    print("[4/6] ATAC pseudobulk", flush=True)
    atac_idx = pd.Index(atac_var)
    rna_idx = pd.Index(rna_var)
    t_col = atac_idx.get_loc(target_peak)
    l_col = atac_idx.get_loc(lead_peak)
    gene_cols = {g: rna_idx.get_loc(g) for g in cis_genes + rng_genes
                 if g in rna_idx}

    ar = CSRReader(args.atac)
    subclass_arr = atac_obs["Subclass"].values
    donor_arr = atac_obs["Donor ID"].astype(str).values
    atac_pb = {}
    l6_cache = None
    for at_sub, _ in CELL_PAIRS:
        rows = np.where(subclass_arr == at_sub)[0]
        udon, code = np.unique(donor_arr[rows], return_inverse=True)
        sub, totals = ar.gather(rows)
        if at_sub == "L6 IT":
            l6_cache = (sub, totals, udon, code)
        g_sum, g_tot = group_sum(sub, totals, code, len(udon))
        df = pd.DataFrame({
            "donor": udon,
            "n_cells": np.bincount(code, minlength=len(udon)),
            "depth": g_tot,
            "acc_target": g_sum[:, t_col].toarray().ravel(),
            "acc_lead": g_sum[:, l_col].toarray().ravel(),
        })
        df["acc_target"] = np.log1p(df["acc_target"] / df["depth"] * 1e6)
        df["acc_lead"] = np.log1p(df["acc_lead"] / df["depth"] * 1e6)
        atac_pb[at_sub] = df
        print(f"  ATAC {at_sub}: {len(df)} donors, "
              f">={args.min_cells} cells: {(df['n_cells'] >= args.min_cells).sum()}")

    # matched random peaks (width- and accessibility-matched) on L6 IT
    _, ws, we = parse_peak(target_peak)
    wlen = we - ws
    cand = [v for v in atac_var if abs(parse_peak(v)[2] - parse_peak(v)[1] - wlen) < wlen * 0.5]
    cand = list(rng.choice(np.array(cand, dtype=object),
                           size=min(args.n_candidate_peaks, len(cand)), replace=False))
    cand_cols = np.array([atac_idx.get_loc(c) for c in cand])
    sub, totals, udon, code = l6_cache
    g_sum, g_tot = group_sum(sub, totals, code, len(udon))
    cpm = np.log1p(g_sum[:, cand_cols].toarray() / g_tot[:, None] * 1e6)
    okc = np.isfinite(cpm).all(axis=0) & (cpm.std(axis=0) > 0)
    cpm, cand = cpm[:, okc], np.array(cand, dtype=object)[okc]
    l6_acc = atac_pb["L6 IT"].set_index("donor").loc[udon, "acc_target"].values
    keep = np.argsort(np.abs(cpm.mean(axis=0) - l6_acc.mean()))[:args.n_random_peaks]
    rand_peaks, rand_cpm = cand[keep], cpm[:, keep]
    print(f"matched {len(rand_peaks)} random peaks (target mean logCPM "
          f"{l6_acc.mean():.2f})")
    pd.DataFrame({"peak": rand_peaks, "mean_logcpm": rand_cpm.mean(axis=0)}
                 ).to_csv(os.path.join(args.tempdir, "peak_expr_random_peaks.tsv"),
                          sep="\t", index=False)
    l6_donors = udon

    # random-peak accessibility in the other cell types
    rp_cols = {p: atac_idx.get_loc(p) for p in rand_peaks}
    for at_sub, _ in CELL_PAIRS[1:]:
        rows = np.where(subclass_arr == at_sub)[0]
        udon, code = np.unique(donor_arr[rows], return_inverse=True)
        sub, totals = ar.gather(rows)
        g_sum, g_tot = group_sum(sub, totals, code, len(udon))
        for p, c in rp_cols.items():
            atac_pb[at_sub][f"rp__{p}"] = np.log1p(
                g_sum[:, c].toarray().ravel() / g_tot * 1e6)
    ar.close()

    # ---------- 5. RNA pseudobulk ----------
    print("[5/6] RNA pseudobulk", flush=True)
    rr = CSRReader(args.rna, group="layers/UMIs")
    rna_subclass = rna_obs["Subclass"].values
    rna_donor = rna_obs["Donor ID"].astype(str).values
    umi = rna_obs["Number of UMIs"].astype(float).values
    gene_names = list(gene_cols)
    gene_col_arr = np.array([gene_cols[g] for g in gene_names])
    rna_pb = {}
    for _, rn_sub in CELL_PAIRS:
        rows = np.where(rna_subclass == rn_sub)[0]
        udon, code = np.unique(rna_donor[rows], return_inverse=True)
        sub, _ = rr.gather(rows)
        g_sum = group_sum(sub[:, gene_col_arr], np.zeros(len(rows)), code, len(udon))[0]
        depth = np.bincount(code, weights=umi[rows], minlength=len(udon))
        df = pd.DataFrame({"donor": udon,
                           "rna_n_cells": np.bincount(code, minlength=len(udon)),
                           "rna_depth": depth})
        gcpm = np.log1p(g_sum.toarray() / depth[:, None] * 1e6)
        for j, gname in enumerate(gene_names):
            df[f"expr__{gname}"] = gcpm[:, j]
        rna_pb[rn_sub] = df
        print(f"  RNA {rn_sub}: {len(df)} donors")
    rr.close()

    # ---------- 6. association tests ----------
    print(f"[6/6] association tests (threads={args.threads})", flush=True)
    cov = (rna_obs.groupby(rna_obs["Donor ID"].astype(str))
           [["Age at Death", "Sex", "Overall AD neuropathological Change"]]
           .first().reset_index())
    cov.columns = ["donor", "age", "sex", "adnc"]

    def merged(at_sub, rn_sub):
        m = atac_pb[at_sub].merge(rna_pb[rn_sub], on="donor").merge(cov, on="donor")
        return m[(m["n_cells"] >= args.min_cells)
                 & (m["rna_n_cells"] >= args.min_cells)]

    tasks = []
    m6 = merged("L6 IT", "L6 IT")
    for at_sub, rn_sub in CELL_PAIRS:
        m = merged(at_sub, rn_sub)
        for pk, col in [("target", "acc_target"), ("lead_negctrl", "acc_lead")]:
            tasks.append(("celltype_specificity", at_sub, f"{pk}:{TARGET_GENE}",
                          m, col, f"expr__{TARGET_GENE}"))
    for gname in gene_names:
        tasks.append(("gene_specificity", "L6 IT", f"target:{gname}",
                      m6, "acc_target", f"expr__{gname}"))
    for at_sub, rn_sub in CELL_PAIRS:
        m = merged(at_sub, rn_sub)
        for j, p in enumerate(rand_peaks):
            if at_sub == "L6 IT":
                mm = m.copy()
                mm["rp_tmp"] = rand_cpm[:, j][np.searchsorted(l6_donors, mm["donor"])]
                tasks.append(("background_peaks", at_sub, f"rand:{p}",
                              mm, "rp_tmp", f"expr__{TARGET_GENE}"))
            else:
                tasks.append(("background_peaks", at_sub, f"rand:{p}",
                              m, f"rp__{p}", f"expr__{TARGET_GENE}"))
    rows_out = run_tests(tasks, args.threads)
    res = pd.DataFrame(rows_out)
    for fam in ["celltype_specificity", "gene_specificity", "background_peaks"]:
        mask = res["family"] == fam
        res.loc[mask, "fdr"] = bh_fdr(res.loc[mask, "p"].values)
    res.to_csv(os.path.join(args.outdir, "peak_expr_tests.tsv"),
               sep="\t", index=False)

    keepcols = ["donor", "n_cells", "rna_n_cells", "acc_target", "acc_lead",
                f"expr__{TARGET_GENE}", "age", "sex", "adnc", "depth", "rna_depth"]
    merged("L6 IT", "L6 IT")[keepcols].to_csv(
        os.path.join(args.outdir, "peak_expr_L6IT_donors.tsv"), sep="\t", index=False)

    # ---------- summary ----------
    def show(fam, celltype=None):
        d = res[res["family"] == fam]
        if celltype:
            d = d[d["celltype"] == celltype]
        d = d.dropna(subset=["p"]).sort_values("p")
        print(f"\n[{fam} {celltype or ''}]")
        for _, r in d.head(12).iterrows():
            print(f"  {r['test']:<44} n={int(r['n']):>3} beta={r['beta']:+.3f} "
                  f"p={r['p']:.2e} fdr={r['fdr']:.3f} rho={r['rho']:+.3f}")
        return d

    show("celltype_specificity")
    show("gene_specificity")
    for ct, _ in CELL_PAIRS:
        d = show("background_peaks", ct)
        print(f"  background |beta| median={np.nanmedian(np.abs(d['beta'])):.3f}, "
              f"min p={d['p'].min():.2e}, frac p<0.05={(d['p'] < 0.05).mean():.2f}")
    print(f"\nL6 IT donors in main test: {len(m6)}")
    print(f"results -> {args.outdir}/peak_expr_tests.tsv")


if __name__ == "__main__":
    main()
