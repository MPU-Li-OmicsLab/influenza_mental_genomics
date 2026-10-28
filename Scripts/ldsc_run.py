# -*- coding: utf-8 -*-

import re
import os
import shutil
import subprocess
import pandas as pd
from pathlib import Path
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

import logging
logging.basicConfig(level = logging.INFO, format = "[%(asctime)s] %(message)s", datefmt = "%H:%M:%S")
log = logging.getLogger(__name__)

DISEASES_DICT = {
    'Mental-Disorders': [
        'Alcohol dependence', 'Alzheimers disease', 'Anorexia nervosa', 'Anxiety disorder', 
        'Attention deficit hyperactivity disorder', 'Autism spectrum disorder', 'Bipolar disorder', 
        'Major depressive disorder', 'Obsessive-compulsive disorder', 'Opiate dependence', 'Parkinson disease', 
        'Phobic disorder', 'Post traumatic stress disorder', 'Schizophrenia', 'Sleep disorder', 
        'Suicide attempt', 'Tobacco use disorder', 'Tourette syndrome'
    ]
}

#
diseases = [disease for diseases in DISEASES_DICT.values() for disease in diseases]

# Exprosure
EXPOSURE_ID = "J10_INFLUENZA" # J10_INFLUENZA or GCST90432107

# Path configuration
BASE_DIR   = '/home/hello/influenza_mental_genomics'
INFO_FILE  = f'{BASE_DIR}/GWAS.infos.xlsx'
DATA_DIR   = f'{BASE_DIR}/data'
PYTHON2    = '/tools/Python-2.7.16/python'
LDSC_PY    = '/tools/ldsc/ldsc.py'
MUNGE_PY   = '/tools/ldsc/munge_sumstats.py'
HAPMAP     = f'{BASE_DIR}/references/w_hm3.snplist'
LD_EUR     = f'{BASE_DIR}/references/eur_w_ld_chr/'
LD_EAS     = f'{BASE_DIR}/references/eas_w_ld_chr/'
OUT_DIR    = f'{BASE_DIR}/results/results/ldsc'
TMP_DIR    = f'{BASE_DIR}/results/results.temp/ldsc'
LOG_FILE   = f'{BASE_DIR}/results/logs/run.ldsc.log'

# Concurrency configuration: None means automatic.
MAX_WORKERS = 8

# Column name mapping
COL_MAP = {
    "snp": "SNP",
    "a1": "EA",
    "a2": "OA",
    "p": "PVAL",
    "frq": "EAF",
    "signed": "BETA",
}


def _is_continuous(row):
    """
    
    Determine whether the phenotype is continuous:
    + If both ncase and ncontrol are empty, and nsample has a specific value assigned → continuous phenotype
    + Otherwise, it is treated as a binary (case-control) phenotype
    """
    ncase_empty = pd.isna(row.get("ncase"))
    ncontrol_empty = pd.isna(row.get("ncontrol"))
    nsample_valid = pd.notna(row.get("nsample"))
    return ncase_empty and ncontrol_empty and nsample_valid


def load_traits(df):

    # Determine phenotype type (binary / continuous)
    df["trait_type"] = df.apply(lambda row: "continuous" if _is_continuous(row) else "binary", axis = 1)

    # Calculate effective sample size (n_eff formula for case-control, nsample directly for continuous)
    def calc_n_eff(row):
        if row["trait_type"] == "continuous":
            return int(round(row["nsample"]))
        else:
            return int(round(4 * row["ncase"] * row["ncontrol"] / (row["ncase"] + row["ncontrol"])))

    df["n_eff"] = df.apply(calc_n_eff, axis = 1)

    # Calculate samp_prev (only needed for binary phenotypes; None for continuous)
    def calc_samp_prev(row):
        if row["trait_type"] == "continuous":
            return None
        else:
            return round(row["ncase"] / row["nsample"], 6)

    df["samp_prev"] = df.apply(calc_samp_prev, axis = 1)

    # Population
    def fix_ancestry(row):
        pop = str(row["population"]).upper()
        if "EUROPEAN" in pop:
            return "EUR"
        elif "EAST ASIAN" in pop or "EAS" in pop:
            return "EAS"
        else:
            return "EUR"

    df["ancestry"] = df.apply(fix_ancestry, axis=1)

    # Construct file path
    def get_filepath(row):
        fname = str(row["data"]) if pd.notna(row["data"]) else f"{row['gwas_id']}.tsv.gz"
        if not fname.endswith(".gz"): fname += ".gz" if not fname.endswith(".tsv") else ""
        return f'{DATA_DIR}/{fname}'
    

    df["filepath"] = df.apply(get_filepath, axis = 1)

    # Convert to list of dictionaries
    traits = []
    for _, row in df.iterrows():
        traits.append({
            "id": row["gwas_id"],
            "filepath": row["filepath"],
            "n_eff": int(row["n_eff"]),
            "samp_prev": None if row["samp_prev"] is None else float(row["samp_prev"]),
            "pop_prev": None,
            "ancestry": row["ancestry"],
            "trait_type": row["trait_type"],
        })
        if verbose: log.info(f"  Load {row['gwas_id']}: Type = {row['trait_type']}, N_eff = {row['n_eff']}, "
                 f"ancestry = {row['ancestry']}, samp_prev = {row['samp_prev']}")

    return traits

def ld_path(ancestry):
    return LD_EUR if ancestry.upper() == "EUR" else LD_EAS

def run(cmd, logfile, mode = "a"):
    if verbose: log.info(f"  LOG → {logfile}")
    Path(logfile).parent.mkdir(parents = True, exist_ok = True)
    with open(logfile, mode, encoding = "utf-8") as fh:
        fh.write(f"\n\n{'='*60}\n")
        fh.write(f"CMD: {' '.join(cmd)}\n")
        fh.write(f"TIME: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        fh.write(f"{'='*60}\n\n")
        ret = subprocess.run(cmd, stdout = fh, stderr = subprocess.STDOUT, text = True)
    if ret.returncode != 0:
        raise RuntimeError(f"Command failed, please check {logfile}")   

def _get_max_workers():
    if MAX_WORKERS is not None:
        return MAX_WORKERS
    return min(multiprocessing.cpu_count(), 8)


def _munge_worker(args):
    """
    
    Munge task for a single trait
    """
    t, dislabel, verbose_local = args
    try:
        cmd = [
            PYTHON2, MUNGE_PY,
            "--sumstats", t["filepath"],
            "--merge-alleles", HAPMAP,
            "--out", f'{TMP_DIR}/{dislabel}/{t["id"]}',
            "--chunksize", "500000",
            "--snp", COL_MAP["snp"],
            "--a1", COL_MAP["a1"],
            "--a2", COL_MAP["a2"],
            "--p", COL_MAP["p"],
            "--signed-sumstats", f"{COL_MAP['signed']},0",
            "--N", str(t["n_eff"]),
            "--frq", COL_MAP["frq"],
        ]

        otfile = f'{TMP_DIR}/{dislabel}/{t["id"]}.sumstats.gz'
        if not os.path.exists(otfile):
            if verbose_local:
                run(cmd, logfile=LOG_FILE)
            else:
                print(' '.join(cmd))
            return ("ok", t["id"], "munge complete")
        else:
            return ("skip", t["id"], "munge output already exists, skipping")
    except Exception as e:
        return ("error", t["id"], f"munge failed: {e}")


def _h2_worker(args):
    """
    
    Run LDSC heritability (h2) estimation for a single trait.
    """
    t, dislabel, verbose_local = args
    try:
        munged = f'{TMP_DIR}/{dislabel}/{t["id"]}.sumstats.gz'
        ld = ld_path(t["ancestry"])

        cmd = [
            PYTHON2, LDSC_PY,
            "--h2", str(munged),
            "--ref-ld-chr", ld,
            "--w-ld-chr", ld,
            "--out", f'{OUT_DIR}/{dislabel}/h2/{t["id"]}',
        ]

        if t["trait_type"] == "binary" and t["samp_prev"] is not None and t["pop_prev"] is not None:
            cmd += ["--samp-prev", str(t["samp_prev"]), "--pop-prev", str(t["pop_prev"])]

        otfile = f'{OUT_DIR}/{dislabel}/h2/{t["id"]}.log'
        if not os.path.exists(otfile):
            if verbose_local:
                run(cmd, logfile=LOG_FILE)
            else:
                print(' '.join(cmd))
            return ("ok", t["id"], "h2 complete")
        else:
            return ("skip", t["id"], "h2 output already exists, skipping")
    except Exception as e:
        return ("error", t["id"], f"h2 failed: {e}")


def _rg_worker(args):
    """
    
    Compute genetic correlation (rg) between a pair of traits using LDSC.
    """
    exposure, t2, dislabel, verbose_local = args
    try:
        m1 = f'{TMP_DIR}/{dislabel}/{exposure["id"]}.sumstats.gz'
        m2 = f'{TMP_DIR}/{dislabel}/{t2["id"]}.sumstats.gz'

        if exposure["ancestry"] != t2["ancestry"]:
            ld = ld_path(exposure["ancestry"])
        else:
            ld = ld_path(exposure["ancestry"])

        cmd = [
            PYTHON2, LDSC_PY,
            "--rg", f"{m1},{m2}",
            "--ref-ld-chr", ld,
            "--w-ld-chr", ld,
            "--out", f'{OUT_DIR}/{dislabel}/{t2["id"]}',
        ]

        otfile = f'{OUT_DIR}/{dislabel}/{t2["id"]}.log'
        if not os.path.exists(otfile):
            if verbose_local:
                run(cmd, logfile=LOG_FILE)
            else:
                print(' '.join(cmd))
            return ("ok", f"{exposure['id']} × {t2['id']}", "rg complete")
        else:
            return ("skip", f"{exposure['id']} × {t2['id']}", "rg output already exists, skipping")
    except Exception as e:
        return ("error", f"{exposure['id']} × {t2['id']}", f"rg failed: {e}")


# ============================================================================
# Step functions (using ProcessPoolExecutor internally for concurrency)
# ============================================================================

def step_munge(dislabel, traits):
    if verbose: log.info("=" * 50)
    if verbose: log.info("Step 1: munge_sumstats (format standardization) [concurrent mode]")
    if verbose: log.info("=" * 50)

    max_workers = _get_max_workers()
    if verbose: log.info(f"  Concurrency: {max_workers}")

    tasks = [(t, dislabel, verbose) for t in traits]

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(_munge_worker, task): task[0]["id"] for task in tasks}
        for future in as_completed(futures):
            status, tid, msg = future.result()
            if status == "error":
                log.error(f"  ✗ {tid}: {msg}")
            elif verbose:
                log.info(f"  ✓ {tid}: {msg}")


def step_h2(dislabel, traits):
    if verbose: log.info("=" * 50)
    if verbose: log.info("Step 2: h2 (heritability estimation) [concurrent mode]")
    if verbose: log.info("=" * 50)

    max_workers = _get_max_workers()
    if verbose: log.info(f"  Concurrency: {max_workers}")

    tasks = [(t, dislabel, verbose) for t in traits]

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(_h2_worker, task): task[0]["id"] for task in tasks}
        for future in as_completed(futures):
            status, tid, msg = future.result()
            if status == "error":
                log.error(f"  ✗ {tid}: {msg}")
            elif verbose:
                log.info(f"  ✓ {tid}: {msg}")


def step_rg(dislabel, traits):
    if verbose: log.info("=" * 50)
    if verbose: log.info("Step 3: rg (genetic correlation estimation) [concurrent mode]")
    if verbose: log.info("=" * 50)

    exposure = next(t for t in traits if t["id"] == EXPOSURE_ID)
    max_workers = _get_max_workers()
    if verbose: log.info(f"  Concurrency: {max_workers}")

    tasks = [(exposure, t2, dislabel, verbose) for t2 in traits if t2["id"] != EXPOSURE_ID]

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(_rg_worker, task): task[1]["id"] for task in tasks}
        for future in as_completed(futures):
            status, pair, msg = future.result()
            if status == "error":
                log.error(f"  ✗ {pair}: {msg}")
            elif verbose:
                log.info(f"  ✓ {pair}: {msg}")


if __name__ == '__main__':
    
    # Enable logging
    verbose = True
    
    # Whether to force re-run
    FORCE = False
    
    # Read information table
    GWASDF = pd.read_excel(INFO_FILE, sheet_name = 'Sheet1')
    
    # Iterate and run LDSC for each disease
    for disease in diseases:
    
        df = GWASDF[((GWASDF['disease'] == disease) & (GWASDF['population'] == 'European')) | (GWASDF['disease'] == 'Influenza')]
        
        dislabel = re.sub('\\s', '_', disease).lower()
        
        if FORCE and os.path.exists(f'{OUT_DIR}/{dislabel}'): shutil.rmtree(f'{OUT_DIR}/{dislabel}')
        if FORCE and os.path.exists(f'{TMP_DIR}/{dislabel}'): shutil.rmtree(f'{TMP_DIR}/{dislabel}')
        if not os.path.exists(f'{OUT_DIR}/{dislabel}/h2'): os.makedirs(f'{OUT_DIR}/{dislabel}/h2')
        if not os.path.exists(f'{TMP_DIR}/{dislabel}'): os.makedirs(f'{TMP_DIR}/{dislabel}')
            
        traits = load_traits(df)
        step_munge(dislabel, traits); print()
        step_h2(dislabel, traits); print()
        step_rg(dislabel, traits); print()
