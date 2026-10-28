#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import gzip
import shutil
import subprocess
import concurrent.futures
from datetime import datetime
import pandas as pd

# ==========================================================
# Step 0：读取 Excel 并清洗数据
# ==========================================================
def clean_dataframe(excel_path, sheet_name):
    """
    
    读取 xlsx 并对所有字段做数据清洗，返回清洗后的 DataFrame
    
    """
    df = pd.read_excel(excel_path, sheet_name = sheet_name, dtype = str)
    df = df[df['keep'].isin(['Primary', 'Replication', 'Sensitivity'])]

    # ---------- 1. 去除所有字符串字段首尾空白 / 统一缺失值 ----------
    for col in df.columns:
        df[col] = df[col].astype(str).str.strip()
        df[col] = df[col].replace({"nan": None, "None": None, "": None})

    # ---------- 2. 状态字段标准化：YES / YESPOS / YESNEG -> 是否纳入 ----------
    df["status"] = df["status"].fillna("")
    df["include_flag"] = df["status"].str.startswith("YES")

    # 效应方向标记：POS(正向) / NEG(反向) / 无
    def parse_direction(s):
        if s.endswith("POS"):
            return "POS"
        if s.endswith("NEG"):
            return "NEG"
        return "NA"

    df["effect_direction"] = df["status"].apply(parse_direction)

    # ---------- 3. 数值字段类型转换（错误值转为 NaN） ----------
    numeric_cols = [
        "nsample", "ncase", "ncontrol", "pmid", "year",
        "score_phenotype", "score_ancestry", "score_sample",
        "score_completeness", "score_overlap", "score_recognition",
    ]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors = "coerce")

    # ncase / ncontrol / nsample 转整数（允许缺失）
    for col in ["nsample", "ncase", "ncontrol"]:
        df[col] = df[col].astype("Int64")

    # ---------- 4. 重新计算 total_score（原表中是公式字符串，需用真实数值重算） ----------
    weights = {
        "score_phenotype": 5,
        "score_ancestry": 5,
        "score_sample": 4,
        "score_completeness": 4,
        "score_overlap": 3,
        "score_recognition": 3,
    }
    weight_sum = sum(weights.values())  # =24

    def recompute_score(row):
        try:
            total = sum(row[col] * w for col, w in weights.items())
            return round(total / weight_sum, 4)
        except Exception:
            return None

    df["total_score_calc"] = df.apply(recompute_score, axis = 1)

    # ---------- 5. EAF 字段标准化：YES / EAF / 空 -> 布尔 ----------
    df["has_eaf"] = df["EAF"].isin(["YES", "EAF"])

    # ---------- 6. 疾病 / trait 名称标准化，用于生成文件夹名 ----------
    def normalize_name(name):
        if name is None:
            return None
        name = re.sub(r"[^\w\-]+", "_", name.strip())
        name = re.sub(r"_+", "_", name).strip("_")
        return name

    df["trait_key"] = df["disease"].apply(normalize_name)

    # ---------- 7. 剔除关键字段缺失或未被纳入的记录 ----------
    before = len(df)
    df = df[df["include_flag"]]
    df = df.dropna(subset = ["gwas_id", "data", "ncase", "ncontrol"])
    after_filter = len(df)

    # ---------- 8. 按 gwas_id 去重（保留第一条） ----------
    df = df.drop_duplicates(subset = ["gwas_id"], keep = "first")
    after_dedup = len(df)

    print(f"[清洗统计] 原始行数 = {before}, 有效纳入行数 = {after_filter}, "
          f"去重后行数 = {after_dedup}")

    return df.reset_index(drop = True)


def save_cleaned_table(df, out_path):
    """ 保存清洗后的数据表，便于核查 """
    df.to_csv(out_path, index = False, encoding = "utf-8-sig")
    print(f"[清洗结果已保存] {out_path}")


# ==========================================================
# Step 1：Convert & QC（对应原 bash 中的 run_convert_qc）
# ==========================================================
def run_convert_qc(python_bin, sumstats_py, data_dir, out_root,
                    trait_key, gwas_id, subdir, ncase, ncontrol):
    """对单个 sumstats 文件执行 csv -> zscore -> qc -> neff -> gzip"""

    outdir = os.path.join(out_root, subdir)
    os.makedirs(outdir, exist_ok = True)

    src_file = os.path.join(data_dir, f"{gwas_id}.tsv.gz")
    csv_out = os.path.join(outdir, f"{gwas_id}.csv")
    nomhc_out = os.path.join(outdir, f"{gwas_id}.noMHC.csv")
    gz_path = nomhc_out + ".gz"

    # ---- 若最终结果文件已存在，直接跳过（断点续传）----
    if os.path.exists(gz_path):
        print(f"[跳过] 已存在结果文件: {gz_path}")
        return

    if not os.path.exists(src_file):
        print(f"[跳过] 源文件不存在: {src_file}")
        return

    # ---- csv 转换 ----
    cmd_csv = [
        python_bin, sumstats_py, "csv",
        "--sumstats", src_file,
        "--out", csv_out,
        "--force", "--auto", "--head", "5",
        "--a1", "EA", "--a2", "OA",
        "--ncase-val", str(ncase),
        "--ncontrol-val", str(ncontrol),
    ]
    subprocess.run(cmd_csv, check = True)

    # ---- zscore | qc | neff 管道 ----
    cmd_zscore = [python_bin, sumstats_py, "zscore", "--sumstats", csv_out]
    cmd_qc = [
        python_bin, sumstats_py, "qc",
        "--exclude-ranges", "6:26000000-34000000",
        "--max-or", "10000.0",
    ]
    cmd_neff = [
        python_bin, sumstats_py, "neff",
        "--drop", "--factor", "4",
        "--out", nomhc_out, "--force",
    ]

    p1 = subprocess.Popen(cmd_zscore, stdout = subprocess.PIPE)
    p2 = subprocess.Popen(cmd_qc, stdin = p1.stdout, stdout = subprocess.PIPE)
    p1.stdout.close()
    subprocess.run(cmd_neff, stdin = p2.stdout, check = True)
    p2.stdout.close()

    # ---- gzip 压缩 ----
    if os.path.exists(nomhc_out):
        with open(nomhc_out, "rb") as f_in, gzip.open(gz_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(nomhc_out)

    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] 完成 Convert+QC: {gwas_id} ({trait_key})")


# ==========================================================
# 主流程
# ==========================================================
def main(excel_path, sheet_name, python_bin, sumstats_py, data_dir, out_root, exposure_gwas_ids, max_jobs):

    os.makedirs(out_root, exist_ok = True)

    # Step 0：清洗数据
    df = clean_dataframe(excel_path, sheet_name)

    # 拆分暴露表型与其它结局表型
    exposure_rows = df[df["gwas_id"].isin(exposure_gwas_ids)]
    outcome_rows = df[~df["gwas_id"].isin(exposure_gwas_ids)]

    tasks = []

    # 暴露表型只需转换一次（统一放到 _exposure_global 目录）
    for _, row in exposure_rows.iterrows():
        tasks.append(dict(
            trait_key = "_exposure_global",
            gwas_id = row["gwas_id"],
            subdir = "_exposure_global",
            ncase = int(row["ncase"]),
            ncontrol = int(row["ncontrol"]),
        ))

    # 其它疾病结局表型：按 trait_key 分子目录
    for _, row in outcome_rows.iterrows():
        tasks.append(dict(
            trait_key = row["trait_key"],
            gwas_id = row["gwas_id"],
            subdir = row["trait_key"],
            ncase = int(row["ncase"]),
            ncontrol = int(row["ncontrol"]),
        ))

    print(f"[任务统计] 共 {len(tasks)} 个 GWAS 数据集待处理，"
          f"并发数 = {max_jobs}")

    # Step 1：并发执行 Convert & QC
    with concurrent.futures.ThreadPoolExecutor(max_workers = max_jobs) as executor:
        futures = [
            executor.submit(
                run_convert_qc, python_bin, sumstats_py, data_dir, out_root,
                t["trait_key"], t["gwas_id"], t["subdir"], t["ncase"], t["ncontrol"]
            )
            for t in tasks
        ]
        for future in concurrent.futures.as_completed(futures):
            future.result()

    print("==========================================")
    print("全部 Convert & QC 任务已完成。")
    print("==========================================")


if __name__ == '__main__':
    
    BASE_DIR = '/home/xab/influenza_mental_genomics'
    BASE_DIR = BASE_DIR if os.path.exists(BASE_DIR) else '/home/hello/influenza_mental_genomics'
    EXCEL_PATH = f"{BASE_DIR}/Merge.infos.M.xlsx"
    SHEET_NAME = "Sheet1"

    PYTHON_BIN = "/tools/Python-3.7.9/python"
    SUMSTATS_PY = "/tools/mixer/python_convert/sumstats.py"

    DATA_DIR = f"{BASE_DIR}/data"
    OUT_ROOT = f"{BASE_DIR}/results/results.temp/mixer"

    EXPOSURE_GWAS_IDs = ["J10_INFLUENZA", "GCST90432107"]
    MAX_JOBS = 16

    main(
        EXCEL_PATH, SHEET_NAME,
        PYTHON_BIN, SUMSTATS_PY, DATA_DIR, OUT_ROOT,
        EXPOSURE_GWAS_IDs, MAX_JOBS,
    )