# -*- coding: utf-8 -*-
"""
Created on Sun Jul 19 01:03:19 2026
 @author Abao Xing
 @email  albertxn7@gmail.com
 This scripts writen by Abao Xing

   ┏┓　　┏┓
  ┏┛┻━━━━┛┻┓
  ┃　　　　  ┃
  ┃　━　　━　 ┃
  ┃　┳┛　┗┳　 ┃
  ┃　　　　　 ┃
  ┃　　　┻　　┃
  ┃　　　　　 ┃
  ┗━━┓　　　┏━┛
  　　┃　　 ┃ 神兽保佑
  　　┃　　 ┃ 代码无BUG！！！
  　　┃　　 ┗━━━━━┓
  　　┃　　　　　　  ┣┓
 　　┃　　　　　　  ┏┛┃
 　　┗┓┓┏━━━━━┳┓┏━━┛
  　　┃┫┫　   ┃┫┫
  　　┗┻┛　   ┗┻┛

"""


import os
import sys
import subprocess
import pandas as pd


# ==== 特例映射：对应 MiXeR-data-convert.py 里的特殊输出规则 ====
# Influenza 是暴露表型（exposure），转换脚本统一把它放到 _exposure_global 目录下，
# 而不是像其他疾病那样各自建 <safe_disease> 子目录。
DISEASE_DIR_OVERRIDES = {
    "influenza": "_exposure_global",
}

NO_MHC_DISEASES = {
    "influenza",
}

VALID_KEEP = ("Primary", "Replication", "Sensitivity")


def load_table(xlsx_path, sheet_name):
    """
    读取性状清单表：只保留 keep 列取值为 Primary / Replication / Sensitivity 的行。
    单变量 fit1 不剔除暴露性状（Influenza 自己也需要 fit1 结果供 fit2/test2 使用）。
    """
    df = pd.read_excel(xlsx_path, sheet_name = sheet_name)
    df.columns = [str(c).strip() for c in df.columns]

    missing = [c for c in ("disease", "gwas_id", "data", "keep") if c not in df.columns]
    if missing:
        sys.exit(f"[ERROR] {os.path.basename(xlsx_path)} [{sheet_name}] 缺少必需列: {missing}")

    for col in ("disease", "gwas_id", "data", "group"):
        df[col] = df[col].astype(str).str.strip()

    n_all = len(df)
    df = df[df["group"].isin(VALID_KEEP)]
    print(f"[INFO] 读取 {n_all} 行，keep 列有效行 {len(df)} 行（过滤掉 {n_all - len(df)} 行）")
    return df


def pick_datasets(df, disease_name, keep):
    """挑选某个疾病名下、指定 keep 的数据集（大小写不敏感匹配疾病名）"""
    sub = df[df["disease"].str.lower() == disease_name.strip().lower()]
    if keep != "ALL":
        sub = sub[sub["group"] == keep]
    return sub


def build_id_from_data(data_filename):
    """
    data 列形如 GCST90449053.tsv.gz -> 返回 GCST90449053
    实际读取路径为 <data_dir>/<disease_subdir>/<ID>.noMHC.csv.gz
    """
    base = data_filename
    for ext in (".tsv.gz", ".csv.gz", ".gz"):
        if base.endswith(ext):
            base = base[: -len(ext)]
            break
    return base


def safe_name(name):
    return name.replace(" ", "_").replace("/", "_")


def get_data_subdir(disease):
    """
    返回该疾病在 data_dir 下实际对应的子目录名。
    默认使用 safe_name(disease)；若在 DISEASE_DIR_OVERRIDES 中有特例（如 Influenza），
    则使用覆盖值（如 _exposure_global）。
    """
    key = disease.strip().lower()
    return DISEASE_DIR_OVERRIDES.get(key, safe_name(disease))


def get_trait_file_suffix(disease):
    """
    返回该疾病对应数据文件的后缀。
    当前 MiXeR-data-convert.py 对所有性状都输出 .noMHC.csv.gz，
    故两个分支返回值相同；保留 NO_MHC_DISEASES 分支以便将来区分处理。
    """
    key = disease.strip().lower()
    return ".noMHC.csv.gz" if key in NO_MHC_DISEASES else ".noMHC.csv.gz"


def build_header_lines(base_dir, ldsr_dir, data_dir, threads, use_docker, docker_image,
                        python_bin, mixer_py, lib):
    """
    生成脚本 header：BASE_DIR / REF / DATA_DIR / THREADS / LIB
    以及 PYTHON / MIXER 变量，兼容 docker 与直连 tools 两种方式。
    """
    lines = []
    lines.append("#!/bin/bash")
    lines.append("set -euo pipefail")
    lines.append("")
    lines.append(f'BASE_DIR="{base_dir}"')
    lines.append(f'REF="{ldsr_dir}"')
    lines.append(f'DATA_DIR="{data_dir}"')
    lines.append(f'THREADS="{threads}"')
    lines.append(f'LIB="{lib}"')
    lines.append("")

    if use_docker:
        lines.append(f'export DOCKER_RUN="docker run -v {base_dir}:{base_dir} -v /references:/references -w {base_dir}"')
        lines.append(f'export MIXER_PY="$DOCKER_RUN {docker_image} python /tools/mixer/precimed/mixer.py"')
        lines.append('PYTHON=""')
        lines.append('MIXER="${MIXER_PY}"')
        lines.append("")
    else:
        lines.append(f'PYTHON="{python_bin}"')
        lines.append(f'MIXER="{mixer_py}"')
        lines.append("")

    return lines


def build_fit1_block(disease, gwas_id, trait_id, n_rep):
    """
    生成 fit1 命令块。

    【断点续跑】
      - 输出 ${ID}.fit.rep${rep}.json 已存在且非空 -> 跳过该 rep
      - --extract 的 rep${rep}.snps 不存在 -> 警告并跳过该 rep（不中断整个脚本）
      - trait 输入文件不存在 -> 警告并结束该数据集（exit 0，避免 set -e 掐断合并脚本）
    """
    safe_disease = safe_name(disease)
    th_dir = f"${{BASE_DIR}}/results/results/mixer/fit1/{safe_disease}_{gwas_id}"

    # 疾病对应的数据子目录 / 文件后缀分别通过 get_data_subdir() / get_trait_file_suffix() 获取，
    # 兼容 Influenza（_exposure_global 目录）这类特例。
    data_subdir = get_data_subdir(disease)
    file_suffix = get_trait_file_suffix(disease)

    lines = []
    lines.append(f"# ==== {disease} ({gwas_id}) : fit1 ====")
    lines.append(f'ID="{trait_id}"')
    lines.append(f'THDIR="{th_dir}"')
    lines.append(f'TRAIT_FILE="${{DATA_DIR}}/{data_subdir}/${{ID}}{file_suffix}"')
    lines.append('mkdir -p "${THDIR}"')

    # ---- 输入检查：缺少 GWAS 输入文件时不必进入循环 ----
    lines.append('if [ ! -s "${TRAIT_FILE}" ]; then')
    lines.append(f'  echo "[WARN] trait file not found: ${{TRAIT_FILE}} , skipping {disease} ({gwas_id})." >&2')
    lines.append('  exit 0')
    lines.append('fi')

    lines.append(f"for rep in $(seq 1 {n_rep})")
    lines.append("do")
    lines.append(f'  echo "[{disease} | ${{ID}}] Processing fit1 replicate $rep"')

    # ---- 断点续跑：该 rep 结果已存在则跳过 ----
    lines.append('  OUT_JSON="${THDIR}/${ID}.fit.rep${rep}.json"')
    lines.append('  if [ -s "${OUT_JSON}" ]; then')
    lines.append('    echo "  -> fit1 rep $rep already exists, skipping."')
    lines.append('    continue')
    lines.append('  fi')

    # ---- 依赖检查：随机 SNP 子集文件必须存在 ----
    lines.append('  EXTRACT_FILE="${REF}/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${rep}.snps"')
    lines.append('  if [ ! -s "${EXTRACT_FILE}" ]; then')
    lines.append('    echo "  -> [WARN] extract snps not found: ${EXTRACT_FILE} , skipping rep $rep." >&2')
    lines.append('    continue')
    lines.append('  fi')

    bim_file = "${REF}/1000G.EUR.QC.@.bim"
    ld_file = "${REF}/1000G.EUR.QC.@.run4.ld"

    lines.append(
        '  ${PYTHON} ${MIXER} fit1 \\\n'
        '        --trait1-file "${TRAIT_FILE}" \\\n'
        '        --out "${OUT_JSON%.json}" \\\n'
        '        --extract "${EXTRACT_FILE}" \\\n'
        f'        --bim-file "{bim_file}" \\\n'
        f'        --ld-file "{ld_file}" \\\n'
        '        --lib "${LIB}" \\\n'
        '        --threads "${THREADS}"'
    )
    lines.append("done")
    lines.append("")
    return lines


def iter_selected_datasets(df, diseases, keep):
    """
    统一遍历所有 (disease, row) 组合，供合并脚本 / 拆分脚本两种模式共用。
    返回 generator，元素为 (disease, gwas_id, trait_id) 三元组。
    """
    for disease in diseases:
        sub = pick_datasets(df, disease, keep)
        if sub.empty:
            print(f"[WARN] 未找到 {disease} 在 keep={keep} 下的数据集，跳过", file = sys.stderr)
            continue

        for _, row in sub.iterrows():
            gwas_id = row.get("gwas_id", "")
            data_file = row.get("data", "")
            if not data_file or data_file.lower() == "nan":
                print(f"[WARN] {disease} ({gwas_id}) 缺少 data 文件名，跳过", file = sys.stderr)
                continue

            trait_id = build_id_from_data(data_file)
            yield disease, gwas_id, trait_id


def get_disease_list(df, exclude_influenza):
    """保持出现顺序去重；可选剔除暴露性状"""
    diseases = list(dict.fromkeys(df["disease"].dropna().tolist()))
    if exclude_influenza:
        diseases = [d for d in diseases if d.strip().lower() != "influenza"]
    return diseases


def generate_script(xlsx_path, sheet, keep, exclude_influenza,
                     base_dir, data_dir, n_rep, threads,
                     python_bin, mixer_py, lib, ldsr_dir,
                     use_docker, docker_image,
                     mode, out_script):
    """生成单个合并的 shell 脚本（所有数据集的 fit1 命令都写在一个文件里）"""

    df = load_table(xlsx_path, sheet)
    diseases = get_disease_list(df, exclude_influenza)

    lines = build_header_lines(base_dir, ldsr_dir, data_dir, threads,
                                use_docker, docker_image,
                                python_bin, mixer_py, lib)

    n_jobs = 0
    for disease, gwas_id, trait_id in iter_selected_datasets(df, diseases, keep):
        block = build_fit1_block(
            disease = disease,
            gwas_id = gwas_id,
            trait_id = trait_id,
            n_rep = n_rep,
        )
        lines.extend(block)
        n_jobs += 1

    script_content = "\n".join(lines)

    if mode == "script":
        with open(out_script, "w", encoding = "utf-8") as f:
            f.write(script_content)
        os.chmod(out_script, 0o755)
        print(f"[OK] 已生成合并脚本: {out_script}")
        print(f"     共 {n_jobs} 个数据集的 fit1 任务，每个重复 {n_rep} 次, threads={threads}")
        print("     已存在的 rep 会自动跳过（支持中断续跑）")
        print(f"     请确认 data_dir ({data_dir}) 下各疾病子目录（含 _exposure_global）中已存在对应数据文件后执行：")
        print(f"     bash {out_script}")
    else:
        print(f"[RUN] 直接执行模式，共 {n_jobs} 个数据集任务")
        proc = subprocess.run(["bash", "-c", script_content])
        sys.exit(proc.returncode)


def generate_split_scripts(xlsx_path, sheet, keep, exclude_influenza,
                            base_dir, data_dir, n_rep, threads,
                            python_bin, mixer_py, lib, ldsr_dir,
                            use_docker, docker_image,
                            shells_dir, nohup_summary_name):
    """
    一个数据集一个 shell 脚本：
      <shells_dir>/<idx>.<disease>_<gwas_id>_fit1.sh
    并打印 / 汇总每个脚本对应的 nohup 后台运行命令。
    """

    df = load_table(xlsx_path, sheet)
    diseases = get_disease_list(df, exclude_influenza)

    os.makedirs(shells_dir, exist_ok = True)

    nohup_cmds = []
    idx = 0

    for disease, gwas_id, trait_id in iter_selected_datasets(df, diseases, keep):
        idx += 1
        safe_disease = safe_name(disease)

        header = build_header_lines(base_dir, ldsr_dir, data_dir, threads,
                                     use_docker, docker_image,
                                     python_bin, mixer_py, lib)
        block = build_fit1_block(
            disease = disease,
            gwas_id = gwas_id,
            trait_id = trait_id,
            n_rep = n_rep,
        )

        script_content = "\n".join(header + block)

        script_name = f"{idx}.{safe_disease}_{gwas_id}_fit1.sh"
        script_path = os.path.join(shells_dir, script_name)

        with open(script_path, "w", encoding = "utf-8") as f:
            f.write(script_content)
        os.chmod(script_path, 0o755)

        log_name = script_name.replace(".sh", ".log")
        log_path = os.path.join(shells_dir, log_name)

        nohup_cmd = f"nohup bash {script_path} > {log_path} 2>&1 &"
        nohup_cmds.append(nohup_cmd)

    # 打印全部 nohup 命令
    print(f"[OK] 已在 {shells_dir} 下生成 {idx} 个独立 fit1 脚本，每个重复 {n_rep} 次, threads={threads}")
    print("     支持中断续跑：重复执行同一脚本只会补跑缺失的 rep")
    print("[INFO] 可直接复制以下命令逐个提交后台运行：")
    print("")
    for cmd in nohup_cmds:
        print(cmd)
    print("")

    # 同时写入一个汇总文件，方便一次性执行
    summary_path = os.path.join(shells_dir, nohup_summary_name)
    with open(summary_path, "w", encoding = "utf-8") as f:
        f.write("#!/bin/bash\n")
        f.write("\n".join(nohup_cmds) + "\n")
    os.chmod(summary_path, 0o755)
    print(f"[OK] 所有 nohup 命令已汇总写入: {summary_path}")
    print(f"     如需一次性全部提交，可执行: bash {summary_path}")


if __name__ == '__main__':

    BASE_DIR = '/home/xab/influenza_mental_genomics'
    BASE_DIR = BASE_DIR if os.path.exists(BASE_DIR) else '/home/hello/influenza_mental_genomics'

    # ==== 输入文件配置 ====
    xlsx_path = f"{BASE_DIR}/GWAS.infos.xlsx"
    sheet = "Sheet1"

    # keep 取值: "Primary" / "Replication" / "Sensitivity" / "ALL"
    keep = "Replication"

    # 是否排除 Influenza（默认不排除：暴露表型自身也需要 fit1 结果）
    exclude_influenza = False

    # ==== 目录 / 工具路径配置 ====
    # data_dir 指向 MiXeR-data-convert.py 实际输出的根目录，
    # 其下按疾病名（safe_name 后的形式）分子目录存放各自的 <ID>.noMHC.csv.gz；
    # Influenza（暴露表型）是特例，存放在 _exposure_global 子目录下，
    # 已通过 get_data_subdir() / get_trait_file_suffix() 在 build_fit1_block() 中处理。
    data_dir = f'{BASE_DIR}/results/results.temp/mixer'
    n_rep = 20
    threads = 12

    # 直连 tools 路径方式
    python_bin = "/tools/Python-3.7.9/python"
    mixer_py = "/tools/mixer/precimed/mixer.py"
    lib = "/tools/mixer/src/build/lib/libbgmg.so"
    ldsr_dir = "/references/MiXeR/1000G_EUR_Phase3_plink"

    # docker 方式
    use_docker = True
    docker_image = "ghcr.io/precimed/gsa-mixer:2.2.1"

    # 输出模式: "script" 只生成合并 shell 脚本；"run" 生成后直接执行合并脚本
    mode = "script"
    out_script = f"{BASE_DIR}/run_fit1_all.sh"

    # ==== 是否为每个数据集单独生成一个 shell 脚本 ====
    split_per_dataset = True
    shells_dir = f"/home/xab/influenza_mental_genomics/scripts/MiXeR-fit1-shells/{keep}"
    shells_dir = shells_dir if os.path.exists(os.path.dirname(os.path.dirname(shells_dir))) \
        else f"/home/hello/influenza_mental_genomics/scripts/MiXeR-fit1-shells/{keep}"
    nohup_summary_name = "run_all_fit1_nohup.sh"

    # ---- 合并脚本（原有功能，保留） ----
    generate_script(
        xlsx_path = xlsx_path,
        sheet = sheet,
        keep = keep,
        exclude_influenza = exclude_influenza,
        base_dir = BASE_DIR,
        data_dir = data_dir,
        n_rep = n_rep,
        threads = threads,
        python_bin = python_bin,
        mixer_py = mixer_py,
        lib = lib,
        ldsr_dir = ldsr_dir,
        use_docker = use_docker,
        docker_image = docker_image,
        mode = mode,
        out_script = out_script,
    )

    # ---- 拆分脚本：一个数据集一个 shell，并打印 nohup 命令 ----
    if split_per_dataset:
        generate_split_scripts(
            xlsx_path = xlsx_path,
            sheet = sheet,
            keep = keep,
            exclude_influenza = exclude_influenza,
            base_dir = BASE_DIR,
            data_dir = data_dir,
            n_rep = n_rep,
            threads = threads,
            python_bin = python_bin,
            mixer_py = mixer_py,
            lib = lib,
            ldsr_dir = ldsr_dir,
            use_docker = use_docker,
            docker_image = docker_image,
            shells_dir = shells_dir,
            nohup_summary_name = nohup_summary_name,
        )