# -*- coding: utf-8 -*-
"""
Created on Thu Jul 24 22:08:00 2026
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


DISEASE_DIR_OVERRIDES = {
    "influenza": "_exposure_global",
}

NO_MHC_DISEASES = {
    "influenza",
}


def load_table(xlsx_path, sheet_name):
    df = pd.read_excel(xlsx_path, sheet_name = sheet_name)
    df = df[~df['gwas_id'].isin([exposure_trait_id])]
    df.columns = [str(c).strip() for c in df.columns]
    for col in ("disease", "gwas_id", "data", "group"):
        if col in df.columns:
            df[col] = df[col].astype(str).str.strip()
    df = df[df['group'].isin(['Primary', 'Replication', 'Sensitivity'])]
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
    这里只取 ID 前缀，实际读取路径为 <data_dir>/<disease_subdir>/<ID>[.noMHC].csv.gz。
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
    """返回该疾病在 data_dir 下实际对应的子目录名"""
    key = disease.strip().lower()
    return DISEASE_DIR_OVERRIDES.get(key, safe_name(disease))


def get_trait_file_suffix(disease):
    """返回该疾病对应数据文件的后缀"""
    key = disease.strip().lower()
    return ".noMHC.csv.gz" if key in NO_MHC_DISEASES else ".noMHC.csv.gz"


def build_header_lines(base_dir, ldsr_dir, data_dir, threads, use_docker, docker_image,
                        python_bin, mixer_py, lib):
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


def build_fit2_block(exposure_disease, exposure_gwas_id, exposure_trait_id,
                      outcome_disease, outcome_gwas_id, outcome_trait_id,
                      n_rep):
    """生成 fit2 命令块，输出目录与 fit1 保持一致的扁平结构"""
    safe_exposure = safe_name(exposure_disease)
    safe_outcome = safe_name(outcome_disease)

    exposure_data_subdir = get_data_subdir(exposure_disease)
    outcome_data_subdir = get_data_subdir(outcome_disease)
    exposure_file_suffix = get_trait_file_suffix(exposure_disease)
    outcome_file_suffix = get_trait_file_suffix(outcome_disease)

    pair_name = f"{safe_exposure}({exposure_trait_id})_vs_{safe_outcome}({outcome_trait_id})"
    th_dir = f"${{BASE_DIR}}/results/results/mixer/fit2/{pair_name}"

    lines = []
    lines.append(f"# ==== {exposure_disease} vs {outcome_disease} ({exposure_gwas_id} vs {outcome_gwas_id}) : fit2 ====")
    lines.append(f'ID1="{exposure_trait_id}"')
    lines.append(f'ID2="{outcome_trait_id}"')
    lines.append(f'PAIR="{pair_name}"')
    lines.append(f'THDIR="{th_dir}"')
    lines.append('mkdir -p "${THDIR}"')
    lines.append(f"for rep in $(seq 1 {n_rep})")
    lines.append("do")
    lines.append(f'  echo "[{exposure_disease} vs {outcome_disease} | ${{PAIR}}] Processing fit2 replicate $rep"')

    lines.append('  OUT_JSON="${THDIR}/${ID1}_vs_${ID2}.rep${rep}.json"')
    lines.append('  if [ -f "${OUT_JSON}" ]; then')
    lines.append('    echo "  -> fit2 rep $rep already exists, skipping."')
    lines.append('    continue')
    lines.append('  fi')

    trait1_file = f"${{DATA_DIR}}/{exposure_data_subdir}/${{ID1}}{exposure_file_suffix}"
    trait2_file = f"${{DATA_DIR}}/{outcome_data_subdir}/${{ID2}}{outcome_file_suffix}"

    trait1_params = f"${{BASE_DIR}}/results/results/mixer/fit1/{safe_exposure}_{exposure_gwas_id}/${{ID1}}.fit.rep${{rep}}.json"
    trait2_params = f"${{BASE_DIR}}/results/results/mixer/fit1/{safe_outcome}_{outcome_gwas_id}/${{ID2}}.fit.rep${{rep}}.json"

    outdir = "${THDIR}"
    extract_file = "${REF}/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${rep}.snps"
    bim_file = "${REF}/1000G.EUR.QC.@.bim"
    ld_file = "${REF}/1000G.EUR.QC.@.run4.ld"

    cmd = (
        '  ${PYTHON} ${MIXER} fit2 \\\n'
        f'        --trait1-file "{trait1_file}" \\\n'
        f'        --trait2-file "{trait2_file}" \\\n'
        f'        --trait1-params-file "{trait1_params}" \\\n'
        f'        --trait2-params-file "{trait2_params}" \\\n'
        '        --out "${OUT_JSON%.json}" \\\n'
        f'        --extract "{extract_file}" \\\n'
        f'        --bim-file "{bim_file}" \\\n'
        f'        --ld-file "{ld_file}" \\\n'
        '        --lib "${LIB}" \\\n'
        '        --threads "${THREADS}"'
    )
    lines.append(cmd)
    lines.append("done")
    lines.append("")
    return lines


def iter_selected_datasets(df, diseases, keep):
    """统一遍历所有 (disease, row) 组合"""
    for disease in diseases:
        sub = pick_datasets(df, disease, keep)
        if sub.empty:
            print(f"[WARN] 未找到 {disease} 在 keep={keep} 下的数据集，跳过", file=sys.stderr)
            continue

        for _, row in sub.iterrows():
            gwas_id = row.get("gwas_id", "")
            data_file = row.get("data", "")
            if not data_file or data_file.lower() == "nan":
                print(f"[WARN] {disease} ({gwas_id}) 缺少 data 文件名，跳过", file=sys.stderr)
                continue

            trait_id = build_id_from_data(data_file)
            yield disease, gwas_id, trait_id


def generate_script(xlsx_path, sheet, keep, exclude_influenza,
                     exposure_disease, exposure_gwas_id, exposure_trait_id,
                     base_dir, data_dir, n_rep, threads,
                     python_bin, mixer_py, lib, ldsr_dir,
                     use_docker, docker_image,
                     mode, out_script):
    """生成单个合并的 shell 脚本"""

    df = load_table(xlsx_path, sheet)

    diseases = list(dict.fromkeys(df["disease"].dropna().tolist()))
    diseases = [d for d in diseases if d.strip().lower() != exposure_disease.strip().lower()]

    lines = build_header_lines(base_dir, ldsr_dir, data_dir, threads,
                                use_docker, docker_image,
                                python_bin, mixer_py, lib)

    n_jobs = 0
    for outcome_disease, outcome_gwas_id, outcome_trait_id in iter_selected_datasets(df, diseases, keep):
        block = build_fit2_block(
            exposure_disease = exposure_disease,
            exposure_gwas_id = exposure_gwas_id,
            exposure_trait_id = exposure_trait_id,
            outcome_disease = outcome_disease,
            outcome_gwas_id = outcome_gwas_id,
            outcome_trait_id = outcome_trait_id,
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
        print(f"     暴露性状: {exposure_disease} ({exposure_gwas_id})")
        print(f"     共 {n_jobs} 个结局性状的 fit2 任务，每个重复 {n_rep} 次, threads={threads}")
        print(f"     bash {out_script}")
    else:
        print(f"[RUN] 直接执行模式，共 {n_jobs} 个数据集任务")
        proc = subprocess.run(["bash", "-c", script_content])
        sys.exit(proc.returncode)


def generate_split_scripts(xlsx_path, sheet, keep, exclude_influenza,
                            exposure_disease, exposure_gwas_id, exposure_trait_id,
                            base_dir, data_dir, n_rep, threads,
                            python_bin, mixer_py, lib, ldsr_dir,
                            use_docker, docker_image,
                            shells_dir, nohup_summary_name):
    """一个结局性状一个 shell 脚本，并打印 / 汇总 nohup 命令"""

    df = load_table(xlsx_path, sheet)

    diseases = list(dict.fromkeys(df["disease"].dropna().tolist()))
    diseases = [d for d in diseases if d.strip().lower() != exposure_disease.strip().lower()]

    os.makedirs(shells_dir, exist_ok = True)

    nohup_cmds = []
    idx = 0

    for outcome_disease, outcome_gwas_id, outcome_trait_id in iter_selected_datasets(df, diseases, keep):
        idx += 1
        safe_exposure = safe_name(exposure_disease)
        safe_outcome = safe_name(outcome_disease)

        header = build_header_lines(base_dir, ldsr_dir, data_dir, threads,
                                     use_docker, docker_image,
                                     python_bin, mixer_py, lib)
        block = build_fit2_block(
            exposure_disease = exposure_disease,
            exposure_gwas_id = exposure_gwas_id,
            exposure_trait_id = exposure_trait_id,
            outcome_disease = outcome_disease,
            outcome_gwas_id = outcome_gwas_id,
            outcome_trait_id = outcome_trait_id,
            n_rep = n_rep,
        )

        script_content = "\n".join(header + block)

        script_name = f"{idx}.{safe_exposure}_vs_{safe_outcome}_{outcome_gwas_id}_fit2.sh"
        script_path = os.path.join(shells_dir, script_name)

        with open(script_path, "w", encoding = "utf-8") as f:
            f.write(script_content)
        os.chmod(script_path, 0o755)

        log_name = script_name.replace(".sh", ".log")
        log_path = os.path.join(shells_dir, log_name)

        nohup_cmd = f"nohup bash {script_path} > {log_path} 2>&1 &"
        nohup_cmds.append(nohup_cmd)

    print(f"[OK] 已在 {shells_dir} 下生成 {idx} 个独立 fit2 脚本")
    print(f"     暴露性状: {exposure_disease} ({exposure_gwas_id})")
    print(f"     每个任务重复 {n_rep} 次, threads={threads}")
    print("[INFO] 可直接复制以下命令逐个提交后台运行：")
    print("")
    for cmd in nohup_cmds:
        print(cmd)
    print("")

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

    # ==== 暴露性状配置 ====
    exposure_disease = "Influenza" # Influenza, Influenza_REP
    exposure_gwas_id = "GCST90432107" # J10_INFLUENZA, GCST90432107
    exposure_trait_id = "GCST90432107" # J10_INFLUENZA, GCST90432107

    # ==== 目录 / 工具路径配置 ====
    data_dir = f'{BASE_DIR}/results/results.temp/mixer'
    n_rep = 20
    threads = 12

    python_bin = "/tools/Python-3.7.9/python"
    mixer_py = "/tools/mixer/precimed/mixer.py"
    lib = "/tools/mixer/src/build/lib/libbgmg.so"
    ldsr_dir = "/references/MiXeR/1000G_EUR_Phase3_plink"

    use_docker = True
    docker_image = "ghcr.io/precimed/gsa-mixer:2.2.1"

    mode = "script"
    out_script = f"{BASE_DIR}/run_fit2_all.sh"

    split_per_dataset = True
    shells_dir = f"/home/xab/influenza_mental_genomics/scripts/MiXeR-fit2-shells/{exposure_trait_id}"
    shells_dir = shells_dir if os.path.exists(os.path.dirname(shells_dir)) else "/home/hello/influenza_mental_genomics/scripts/MiXeR-fit2-shells/{exposure_trait_id}"
    nohup_summary_name = "run_all_fit2_nohup.sh"

    generate_script(
        xlsx_path = xlsx_path,
        sheet = sheet,
        keep = keep,
        exclude_influenza = False,
        exposure_disease = exposure_disease,
        exposure_gwas_id = exposure_gwas_id,
        exposure_trait_id = exposure_trait_id,
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

    if split_per_dataset:
        generate_split_scripts(
            xlsx_path = xlsx_path,
            sheet = sheet,
            keep = keep,
            exclude_influenza = False,
            exposure_disease = exposure_disease,
            exposure_gwas_id = exposure_gwas_id,
            exposure_trait_id = exposure_trait_id,
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