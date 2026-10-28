# -*- coding: utf-8 -*-
"""
Created on Thu Jul 31 11:09:00 2026
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


# ==== 特例映射：与 fit1/fit2/test1/test2 脚本保持一致 ====
DISEASE_DIR_OVERRIDES = {
    "influenza": "_exposure_global",
}

NO_MHC_DISEASES = {
    "influenza",
}


def load_table(xlsx_path, sheet_name, keep, exclude_ids = None):
    df = pd.read_excel(xlsx_path, sheet_name = sheet_name)
    df.columns = [str(c).strip() for c in df.columns]
    for col in ("disease", "gwas_id", "data", "group"):
        if col in df.columns:
            df[col] = df[col].astype(str).str.strip()

    if exclude_ids:
        exclude_ids = [str(x).strip() for x in exclude_ids if str(x).strip()]
        if "gwas_id" in df.columns and exclude_ids:
            df = df[~df["gwas_id"].isin(exclude_ids)]

    df = df[df["group"].isin(["Primary", "Replication", "Sensitivity"])]
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
    key = disease.strip().lower()
    return DISEASE_DIR_OVERRIDES.get(key, safe_name(disease))


def get_trait_file_suffix(disease):
    key = disease.strip().lower()
    return ".noMHC.csv.gz" if key in NO_MHC_DISEASES else ".noMHC.csv.gz"


def build_header_lines(base_dir, use_docker, docker_image,
                        python_bin, mfigure_py):
    """
    生成脚本 header：BASE_DIR
    以及 PYTHON / MFIGURE 变量，兼容 docker 与直连 tools 两种方式。
    """
    lines = []
    lines.append("#!/bin/bash")
    lines.append("set -euo pipefail")
    lines.append("")
    lines.append(f'BASE_DIR="{base_dir}"')
    lines.append("")

    if use_docker:
        lines.append(f'export DOCKER_RUN="docker run -v {base_dir}:{base_dir} -v /references:/references -w {base_dir}"')
        lines.append(f'export MFIGURE="$DOCKER_RUN {docker_image} python {mfigure_py}"')
        lines.append('PYTHON=""')
        lines.append("")
    else:
        lines.append(f'PYTHON="{python_bin}"')
        lines.append(f'MFIGURE="{mfigure_py}"')
        lines.append("")

    return lines


def check_pair_completed(base_dir, pair_name, exposure_trait_id, outcome_trait_id, n_rep, stage):
    """
    检查某个 pair 的 fit2/test2 是否全部跑完。
    stage: 'fit2' 或 'test2'
    遍历 rep 1..n_rep，检查 .json 参数文件是否全部存在。
    """
    stage_dir = os.path.join(base_dir, "results", "results", "mixer", stage, pair_name)
    prefix = f"{exposure_trait_id}_vs_{outcome_trait_id}"

    for rep in range(1, n_rep + 1):
        json_file = os.path.join(stage_dir, f"{prefix}.rep{rep}.json")
        if not os.path.isfile(json_file):
            return False
    return True


def build_figures_block(exposure_disease, exposure_gwas_id, exposure_trait_id,
                        outcome_disease, outcome_gwas_id, outcome_trait_id,
                        exts, convert_pdf):
    """
    生成结果整合与可视化的命令块。
    1. combine fit2 replicates
    2. combine test2 replicates
    3. mixer_figures two 生成双性状图（支持多格式 ext）
    4. 可选：svg -> pdf 转换
    """
    safe_exposure = safe_name(exposure_disease)
    safe_outcome = safe_name(outcome_disease)
    pair_name = f"{safe_exposure}({exposure_trait_id})_vs_{safe_outcome}({outcome_trait_id})"

    fit2_dir = f"${{BASE_DIR}}/results/results/mixer/fit2/{pair_name}"
    test2_dir = f"${{BASE_DIR}}/results/results/mixer/test2/{pair_name}"
    redir = f"${{BASE_DIR}}/results/results/mixer/figures/{pair_name}"

    ext_str = " ".join(exts)

    lines = []
    lines.append(f"# ==== {exposure_disease} vs {outcome_disease} ({exposure_gwas_id} vs {outcome_gwas_id}) : figures ====")
    lines.append(f'PAIR="{pair_name}"')
    lines.append(f'FIT2_DIR="{fit2_dir}"')
    lines.append(f'TEST2_DIR="{test2_dir}"')
    lines.append(f'REDIR="{redir}"')
    lines.append('mkdir -p "${REDIR}"')
    lines.append("")

    lines.append(f'echo "[${{PAIR}}] Combining fit2 replicates..."')
    lines.append(
        '  ${PYTHON} ${MFIGURE} combine \\\n'
        f'        --json "${{FIT2_DIR}}/{exposure_trait_id}_vs_{outcome_trait_id}.rep@.json" \\\n'
        f'        --out "${{REDIR}}/${{PAIR}}.fit"'
    )

    lines.append(f'echo "[${{PAIR}}] Combining test2 replicates..."')
    lines.append(
        '  ${PYTHON} ${MFIGURE} combine \\\n'
        f'        --json "${{TEST2_DIR}}/{exposure_trait_id}_vs_{outcome_trait_id}.rep@.json" \\\n'
        f'        --out "${{REDIR}}/${{PAIR}}.test"'
    )

    lines.append(f'echo "[${{PAIR}}] Generating two-trait figures ({ext_str})..."')
    lines.append(
        '  ${PYTHON} ${MFIGURE} two \\\n'
        f'        --json-fit "${{REDIR}}/${{PAIR}}.fit.json" \\\n'
        f'        --json-test "${{REDIR}}/${{PAIR}}.test.json" \\\n'
        f'        --out "${{REDIR}}/${{PAIR}}" \\\n'
        f'        --statistic mean std \\\n'
        f'        --ext {ext_str}'
    )

    # ==== 额外：svg -> pdf 转换（因为 mixer_figures 本身不支持 pdf）====
    if convert_pdf and "svg" in exts:
        lines.append(f'echo "[${{PAIR}}] Converting SVG to PDF..."')
        lines.append(
            'for svg_file in "${REDIR}"/*.svg; do\n'
            '  [ -e "$svg_file" ] || continue\n'
            '  pdf_file="${svg_file%.svg}.pdf"\n'
            '  if command -v inkscape &> /dev/null; then\n'
            '    inkscape "$svg_file" --export-filename="$pdf_file" > /dev/null 2>&1\n'
            '  elif command -v rsvg-convert &> /dev/null; then\n'
            '    rsvg-convert -f pdf -o "$pdf_file" "$svg_file"\n'
            '  else\n'
            '    echo "[WARN] Neither inkscape nor rsvg-convert found, skip PDF conversion for $svg_file"\n'
            '  fi\n'
            'done'
        )

    lines.append("")
    return lines


def iter_selected_datasets(df, diseases, keep):
    """
    统一遍历所有 (disease, row) 组合。
    返回 generator，元素为 (disease, gwas_id, trait_id) 三元组。
    """
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
                     base_dir, n_rep,
                     use_docker, docker_image,
                     python_bin, mfigure_py,
                     exts, convert_pdf,
                     mode, out_script):
    """
    生成单个合并的 shell 脚本（所有已完成 pair 的 figures 命令都写在一个文件里）。
    运行前会自动检测 fit2/test2 是否全部完成，未完成的 pair 会被跳过并打印提示。
    """

    df = load_table(xlsx_path, sheet, keep, exclude_ids=[exposure_trait_id, exposure_gwas_id])

    # 结局性状列表（排除暴露性状本身）
    diseases = list(dict.fromkeys(df["disease"].dropna().tolist()))
    diseases = [d for d in diseases if d.strip().lower() != exposure_disease.strip().lower()]

    lines = build_header_lines(base_dir, use_docker, docker_image,
                                python_bin, mfigure_py)

    n_jobs = 0
    skipped = []

    for outcome_disease, outcome_gwas_id, outcome_trait_id in iter_selected_datasets(df, diseases, keep):
        safe_exposure = safe_name(exposure_disease)
        safe_outcome = safe_name(outcome_disease)
        pair_name = f"{safe_exposure}({exposure_trait_id})_vs_{safe_outcome}({outcome_trait_id})"

        # ==== 核心检测：只跑已经全部跑完的 pair ====
        fit2_ok = check_pair_completed(base_dir, pair_name,
                                       exposure_trait_id, outcome_trait_id, n_rep, 'fit2')
        test2_ok = check_pair_completed(base_dir, pair_name,
                                        exposure_trait_id, outcome_trait_id, n_rep, 'test2')

        if not fit2_ok:
            skipped.append(f"{pair_name} (fit2 incomplete)")
            print(f"[SKIP] {pair_name} fit2 结果未全部就绪 ({n_rep} reps)，跳过", file=sys.stderr)
            continue
        if not test2_ok:
            skipped.append(f"{pair_name} (test2 incomplete)")
            print(f"[SKIP] {pair_name} test2 结果未全部就绪 ({n_rep} reps)，跳过", file=sys.stderr)
            continue

        block = build_figures_block(
            exposure_disease = exposure_disease,
            exposure_gwas_id = exposure_gwas_id,
            exposure_trait_id = exposure_trait_id,
            outcome_disease = outcome_disease,
            outcome_gwas_id = outcome_gwas_id,
            outcome_trait_id = outcome_trait_id,
            exts = exts,
            convert_pdf = convert_pdf,
        )
        lines.extend(block)
        n_jobs += 1

    # 汇总：将所有 svg/png/pdf 复制到统一目录
    summary_dir = f"${{BASE_DIR}}/results/results/mixer/figures_summary"
    lines.append("# ==== 汇总所有图片到统一目录 ====")
    lines.append(f'SUMMARY_DIR="{summary_dir}"')
    lines.append('mkdir -p "${SUMMARY_DIR}"')
    lines.append('find "${BASE_DIR}/results/results/mixer/figures" -type f \\( -name "*.png" -o -name "*.svg" -o -name "*.pdf" \\) -exec cp -f {} "${SUMMARY_DIR}" \\;')
    lines.append('echo "[OK] All figures copied to ${SUMMARY_DIR}"')
    lines.append("")

    script_content = "\n".join(lines)

    if mode == "script":
        with open(out_script, "w", encoding = "utf-8") as f:
            f.write(script_content)
        os.chmod(out_script, 0o755)
        print(f"[OK] 已生成合并脚本: {out_script}")
        print(f"     暴露性状: {exposure_disease} ({exposure_gwas_id})")
        print(f"     输出格式: {exts}" + (" + pdf(convert)" if convert_pdf else ""))
        print(f"     共 {n_jobs} 个pair已就绪并加入脚本，跳过 {len(skipped)} 个未完成的pair")
        if skipped:
            print("     [INFO] 以下pair尚未完成，已跳过：")
            for s in skipped:
                print(f"            - {s}")
        print("     请确认上述结果后执行：")
        print(f"     bash {out_script}")
    else:
        print(f"[RUN] 直接执行模式，共 {n_jobs} 个pair")
        proc = subprocess.run(["bash", "-c", script_content])
        sys.exit(proc.returncode)


def generate_split_scripts(xlsx_path, sheet, keep, exclude_influenza,
                            exposure_disease, exposure_gwas_id, exposure_trait_id,
                            base_dir, n_rep,
                            use_docker, docker_image,
                            python_bin, mfigure_py,
                            exts, convert_pdf,
                            shells_dir, nohup_summary_name):
    """
    一个已完成 pair 一个 shell 脚本：
      <shells_dir>/<idx>.<exposure>_vs_<outcome>_figures.sh
    并打印 / 汇总每个脚本对应的 nohup 后台运行命令。
    未完成的 pair 会被跳过。
    """

    df = load_table(xlsx_path, sheet, keep, exclude_ids=[exposure_trait_id, exposure_gwas_id])

    diseases = list(dict.fromkeys(df["disease"].dropna().tolist()))
    diseases = [d for d in diseases if d.strip().lower() != exposure_disease.strip().lower()]

    os.makedirs(shells_dir, exist_ok = True)

    nohup_cmds = []
    idx = 0
    skipped = []

    for outcome_disease, outcome_gwas_id, outcome_trait_id in iter_selected_datasets(df, diseases, keep):
        safe_exposure = safe_name(exposure_disease)
        safe_outcome = safe_name(outcome_disease)
        pair_name = f"{safe_exposure}({exposure_trait_id})_vs_{safe_outcome}({outcome_trait_id})"

        # 检测是否全部完成
        fit2_ok = check_pair_completed(base_dir, pair_name,
                                       exposure_trait_id, outcome_trait_id, n_rep, 'fit2')
        test2_ok = check_pair_completed(base_dir, pair_name,
                                        exposure_trait_id, outcome_trait_id, n_rep, 'test2')

        if not fit2_ok:
            skipped.append(f"{pair_name} (fit2 incomplete)")
            print(f"[SKIP] {pair_name} fit2 结果未全部就绪，跳过", file=sys.stderr)
            continue
        if not test2_ok:
            skipped.append(f"{pair_name} (test2 incomplete)")
            print(f"[SKIP] {pair_name} test2 结果未全部就绪，跳过", file=sys.stderr)
            continue

        idx += 1

        header = build_header_lines(base_dir, use_docker, docker_image,
                                     python_bin, mfigure_py)
        block = build_figures_block(
            exposure_disease = exposure_disease,
            exposure_gwas_id = exposure_gwas_id,
            exposure_trait_id = exposure_trait_id,
            outcome_disease = outcome_disease,
            outcome_gwas_id = outcome_gwas_id,
            outcome_trait_id = outcome_trait_id,
            exts = exts,
            convert_pdf = convert_pdf,
        )

        script_content = "\n".join(header + block)

        script_name = f"{idx}.{safe_exposure}_vs_{safe_outcome}_{outcome_gwas_id}_figures.sh"
        script_path = os.path.join(shells_dir, script_name)

        with open(script_path, "w", encoding = "utf-8") as f:
            f.write(script_content)
        os.chmod(script_path, 0o755)

        log_name = script_name.replace(".sh", ".log")
        log_path = os.path.join(shells_dir, log_name)

        nohup_cmd = f"nohup bash {script_path} > {log_path} 2>&1 &"
        nohup_cmds.append(nohup_cmd)

    # 打印全部 nohup 命令
    print(f"[OK] 已在 {shells_dir} 下生成 {idx} 个独立 figures 脚本")
    print(f"     暴露性状: {exposure_disease} ({exposure_gwas_id})")
    print(f"     输出格式: {exts}" + (" + pdf(convert)" if convert_pdf else ""))
    print(f"     跳过 {len(skipped)} 个未完成的pair")
    if skipped:
        print("     [INFO] 以下pair尚未完成，已跳过：")
        for s in skipped:
            print(f"            - {s}")
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

    # ==== 暴露性状配置 ====
    # 暴露性状固定，与所有结局性状做双变量分析后的结果整合
    exposure_disease = "Influenza"         # Influenza
    exposure_gwas_id = "GCST90432107"     # J10_INFLUENZA, GCST90432107
    exposure_trait_id = "GCST90432107"    # J10_INFLUENZA, GCST90432107
    
    # ==== 目录 / 工具路径配置 ====
    # n_rep 必须与 fit2/test2 脚本中的重复次数一致，用于检测是否全部完成
    n_rep = 20

    # 直连 tools 路径方式（docker=False 时生效）
    python_bin = "/tools/Python-3.7.9/python"
    mfigure_py = "/tools/mixer/precimed/mixer_figures.py"

    # docker 方式
    use_docker = True
    docker_image = "ghcr.io/precimed/gsa-mixer:2.2.1"

    # ==== figures 特有配置 ====
    # 输出图片格式：mixer_figures 支持 png / svg，可传多个
    # 例如 ["svg"] 只出 svg；["png", "svg"] 同时出两种
    exts = ["png", "svg"]

    # 是否额外将 svg 转换为 pdf（需要容器内有 inkscape 或 rsvg-convert）
    # 如果 docker 镜像里没有这两个工具，可以设为 False，手动在外面转
    convert_pdf = False

    # 输出模式: "script" 只生成合并 shell 脚本；"run" 生成后直接执行合并脚本
    mode = "script"
    out_script = f"{BASE_DIR}/run_figures_all.sh"

    # ==== 是否为每个 pair 单独生成一个 shell 脚本 ====
    split_per_dataset = True
    shells_dir = f"/home/xab/influenza_mental_genomics/scripts/MiXeR-figures-shells/{exposure_trait_id}"
    shells_dir = shells_dir if os.path.exists(os.path.dirname(os.path.dirname(shells_dir))) \
        else f"/home/hello/influenza_mental_genomics/scripts/MiXeR-figures-shells/{exposure_trait_id}"
    nohup_summary_name = "run_all_figures_nohup.sh"

    # ---- 合并脚本（原有功能，保留） ----
    generate_script(
        xlsx_path = xlsx_path,
        sheet = sheet,
        keep = keep,
        exclude_influenza = False,
        exposure_disease = exposure_disease,
        exposure_gwas_id = exposure_gwas_id,
        exposure_trait_id = exposure_trait_id,
        base_dir = BASE_DIR,
        n_rep = n_rep,
        use_docker = use_docker,
        docker_image = docker_image,
        python_bin = python_bin,
        mfigure_py = mfigure_py,
        exts = exts,
        convert_pdf = convert_pdf,
        mode = mode,
        out_script = out_script,
    )

    # ---- 拆分脚本：一个已完成 pair 一个 shell，并打印 nohup 命令 ----
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
            n_rep = n_rep,
            use_docker = use_docker,
            docker_image = docker_image,
            python_bin = python_bin,
            mfigure_py = mfigure_py,
            exts = exts,
            convert_pdf = convert_pdf,
            shells_dir = shells_dir,
            nohup_summary_name = nohup_summary_name,
        )