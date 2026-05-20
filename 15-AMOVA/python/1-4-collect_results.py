"""
读取所有场景的 AMOVA 结果 TSV，生成与 Arlequin 预期格式一致的 Markdown 表格。

输出格式示例：
  | Groupings | Number of Populations | Number of Groups | Among Groups |
  | Among Populations within Groups | Within Populations | P-Value | Other |
"""

import argparse
import logging
import math
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

log = logging.getLogger(__name__)


# ──────────────────────────────────────────────────────────────────────
# 辅助函数
# ──────────────────────────────────────────────────────────────────────

def _load_scenarios(scenarios_path: str) -> list[dict]:
    """解析 scenarios.yaml，返回场景列表。"""
    with open(scenarios_path, encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    return cfg.get("scenarios", [])


def _fmt_pct(val: float) -> str:
    """格式化百分比：保留2位小数；NaN 显示为 '-'。"""
    if val is None or (isinstance(val, float) and math.isnan(val)):
        return "-"
    return f"{float(val):.2f}"


def _fmt_pvalue(val: float) -> str:
    """
    格式化 P 值：
    - P < 0.001 显示为 '< 0.001'
    - P < 0.05  保留5位小数
    - 否则保留4位小数
    """
    if val is None or (isinstance(val, float) and math.isnan(val)):
        return "-"
    v = float(val)
    if v < 0.001:
        return "< 0.001"
    if v < 0.05:
        return f"{v:.5f}"
    return f"{v:.4f}"


def _fmt_fstat(val: float) -> str:
    """格式化 F 统计量：保留5位小数；NaN 显示为空字符串。"""
    if val is None or (isinstance(val, float) and math.isnan(val)):
        return ""
    return f"{float(val):.5f}"


def _build_row(label: str, result: dict) -> dict:
    """将单个场景结果转换为表格行字典。"""
    two_level = result.get("two_level", True)
    n_pops = int(result.get("n_pops", 0))
    n_groups = int(result.get("n_groups", 1))

    # 百分比列
    pct_among_groups = result.get("pct_among_groups", float("nan"))
    pct_among_pops = result.get("pct_among_pops", float("nan"))
    pct_within = result.get("pct_within", float("nan"))

    # 2级 AMOVA 无大组间方差
    if two_level or n_groups <= 1:
        col_among_groups = "-"
    else:
        col_among_groups = _fmt_pct(pct_among_groups)

    col_among_pops = _fmt_pct(pct_among_pops)
    col_within = _fmt_pct(pct_within)
    col_pvalue = _fmt_pvalue(result.get("p_value"))

    # "Other" 列：FSC、FST、FCT（3级 AMOVA 时显示）
    fst = result.get("fst", float("nan"))
    fsc = result.get("fsc", float("nan"))
    fct = result.get("fct", float("nan"))

    if two_level or n_groups <= 1:
        col_other = "-"
    else:
        parts = []
        if _fmt_fstat(fsc):
            parts.append(f"FSC: {_fmt_fstat(fsc)}")
        if _fmt_fstat(fst):
            parts.append(f"FST: {_fmt_fstat(fst)}")
        if _fmt_fstat(fct):
            parts.append(f"FCT: {_fmt_fstat(fct)}")
        col_other = "<br>".join(parts) if parts else "-"

    return {
        "Groupings": label,
        "Number of Populations": n_pops,
        "Number of Groups": n_groups,
        "Among Groups": col_among_groups,
        "Among Populations within Groups": col_among_pops,
        "Within Populations": col_within,
        "P-Value": col_pvalue,
        "Other": col_other,
    }


def _rows_to_markdown(rows: list[dict]) -> str:
    """将行列表转换为 Markdown 表格字符串。"""
    if not rows:
        return ""

    headers = list(rows[0].keys())
    lines = []

    # 表头
    lines.append("| " + " | ".join(headers) + " |")
    # 分隔线
    sep = [":---" if h == "Groupings" else "---:" for h in headers]
    lines.append("| " + " | ".join(sep) + " |")
    # 数据行
    for row in rows:
        cells = [str(row.get(h, "")) for h in headers]
        lines.append("| " + " | ".join(cells) + " |")

    return "\n".join(lines)


# ──────────────────────────────────────────────────────────────────────
# 核心业务函数
# ──────────────────────────────────────────────────────────────────────

def run(
    results_dir: str,
    scenarios_path: str,
    output_md: str,
) -> None:
    """
    汇总所有场景结果，生成 Markdown 报告。

    Args:
        results_dir    — 场景结果目录（每个场景一个子目录，含 amova_result.tsv）
        scenarios_path — scenarios.yaml 文件路径
        output_md      — 输出 Markdown 文件路径
    """
    log.info("读取场景配置: %s", scenarios_path)
    scenarios = _load_scenarios(scenarios_path)

    rows = []
    missing = []

    for scenario in scenarios:
        name = scenario["name"]
        label = scenario.get("label", name)
        result_path = Path(results_dir) / name / "amova_result.tsv"

        if not result_path.exists():
            log.warning("场景 '%s' 结果文件不存在: %s", name, result_path)
            missing.append(name)
            continue

        df = pd.read_csv(result_path, sep="\t")
        if df.empty:
            log.warning("场景 '%s' 结果文件为空: %s", name, result_path)
            continue

        result = df.iloc[0].to_dict()
        rows.append(_build_row(label, result))
        log.info("已加载场景: %s（FST=%.5f，P=%s）",
                 label, float(result.get("fst", float("nan"))),
                 _fmt_pvalue(result.get("p_value")))

    if not rows:
        raise ValueError("所有场景均无有效结果，无法生成报告")

    if missing:
        log.warning("以下场景结果缺失，未纳入报告: %s", missing)

    # 生成 Markdown 内容
    md_table = _rows_to_markdown(rows)
    md_content = (
        "Table 1. AMOVA results based on different groups.\n\n"
        + md_table
        + "\n"
    )

    # 写入文件
    Path(output_md).parent.mkdir(parents=True, exist_ok=True)
    with open(output_md, "w", encoding="utf-8") as fh:
        fh.write(md_content)

    log.info("已生成 Markdown 报告: %s（%d 行）", output_md, len(rows))


# ──────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="汇总所有 AMOVA 场景结果，生成 Markdown 报告表格"
    )
    p.add_argument("--results-dir", required=True,
                   help="场景结果目录（各场景子目录含 amova_result.tsv）")
    p.add_argument("--scenarios", required=True,
                   help="scenarios.yaml 文件路径")
    p.add_argument("--output-md", required=True,
                   help="输出 Markdown 报告路径")
    return p


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    args = build_parser().parse_args(argv)
    run(
        results_dir=args.results_dir,
        scenarios_path=args.scenarios,
        output_md=args.output_md,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
