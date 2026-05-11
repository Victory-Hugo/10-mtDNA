#!/usr/bin/env python3
"""根据 fastHaN JSON 和 metadata 指定列生成 tcsBU 可视化配置。"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import re
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

log = logging.getLogger(__name__)

DEFAULT_COLOR_PALETTE = "#274753,#297270,#299d8f,#8ab07c,#e7c66b,#f3a361,#e66d50"


class VisualConfigError(ValueError):
    """可视化配置生成失败。"""


def safe_name(name: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", name.strip())
    return safe.strip("_") or "column"


def split_csv(value: str) -> List[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def parse_hex_color(color: str) -> Tuple[int, int, int]:
    color = color.strip()
    if not re.fullmatch(r"#[0-9A-Fa-f]{6}", color):
        raise VisualConfigError(f"颜色必须是 #RRGGBB 十六进制格式: {color}")
    return int(color[1:3], 16), int(color[3:5], 16), int(color[5:7], 16)


def rgb_to_hex(rgb: Tuple[int, int, int]) -> str:
    return "#{:02x}{:02x}{:02x}".format(*rgb)


def interpolate_palette(palette: str, n_colors: int) -> List[str]:
    anchors = [parse_hex_color(color) for color in split_csv(palette)]
    if not anchors:
        raise VisualConfigError("fasthan_color_palette 至少需要 1 个颜色")
    if n_colors <= 0:
        return []
    if n_colors == 1:
        return [rgb_to_hex(anchors[0])]
    if len(anchors) == 1:
        return [rgb_to_hex(anchors[0]) for _ in range(n_colors)]

    result: List[str] = []
    max_anchor = len(anchors) - 1
    for i in range(n_colors):
        position = i * max_anchor / (n_colors - 1)
        left = int(position)
        right = min(left + 1, max_anchor)
        fraction = position - left
        rgb = tuple(
            round(anchors[left][channel] + (anchors[right][channel] - anchors[left][channel]) * fraction)
            for channel in range(3)
        )
        result.append(rgb_to_hex(rgb))
    return result


def load_metadata(path: Path, sample_id_column: str) -> Tuple[List[str], Dict[str, Dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise VisualConfigError(f"metadata 表头为空: {path}")
        if sample_id_column not in reader.fieldnames:
            raise VisualConfigError(f"metadata 缺少样本列: {sample_id_column}")
        rows = {}
        for row in reader:
            sample_id = (row.get(sample_id_column) or "").strip()
            if sample_id:
                rows[sample_id] = row
    return list(reader.fieldnames), rows


def load_hap_samples(json_path: Path) -> List[str]:
    data = json.loads(json_path.read_text(encoding="utf-8"))
    if "nodes" not in data:
        raise VisualConfigError(f"fastHaN JSON 缺少 nodes: {json_path}")
    samples: List[str] = []
    for node in data["nodes"]:
        title = str(node.get("title2", ""))
        for sample_id in title.split(";"):
            sample_id = sample_id.strip()
            if sample_id and not sample_id.startswith("IN"):
                samples.append(sample_id)
    return samples


def write_configs(
    samples: Sequence[str],
    metadata: Dict[str, Dict[str, str]],
    column: str,
    output_prefix: Path,
    unknown_label: str,
    color_palette: str,
) -> Tuple[Path, Path, Path]:
    hapconf = output_prefix.with_name(output_prefix.name + "_hapconf.csv")
    groupconf = output_prefix.with_name(output_prefix.name + "_groupconf.csv")
    updated_groupconf = output_prefix.with_name(output_prefix.name + "_groupconf_updated.csv")

    groups: List[str] = []
    sample_groups: List[Tuple[str, str]] = []
    for sample_id in samples:
        value = (metadata.get(sample_id, {}).get(column) or "").strip() or unknown_label
        value = value.replace(" ", "_")
        sample_groups.append((sample_id, value))
        if value not in groups and value != "Default":
            groups.append(value)

    with hapconf.open("w", encoding="utf-8", newline="\n") as handle:
        for sample_id, group in sample_groups:
            handle.write(f"{sample_id};{group}\n")

    interpolated_colors = interpolate_palette(color_palette, len(groups))

    with groupconf.open("w", encoding="utf-8", newline="\n") as raw_handle, updated_groupconf.open(
        "w", encoding="utf-8", newline="\n"
    ) as updated_handle:
        color_counts: Dict[str, int] = {}
        for group, color in zip(groups, interpolated_colors):
            raw_handle.write(f"{group};#999999;none\n")
            value = "none"
            if color in color_counts:
                color_counts[color] += 1
                value = f"lines-{color_counts[color] * 2 - 1}"
            else:
                color_counts[color] = 1
            updated_handle.write(f"{group};{color};{value}\n")

    return hapconf, groupconf, updated_groupconf


def run(
    json_file: str,
    metadata: str,
    sample_id_column: str,
    columns: str,
    output_dir: str,
    prefix: str,
    unknown_label: str = "Unknown",
    color_palette: str = DEFAULT_COLOR_PALETTE,
) -> int:
    json_path = Path(json_file)
    metadata_path = Path(metadata)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    if not json_path.is_file():
        raise VisualConfigError(f"fastHaN JSON 不存在: {json_path}")
    if not metadata_path.is_file():
        raise VisualConfigError(f"metadata 不存在: {metadata_path}")

    fieldnames, meta_by_id = load_metadata(metadata_path, sample_id_column)
    target_columns = split_csv(columns)
    missing_columns = [column for column in target_columns if column not in fieldnames]
    if missing_columns:
        raise VisualConfigError(f"metadata 缺少可视化列: {', '.join(missing_columns)}")

    samples = load_hap_samples(json_path)
    for column in target_columns:
        output_prefix = output_path / f"{safe_name(prefix)}_{safe_name(column)}"
        write_configs(samples, meta_by_id, column, output_prefix, unknown_label, color_palette)
        log.info("生成 fastHaN 可视化配置: %s / %s", prefix, column)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-file", required=True, help="fastHaN 输出 JSON")
    parser.add_argument("--metadata", required=True, help="metadata TSV")
    parser.add_argument("--sample-id-column", required=True, help="样本 ID 列名")
    parser.add_argument("--columns", required=True, help="逗号分隔的 metadata 列名")
    parser.add_argument("--output-dir", required=True, help="输出目录")
    parser.add_argument("--prefix", required=True, help="输出文件名前缀")
    parser.add_argument("--unknown-label", default="Unknown", help="缺失分组标签")
    parser.add_argument(
        "--color-palette",
        default=DEFAULT_COLOR_PALETTE,
        help="逗号分隔的 #RRGGBB 颜色锚点，会按分组数量线性插值",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(**vars(args))
    except VisualConfigError as exc:
        parser.error(str(exc))
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
