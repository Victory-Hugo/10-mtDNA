from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


class _ColorFormatter(logging.Formatter):
    _colors = {
        logging.DEBUG: "\033[36m",
        logging.INFO: "\033[32m",
        logging.WARNING: "\033[33m",
        logging.ERROR: "\033[31m",
        logging.CRITICAL: "\033[41m\033[97m",
    }
    _reset = "\033[0m"

    def format(self, record: logging.LogRecord) -> str:
        message = super().format(record)
        color = self._colors.get(record.levelno, "")
        if not color:
            return message
        return f"{color}{message}{self._reset}"


def configure_logging(log_path: Path, level: str = "INFO") -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_level = logging.getLevelName(level.upper())

    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setFormatter(
        logging.Formatter("%(asctime)s | %(levelname)s | %(name)s | %(message)s")
    )

    console_handler = logging.StreamHandler()
    console_handler.setFormatter(
        _ColorFormatter("%(asctime)s | %(levelname)s | %(message)s", datefmt="%H:%M:%S")
    )

    logging.basicConfig(level=log_level, handlers=[file_handler, console_handler], force=True)


def _read_table(table_path: Path) -> pd.DataFrame:
    logger.info("读取样本表: %s", table_path)
    if table_path.suffix.lower() == ".tsv":
        df = pd.read_csv(table_path, sep="\t")
        logger.info("样本表读取完成: %s 行, %s 列", df.shape[0], df.shape[1])
        return df
    if table_path.suffix.lower() == ".csv":
        df = pd.read_csv(table_path)
        logger.info("样本表读取完成: %s 行, %s 列", df.shape[0], df.shape[1])
        return df
    raise ValueError("样本信息表格仅支持 .csv 或 .tsv 格式")


def run(
    input_path: str,
    output_path: str,
    group_cols: list[str],
) -> None:
    input_path = Path(input_path)
    output_path = Path(output_path)

    logger.info("开始过滤样本表")
    logger.info("输入: %s", input_path)
    logger.info("输出: %s", output_path)
    logger.info("分组列: %s", ", ".join(group_cols))

    sample_df = _read_table(input_path)

    missing_cols = [col for col in group_cols if col not in sample_df.columns]
    if missing_cols:
        raise KeyError(
            "样本表缺少列: "
            f"{', '.join(missing_cols)}。可用列: {list(sample_df.columns)}"
        )

    logger.info("分组列数量: %s", len(group_cols))

    filtered = sample_df.copy()
    for col in group_cols:
        counts = filtered[col].value_counts(dropna=True)
        logger.info("%s 可用群体数: %s", col, len(counts))
        filtered = filtered[filtered[col].map(counts).fillna(0).ge(2)].copy()

    logger.info("过滤前样本数: %s", len(sample_df))
    logger.info("过滤后样本数: %s", len(filtered))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    filtered.to_csv(output_path, sep="\t", index=False)
    logger.info("过滤结果已写入: %s", output_path)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="过滤样本表：确保各分组列群体数 >= 2")
    parser.add_argument("--input", required=True, help="输入样本表 (.csv/.tsv)")
    parser.add_argument("--output", required=True, help="输出过滤后样本表 (.tsv)")
    parser.add_argument(
        "--group-cols",
        required=True,
        help="分组列名称，逗号分隔",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    configure_logging(Path(args.output).with_suffix(".log"))
    group_cols = [col.strip() for col in args.group_cols.split(",") if col.strip()]
    run(
        input_path=args.input,
        output_path=args.output,
        group_cols=group_cols,
    )


if __name__ == "__main__":
    main()
