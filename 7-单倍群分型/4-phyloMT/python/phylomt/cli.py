"""命令行入口。"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

from .tree import DEFAULT_TREES_DIR, list_installed_trees, load_tree_bundle


DEFAULT_HITS = 50
DEFAULT_HET_LEVEL = 0.9


class ConsoleUI:
    """运行过程终端输出。"""

    def __init__(self) -> None:
        self.use_color = self._detect_color()
        self.use_unicode = self._detect_unicode()
        self.started_at = time.monotonic()
        self.reset = "\033[0m" if self.use_color else ""
        self.bold = "\033[1m" if self.use_color else ""
        self.dim = "\033[2m" if self.use_color else ""
        self.blue = "\033[38;5;39m" if self.use_color else ""
        self.cyan = "\033[38;5;45m" if self.use_color else ""
        self.green = "\033[38;5;42m" if self.use_color else ""
        self.yellow = "\033[38;5;220m" if self.use_color else ""
        self.red = "\033[38;5;196m" if self.use_color else ""
        self.magenta = "\033[38;5;141m" if self.use_color else ""
        self.rule_char = "-" if not self.use_unicode else "─"
        self.edge = "|" if not self.use_unicode else "│"
        self.join = "-" if not self.use_unicode else "┄"
        self.stage_icon = ">" if not self.use_unicode else "▶"
        self.ok_icon = "+" if not self.use_unicode else "✔"

    @staticmethod
    def _detect_color() -> bool:
        if os.environ.get("NO_COLOR") or os.environ.get("TERM") == "dumb":
            return False
        return sys.stderr.isatty()

    @staticmethod
    def _detect_unicode() -> bool:
        if os.environ.get("PHYLO_UI_ASCII_ONLY"):
            return False
        locale_hint = " ".join(
            value for value in (os.environ.get("LC_ALL"), os.environ.get("LC_CTYPE"), os.environ.get("LANG")) if value
        )
        return "UTF-8" in locale_hint.upper() or "UTF8" in locale_hint.upper()

    def paint(self, text: str, color: str = "", *, bold: bool = False, dim: bool = False) -> str:
        prefix = ""
        if bold:
            prefix += self.bold
        if dim:
            prefix += self.dim
        prefix += color
        return f"{prefix}{text}{self.reset}" if prefix else text

    def write(self, message: str = "") -> None:
        print(message, file=sys.stderr)

    def rule(self) -> None:
        self.write(self.paint(self.rule_char * 72, self.dim))

    def logo(self) -> None:
        if self.use_unicode:
            lines = [
                (self.cyan, "                 ╭──●"),
                (self.cyan, "           ╭─────┤"),
                (self.blue, "      ╭────┤     ╰──●"),
                (self.blue, "  ╭───┤"),
                (self.magenta, "──┤   ╰────────╮"),
                (self.magenta, "  │            ├──●"),
                (self.blue, "  │   phyloMT  │"),
                (self.blue, "──┤            ├──●"),
                (self.cyan, "  │   mtDNA    ╯"),
                (self.cyan, "  ╰────────────────●"),
            ]
        else:
            lines = [
                ("", "                  ,-o"),
                ("", "            ,-----|"),
                ("", "       ,----|     '-o"),
                ("", "   ,---|"),
                ("", "---|   '--------."),
                ("", "   |            |-o"),
                ("", "   |   phyloMT  |"),
                ("", "---|            |-o"),
                ("", "   |   mtDNA    '"),
                ("", "   '---------------o"),
            ]
        for color, line in lines:
            self.write(self.paint(line, color, bold=True))

    def section(self, title: str, subtitle: str = "") -> None:
        self.rule()
        self.write(self.paint(f"{self.edge} {title}", self.blue, bold=True))
        if subtitle:
            self.write(f"{self.paint(self.edge, self.dim)} {subtitle}")
        self.rule()

    def stage(self, step: str, title: str) -> None:
        self.section(f"{self.stage_icon} {step}", title)

    def kv(self, key: str, value: str) -> None:
        key_text = f"{key + ':':18s}"
        self.write(f"{self.paint(self.join, self.dim)} {self.paint(key_text, bold=True)} {value}")

    def ok(self, message: str) -> None:
        elapsed = time.monotonic() - self.started_at
        self.write(self.paint(f"{self.ok_icon} {message} ({elapsed:.1f}s)", self.green))

    def warn(self, message: str) -> None:
        self.write(self.paint(f"! {message}", self.yellow))

    def error(self, message: str) -> None:
        self.write(self.paint(f"x {message}", self.red, bold=True))


def main() -> int:
    """命令行主入口。"""

    parser = build_parser()
    if any(arg in {"-h", "--help"} for arg in sys.argv[1:]):
        print_pretty_help()
        return 0

    args = parser.parse_args()
    try:
        return run_classify(args)
    except ValueError as exc:
        print_pretty_error(str(exc))
        return 2


def build_parser() -> argparse.ArgumentParser:
    """构建命令行解析器。"""

    parser = argparse.ArgumentParser(prog="phyloMT", description="基于 Python 的 mtDNA 单倍群分型工具", add_help=False)
    parser.add_argument("-h", "--help", action="store_true", help="显示帮助信息并退出")
    parser.add_argument("--tree", help="树 ID，例如 phylotree-rcrs@17.2；可先用 --list-trees 查看")
    parser.add_argument("--input", help="输入文件路径")
    parser.add_argument("--output", help="输出文件路径")
    parser.add_argument("--hits", type=int, default=None, help=f"输出前 N 个候选，默认 {DEFAULT_HITS}")
    parser.add_argument("--het-level", type=float, default=None, help=f"VCF 杂合阈值，默认 {DEFAULT_HET_LEVEL}")
    parser.add_argument("--extended-report", action="store_true", help="输出扩展报告")
    parser.add_argument("--write-fasta", action="store_true", help="重建 FASTA 输出")
    parser.add_argument("--threads", type=int, default=1, help="halign4 比对线程数，默认 1")
    parser.add_argument("--chip", action="store_true", help="芯片/微阵列模式：仅在 VCF 覆盖位点评估期望变异")
    parser.add_argument("--trees-dir", default=None, help="树目录 默认情况下为项目 data/trees 目录")
    parser.add_argument("--list-trees", action="store_true", help="列出 trees-dir 下所有可用树并退出")
    return parser


def print_usage_block(ui: ConsoleUI) -> None:
    """输出命令用法。"""

    ui.section("用法", "非交互式 mtDNA 单倍群分型命令")
    ui.write(ui.paint("python3 phyloMT.py \\", ui.cyan, bold=True))
    ui.write("  --tree <树ID> \\")
    ui.write("  --input <输入文件> \\")
    ui.write("  --output <输出TSV> [选项]")


def print_argument_block(ui: ConsoleUI) -> None:
    """输出参数说明。"""

    ui.section("必需参数")
    ui.kv("--tree", "树 ID，例如 phylotree-rsrs@17.1；可用 --list-trees 查看")
    ui.kv("--input", "输入文件路径，支持 FASTA / HSD / VCF")
    ui.kv("--output", "输出 TSV 文件路径")

    ui.section("常用选项")
    ui.kv("--hits", f"输出前 N 个候选，默认 {DEFAULT_HITS}")
    ui.kv("--threads", "halign4 比对线程数，默认 1")
    ui.kv("--het-level", f"VCF 杂合阈值，默认 {DEFAULT_HET_LEVEL}")
    ui.kv("--extended-report", "同时输出 .extended.tsv 扩展报告")
    ui.kv("--write-fasta", "根据样本变异重建 FASTA")
    ui.kv("--chip", "芯片/微阵列 VCF 模式")
    ui.kv("--trees-dir", "手动指定树目录，默认 data/trees")
    ui.kv("--list-trees", "列出可用树并退出")
    ui.kv("-h, --help", "显示当前帮助界面")


def print_example_block(ui: ConsoleUI) -> None:
    """输出示例命令。"""

    ui.section("示例")
    ui.write("python3 phyloMT.py \\")
    ui.write("  --tree phylotree-rsrs@17.1 \\")
    ui.write("  --input input/hsd/example.hsd \\")
    ui.write("  --output output/example_hsd_output.tsv \\")
    ui.write("  --hits 1 --threads 8 --extended-report")
    ui.write("")


def print_pretty_help() -> None:
    """输出美化帮助。"""

    ui = ConsoleUI()
    ui.logo()
    ui.section("phyloMT mtDNA Haplogroup Classifier", "基于 Python 的 mtDNA 单倍群分型工具")
    print_usage_block(ui)
    print_argument_block(ui)
    print_example_block(ui)


def print_pretty_error(message: str) -> None:
    """输出美化错误信息。"""

    ui = ConsoleUI()
    ui.logo()
    ui.section("phyloMT mtDNA Haplogroup Classifier", "命令参数不完整或不合法")
    ui.error(message)
    print_usage_block(ui)
    ui.warn("使用 python3 phyloMT.py --help 查看完整帮助；使用 --list-trees 查看可用树。")


def resolve_trees_dir(trees_dir: str | None) -> Path:
    """确定树目录。"""

    if trees_dir:
        return Path(trees_dir).resolve()
    return DEFAULT_TREES_DIR.resolve()


def split_tree_id(full_id: str) -> tuple[str, str]:
    """拆分树 ID 和版本。"""

    if "@" not in full_id:
        raise ValueError("树 ID 必须写成 <id@version>。")
    tree_id, version = full_id.split("@", 1)
    return tree_id, version


def validate_classify_args(args: argparse.Namespace) -> None:
    """校验分类命令所需参数。"""

    missing: list[str] = []
    if not args.tree:
        missing.append("--tree")
    if not args.input:
        missing.append("--input")
    if not args.output:
        missing.append("--output")
    if missing:
        missing_text = ", ".join(missing)
        raise ValueError(f"执行分类时必须提供 {missing_text}。")


def print_available_trees(trees_dir: Path) -> None:
    """输出可用树列表。"""

    entries = list_installed_trees(trees_dir)
    if not entries:
        print(f"未在 {trees_dir} 中找到可用树。")
        return

    print(f"可用树目录: {trees_dir}")
    for entry in entries:
        print(entry["full_id"])


def extended_report_path(output_path: Path) -> Path:
    """匹配 report.write_extended_report 的输出命名。"""

    return output_path.with_name(f"{output_path.stem}.extended{output_path.suffix or '.txt'}")


def resolve_installed_tree_dir(tree_id: str, version: str, trees_dir: Path) -> Path:
    """确定已安装树目录。"""

    tree_dir = trees_dir / tree_id / version
    if (tree_dir / "tree.yaml").exists():
        return tree_dir

    installed = [entry["full_id"] for entry in list_installed_trees(trees_dir)]
    if installed:
        installed_text = ", ".join(installed)
        raise ValueError(f"未在 {trees_dir} 中找到 {tree_id}@{version}。当前可用树: {installed_text}")
    raise ValueError(f"未在 {trees_dir} 中找到可用树。")


def run_classify(args: argparse.Namespace) -> int:
    """执行分类。"""

    from .classifier import classify_samples
    from .io.readers import load_samples
    from .report import write_classification_report, write_extended_report, write_reconstructed_fasta

    ui = ConsoleUI()
    trees_dir = resolve_trees_dir(args.trees_dir)
    if args.list_trees:
        print_available_trees(trees_dir)
        return 0

    validate_classify_args(args)
    output_path = Path(args.output).resolve()
    hits = args.hits or DEFAULT_HITS
    het_level = args.het_level if args.het_level is not None else DEFAULT_HET_LEVEL

    ui.logo()
    ui.section("phyloMT mtDNA Haplogroup Classifier", "基于本地系统发育树进行单倍群分型")
    ui.kv("Tree", args.tree)
    ui.kv("Input", str(Path(args.input).resolve()))
    ui.kv("Output", str(output_path))
    ui.kv("Hits", str(hits))
    ui.kv("Threads", str(args.threads))
    ui.kv("Het level", str(het_level))
    ui.kv("Extended", str(args.extended_report))
    ui.kv("Write FASTA", str(args.write_fasta))
    ui.kv("Chip mode", str(args.chip))

    ui.stage("1/4", "加载系统发育树资源")
    tree_id, version = split_tree_id(args.tree)
    tree_dir = resolve_installed_tree_dir(tree_id, version, trees_dir)
    tree_bundle = load_tree_bundle(tree_dir)
    ui.kv("Tree dir", str(tree_dir))

    ui.stage("2/4", "读取并标准化输入样本")
    samples = load_samples(
        Path(args.input).resolve(),
        tree_bundle.reference_sequence,
        tree_bundle.rules,
        het_level,
        threads=args.threads,
        chip=args.chip,
    )
    ui.kv("Samples", str(len(samples)))

    ui.stage("3/4", "计算候选单倍群")
    all_hits = classify_samples(samples, tree_bundle, hits, threads=args.threads)

    ui.stage("4/4", "写出结果文件")
    write_classification_report(output_path, all_hits)
    ui.kv("Report", str(output_path))
    if args.extended_report:
        write_extended_report(output_path, all_hits)
        ui.kv("Extended report", str(extended_report_path(output_path)))
    if args.write_fasta:
        write_reconstructed_fasta(output_path, samples, tree_bundle)
        ui.kv("FASTA", str(output_path.with_suffix(".fasta")))
    ui.ok("Classification completed")
    return 0
