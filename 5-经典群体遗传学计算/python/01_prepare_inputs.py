import argparse
import gzip
import logging
from pathlib import Path

import pandas as pd


log = logging.getLogger(__name__)


def detect_separator(path: str) -> str:
    suffix = Path(path).suffix.lower()
    return "," if suffix == ".csv" else "\t"


def read_sample_table(path: str, sample_id_column: str, group_column: str) -> pd.DataFrame:
    separator = detect_separator(path)
    header_df = pd.read_csv(path, sep=separator, dtype=str)
    if sample_id_column in header_df.columns and group_column in header_df.columns:
        df = header_df.copy()
    else:
        df = pd.read_csv(path, sep=separator, header=None, dtype=str)
        if df.shape[1] < 2:
            raise ValueError("Sample table must contain at least two columns.")
        rename_map = {df.columns[0]: sample_id_column, df.columns[1]: group_column}
        df = df.rename(columns=rename_map)
    df[sample_id_column] = df[sample_id_column].astype(str).str.strip()
    df[group_column] = df[group_column].astype(str).str.strip()
    df = df[(df[sample_id_column] != "") & (df[group_column] != "")]
    df = df.drop_duplicates(subset=[sample_id_column])
    return df


def read_vcf_samples(vcf_path: str) -> list[str]:
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
    raise ValueError("VCF header line starting with #CHROM was not found.")


def run(
    input_vcf: str,
    sample_table: str,
    sample_id_column: str,
    group_column: str,
    min_group_size: int,
    matched_samples_output: str,
    group_counts_output: str,
    summary_output: str,
) -> dict[str, int]:
    vcf_path = Path(input_vcf)
    if not vcf_path.exists():
        raise FileNotFoundError(f"Input VCF not found: {input_vcf}")
    if not vcf_path.name.endswith(".vcf.gz"):
        raise ValueError("Input VCF must be bgzip-compressed with .vcf.gz suffix.")

    sample_df = read_sample_table(sample_table, sample_id_column, group_column)
    vcf_samples = read_vcf_samples(input_vcf)
    vcf_sample_set = set(vcf_samples)

    matched_df = sample_df[sample_df[sample_id_column].isin(vcf_sample_set)].copy()
    matched_df = matched_df.dropna(subset=[sample_id_column, group_column])

    group_counts = (
        matched_df.groupby(group_column, as_index=False)
        .size()
        .rename(columns={"size": "n_samples"})
        .sort_values(["n_samples", group_column], ascending=[False, True])
    )
    retained_groups = set(group_counts.loc[group_counts["n_samples"] >= min_group_size, group_column])
    matched_df = matched_df[matched_df[group_column].isin(retained_groups)].copy()
    group_counts = group_counts[group_counts[group_column].isin(retained_groups)].copy()

    Path(matched_samples_output).parent.mkdir(parents=True, exist_ok=True)
    Path(group_counts_output).parent.mkdir(parents=True, exist_ok=True)
    Path(summary_output).parent.mkdir(parents=True, exist_ok=True)

    matched_df.to_csv(matched_samples_output, sep="\t", index=False)
    group_counts.to_csv(group_counts_output, sep="\t", index=False)

    summary = pd.DataFrame(
        [
            {
                "vcf_sample_count": len(vcf_samples),
                "sample_table_rows": len(sample_df),
                "matched_sample_count": len(matched_df),
                "matched_group_count": matched_df[group_column].nunique(),
                "min_group_size": min_group_size,
                "index_present": int(Path(f"{input_vcf}.tbi").exists()),
            }
        ]
    )
    summary.to_csv(summary_output, sep="\t", index=False)

    log.info(
        "Prepared inputs: %d VCF samples, %d matched samples across %d groups",
        len(vcf_samples),
        len(matched_df),
        matched_df[group_column].nunique(),
    )
    return {
        "vcf_sample_count": len(vcf_samples),
        "matched_sample_count": len(matched_df),
        "matched_group_count": int(matched_df[group_column].nunique()),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Validate sample metadata and intersect with VCF samples.")
    parser.add_argument("--input-vcf", required=True)
    parser.add_argument("--sample-table", required=True)
    parser.add_argument("--sample-id-column", required=True)
    parser.add_argument("--group-column", required=True)
    parser.add_argument("--min-group-size", type=int, default=1)
    parser.add_argument("--matched-samples-output", required=True)
    parser.add_argument("--group-counts-output", required=True)
    parser.add_argument("--summary-output", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    args = build_parser().parse_args(argv)
    run(
        input_vcf=args.input_vcf,
        sample_table=args.sample_table,
        sample_id_column=args.sample_id_column,
        group_column=args.group_column,
        min_group_size=args.min_group_size,
        matched_samples_output=args.matched_samples_output,
        group_counts_output=args.group_counts_output,
        summary_output=args.summary_output,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
