import argparse
import logging
import shutil
import sys
from pathlib import Path

import allel
import numcodecs
import pandas as pd


log = logging.getLogger(__name__)


def run(
    input_vcf: str,
    matched_samples: str,
    sample_id_column: str,
    zarr_path: str,
    overwrite: int,
    chunk_length: int,
    chunk_width: int,
    alt_number: int,
) -> dict[str, str]:
    input_path = Path(input_vcf)
    matched_path = Path(matched_samples)
    output_path = Path(zarr_path)

    if not input_path.exists():
        raise FileNotFoundError(f"Input VCF not found: {input_vcf}")
    if not matched_path.exists():
        raise FileNotFoundError(f"Matched sample table not found: {matched_samples}")

    sample_df = pd.read_csv(matched_path, sep="\t", dtype=str)
    sample_list = sample_df[sample_id_column].dropna().astype(str).tolist()
    if not sample_list:
        raise ValueError("No matched samples available for VCF import.")

    if output_path.exists():
        if overwrite:
            shutil.rmtree(output_path)
            log.info("Removed existing Zarr store: %s", output_path)
        else:
            log.info("Zarr store already exists and overwrite=0, skipping import")
            return {"zarr_path": str(output_path), "status": "skipped"}

    output_path.parent.mkdir(parents=True, exist_ok=True)
    compressor = numcodecs.Blosc(cname="zstd", clevel=5, shuffle=numcodecs.Blosc.BITSHUFFLE)
    fields = ["samples", "variants/CHROM", "variants/POS", "variants/REF", "variants/ALT", "calldata/GT"]

    allel.vcf_to_zarr(
        str(input_path),
        str(output_path),
        group="/",
        fields=fields,
        samples=sample_list,
        overwrite=True,
        log=sys.stdout,
        compressor=compressor,
        chunk_length=chunk_length,
        chunk_width=chunk_width,
        alt_number=alt_number,
    )
    log.info("Imported VCF subset into Zarr: %s", output_path)
    return {"zarr_path": str(output_path), "status": "created"}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Import selected VCF samples into a Zarr store.")
    parser.add_argument("--input-vcf", required=True)
    parser.add_argument("--matched-samples", required=True)
    parser.add_argument("--sample-id-column", required=True)
    parser.add_argument("--zarr-path", required=True)
    parser.add_argument("--overwrite", type=int, default=0)
    parser.add_argument("--chunk-length", type=int, default=2048)
    parser.add_argument("--chunk-width", type=int, default=128)
    parser.add_argument("--alt-number", type=int, default=4)
    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    args = build_parser().parse_args(argv)
    run(
        input_vcf=args.input_vcf,
        matched_samples=args.matched_samples,
        sample_id_column=args.sample_id_column,
        zarr_path=args.zarr_path,
        overwrite=args.overwrite,
        chunk_length=args.chunk_length,
        chunk_width=args.chunk_width,
        alt_number=args.alt_number,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
