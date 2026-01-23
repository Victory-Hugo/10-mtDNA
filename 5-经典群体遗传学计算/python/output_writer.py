from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable
import logging

import numpy as np
import pandas as pd

from metrics import SfsResult

logger = logging.getLogger(__name__)


def save_population_metrics(
    output_dir: Path,
    pi_results: Dict[str, float],
    theta_w_results: Dict[str, float],
    tajima_d_results: Dict[str, float],
    sfs_folded: Dict[str, SfsResult],
    sfs_unfolded: Dict[str, SfsResult],
    population_counts: pd.DataFrame,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("写入群体指标: %s", output_dir)

    population_counts.to_csv(output_dir / "population_counts.csv", index=False)

    pd.DataFrame({"population": list(pi_results.keys()), "pi": list(pi_results.values())}).to_csv(
        output_dir / "pi.csv", index=False
    )
    pd.DataFrame({"population": list(theta_w_results.keys()), "theta_w": list(theta_w_results.values())}).to_csv(
        output_dir / "theta_w.csv", index=False
    )
    pd.DataFrame({"population": list(tajima_d_results.keys()), "tajima_d": list(tajima_d_results.values())}).to_csv(
        output_dir / "tajima_d.csv", index=False
    )

    def _iter_sfs_rows(
        sfs_data: Dict[str, SfsResult],
        folded: bool,
    ) -> Iterable[dict]:
        for pop, payload in sfs_data.items():
            values = payload["sfs"]
            n_samples = payload["n_samples"]
            n_segregating = payload["n_segregating"]
            if values is None:
                yield {
                    "population": pop,
                    "folded": folded,
                    "bin": np.nan,
                    "count": np.nan,
                    "n_samples": n_samples,
                    "n_segregating": n_segregating,
                }
            else:
                for idx, count in enumerate(values):
                    yield {
                        "population": pop,
                        "folded": folded,
                        "bin": idx,
                        "count": float(count),
                        "n_samples": n_samples,
                        "n_segregating": n_segregating,
                    }

    sfs_long = pd.DataFrame(
        list(_iter_sfs_rows(sfs_folded, True)) + list(_iter_sfs_rows(sfs_unfolded, False))
    )
    sfs_long.to_csv(output_dir / "sfs_long.csv", index=False)
    logger.info("群体指标写入完成: %s", output_dir)


def save_pairwise_metrics(
    output_dir: Path,
    fst: pd.DataFrame,
    dxy: pd.DataFrame,
    fst_p: pd.DataFrame | None = None,
    dxy_p: pd.DataFrame | None = None,
    pi_bootstrap: pd.DataFrame | None = None,
    theta_w_bootstrap: pd.DataFrame | None = None,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("写入成对指标: %s", output_dir)
    fst.to_csv(output_dir / "fst.csv")
    dxy.to_csv(output_dir / "dxy.csv")
    if fst_p is not None:
        fst_p.to_csv(output_dir / "fst_pvalue.csv")
    if dxy_p is not None:
        dxy_p.to_csv(output_dir / "dxy_pvalue.csv")
    if pi_bootstrap is not None:
        pi_bootstrap.to_csv(output_dir / "pi_bootstrap.csv", index=False)
    if theta_w_bootstrap is not None:
        theta_w_bootstrap.to_csv(output_dir / "theta_w_bootstrap.csv", index=False)
    logger.info("成对指标写入完成: %s", output_dir)
