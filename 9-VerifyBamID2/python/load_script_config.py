#!/usr/bin/env python3
import argparse
from pathlib import Path

import yaml


def load_config(config_path: Path) -> dict:
    with config_path.open("r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh) or {}
    return cfg.get("script", {})


def to_shell_assignments(cfg: dict) -> str:
    lines = []
    for key, value in cfg.items():
        if isinstance(value, bool):
            val = "true" if value else "false"
        else:
            val = str(value)
        lines.append(f"{key}={val!r}")
    return "\n".join(lines)


def run(config: str) -> str:
    cfg = load_config(Path(config))
    return to_shell_assignments(cfg)


def main() -> None:
    parser = argparse.ArgumentParser(description="Load script config as shell assignments.")
    parser.add_argument("--config", required=True, help="Config YAML path")
    args = parser.parse_args()
    print(run(args.config))


if __name__ == "__main__":
    main()
