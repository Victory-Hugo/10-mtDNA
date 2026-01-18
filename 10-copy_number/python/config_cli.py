#!/usr/bin/env python3
import argparse
import sys
from typing import Any, Dict

import yaml


def load_config(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def run(config_path: str, key: str, default: Any = None) -> int:
    config = load_config(config_path)
    if key not in config:
        if default is None:
            sys.stderr.write(f"Missing key in config: {key}\n")
            return 1
        value = default
    else:
        value = config[key]
    if isinstance(value, (dict, list)):
        sys.stderr.write(f"Key is not a scalar value: {key}\n")
        return 1
    sys.stdout.write(str(value))
    return 0


def main() -> None:
    parser = argparse.ArgumentParser(description="Read scalar config values.")
    parser.add_argument("--config", required=True, help="Path to YAML config")
    parser.add_argument("--key", required=True, help="Config key to read")
    parser.add_argument("--default", help="Default value if key is missing")
    args = parser.parse_args()
    raise SystemExit(run(args.config, args.key, args.default))


if __name__ == "__main__":
    main()
