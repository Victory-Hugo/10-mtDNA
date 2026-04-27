"""phyloMT 根目录入口脚本。"""

from __future__ import annotations

import sys
from pathlib import Path


def main() -> int:
    """支持直接执行 python3 phyloMT.py ..."""

    project_dir = Path(__file__).resolve().parent
    python_dir = project_dir / "python"
    if str(python_dir) not in sys.path:
        sys.path.insert(0, str(python_dir))

    from phylomt.cli import main as cli_main

    return cli_main()


if __name__ == "__main__":
    raise SystemExit(main())
