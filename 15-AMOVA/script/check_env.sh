#!/usr/bin/env bash
# 检查运行环境依赖（Python 解释器 + 必要的 Python 包）

set -euo pipefail

CONFIG_PATH="${1:-}"
if [[ -z "$CONFIG_PATH" ]]; then
    echo "[ERROR] 用法: bash script/check_env.sh <config_path>" >&2
    exit 1
fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/load_config.sh" "$CONFIG_PATH"

PYTHON_BIN="$CFG_TOOLS_PYTHON_BIN"

# 检查 Python 解释器
if [[ ! -x "$PYTHON_BIN" ]]; then
    echo "[ERROR] Python 解释器不存在或不可执行: $PYTHON_BIN" >&2
    exit 1
fi

# 检查必要的 Python 包
REQUIRED_PACKAGES=("pysam" "numpy" "pandas" "yaml")
for pkg in "${REQUIRED_PACKAGES[@]}"; do
    if ! "$PYTHON_BIN" -c "import $pkg" 2>/dev/null; then
        echo "[ERROR] 缺少 Python 包: $pkg（请在 BigLin 环境中安装）" >&2
        exit 1
    fi
done

echo "[OK] 环境检查通过（Python: $PYTHON_BIN）"
