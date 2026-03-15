#!/usr/bin/env bash
set -euo pipefail

require_cmd() {
    local cmd="$1"
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found: $cmd" >&2
        exit 1
    fi
}

require_cmd bash
require_cmd awk
require_cmd sed
require_cmd grep
require_cmd conda

echo "[OK] Base environment check passed."
