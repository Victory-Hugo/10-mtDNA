#!/usr/bin/env bash

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
PIPE_ENTRY="${PROJECT_DIR}/pipe/1-run-pipe.sh"

exec "$PIPE_ENTRY" "$@"
