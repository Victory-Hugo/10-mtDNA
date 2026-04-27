#!/usr/bin/env bash
set -euo pipefail

CONFIG_FILE="${1:?Usage: source script/load_config.sh conf/Config.yaml}"

if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "配置文件不存在: $CONFIG_FILE" >&2
  return 1 2>/dev/null || exit 1
fi

# 解析本项目使用的简单 YAML 子集：顶层 section + 二级 key + 标量值。
# 变量名格式为 SECTION_KEY，例如 project.base_dir -> PROJECT_BASE_DIR。
while IFS= read -r raw_line || [[ -n "$raw_line" ]]; do
  line="${raw_line%%#*}"
  [[ -z "${line//[[:space:]]/}" ]] && continue

  if [[ "$line" =~ ^([A-Za-z0-9_]+):[[:space:]]*$ ]]; then
    section="${BASH_REMATCH[1]}"
    continue
  fi

  if [[ "$line" =~ ^[[:space:]]{2}([A-Za-z0-9_]+):[[:space:]]*(.*)$ ]]; then
    key="${BASH_REMATCH[1]}"
    value="${BASH_REMATCH[2]}"
    value="${value#"${value%%[![:space:]]*}"}"
    value="${value%"${value##*[![:space:]]}"}"
    value="${value%\"}"
    value="${value#\"}"
    value="${value%\'}"
    value="${value#\'}"

    var_name="$(printf '%s_%s' "$section" "$key" | tr '[:lower:]' '[:upper:]')"
    printf -v "$var_name" '%s' "$value"
    export "$var_name"
  fi
done < "$CONFIG_FILE"

resolve_path() {
  local path_value="${1:?path required}"
  if [[ "$path_value" = /* ]]; then
    printf '%s\n' "$path_value"
  else
    printf '%s/%s\n' "$PROJECT_BASE_DIR" "$path_value"
  fi
}
