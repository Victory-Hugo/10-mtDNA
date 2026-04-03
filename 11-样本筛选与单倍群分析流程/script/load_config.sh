#!/usr/bin/env bash

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "请使用: source script/load_config.sh conf/Config.yaml" >&2
    exit 1
fi

load_config_yaml() {
    local config_file="$1"
    if [[ -z "$config_file" ]]; then
        echo "缺少配置文件路径" >&2
        return 1
    fi
    if [[ ! -f "$config_file" ]]; then
        echo "配置文件不存在: $config_file" >&2
        return 1
    fi

    local assignments
    assignments="$(
        awk '
        function trim(value) {
            sub(/^[[:space:]]+/, "", value)
            sub(/[[:space:]]+$/, "", value)
            return value
        }
        function shell_quote(value, escaped) {
            gsub(/\047/, "\047\\\047\047", value)
            return "\047" value "\047"
        }
        {
            gsub(/\r/, "", $0)
            if ($0 ~ /^[[:space:]]*#/ || $0 ~ /^[[:space:]]*$/) {
                next
            }
            if ($0 ~ /^[A-Za-z0-9_]+:[[:space:]]*$/) {
                split($0, parts, ":")
                section = toupper(parts[1])
                next
            }
            if ($0 ~ /^  [A-Za-z0-9_]+:[[:space:]]*/) {
                line = substr($0, 3)
                split(line, parts, ":")
                key = toupper(parts[1])
                value = substr(line, index(line, ":") + 1)
                sub(/[[:space:]]+#.*$/, "", value)
                value = trim(value)
                if ((value ~ /^".*"$/) || (value ~ /^\047.*\047$/)) {
                    value = substr(value, 2, length(value) - 2)
                }
                printf "%s_%s=%s\n", section, key, shell_quote(value)
            }
        }
        ' "$config_file"
    )" || return 1

    eval "$assignments"
}

load_config_yaml "$1"
