#!/usr/bin/env bash
# 读取简单 YAML 配置，并导出为 SECTION_KEY 形式的变量。
# 支持范围：顶层 section + 二级 key + 标量值；不支持数组和深层嵌套。

strip_yaml_comment() {
    local input="$1"
    local output=""
    local ch=""
    local in_single=0
    local in_double=0
    local i

    for ((i = 0; i < ${#input}; i++)); do
        ch="${input:i:1}"
        if [[ "$ch" == "'" && "$in_double" -eq 0 ]]; then
            if [[ "$in_single" -eq 0 ]]; then
                in_single=1
            else
                in_single=0
            fi
            output+="$ch"
        elif [[ "$ch" == '"' && "$in_single" -eq 0 ]]; then
            if [[ "$in_double" -eq 0 ]]; then
                in_double=1
            else
                in_double=0
            fi
            output+="$ch"
        elif [[ "$ch" == "#" && "$in_single" -eq 0 && "$in_double" -eq 0 ]]; then
            break
        else
            output+="$ch"
        fi
    done
    printf '%s\n' "$output"
}

load_config() {
    local config_file="$1"
    local section=""
    local line key value var_name

    if [[ ! -f "$config_file" ]]; then
        echo "[load_config] 配置文件不存在: $config_file" >&2
        return 1
    fi

    while IFS= read -r line || [[ -n "$line" ]]; do
        line="$(strip_yaml_comment "$line")"
        line="$(printf '%s' "$line" | sed 's/[[:space:]]*$//')"
        [[ -z "$line" ]] && continue

        if [[ "$line" =~ ^([A-Za-z0-9_]+):[[:space:]]*$ ]]; then
            section="${BASH_REMATCH[1]}"
            continue
        fi

        if [[ "$line" =~ ^[[:space:]]{2}([A-Za-z0-9_]+):[[:space:]]*(.*)$ ]]; then
            key="${BASH_REMATCH[1]}"
            value="${BASH_REMATCH[2]}"
            value="$(printf '%s' "$value" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')"
            value="${value%\"}"
            value="${value#\"}"
            value="${value%\'}"
            value="${value#\'}"
            var_name="$(printf '%s_%s' "$section" "$key" | tr '[:lower:]' '[:upper:]')"
            printf -v "$var_name" '%s' "$value"
            export "$var_name"
        fi
    done < "$config_file"
}
