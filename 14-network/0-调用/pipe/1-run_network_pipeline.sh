#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'USAGE'
用法:
  bash pipe/1-run_network_pipeline.sh [-c conf/Config.yaml]

说明:
  读取配置文件，整理 FASTA/metadata，并按 runtime.software 运行 fastHaN、McAN 或二者。
USAGE
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIG_FILE="$DEFAULT_BASE_DIR/conf/Config.yaml"

while getopts ":c:h" opt; do
    case "$opt" in
        c) CONFIG_FILE="$OPTARG" ;;
        h) usage; exit 0 ;;
        *) usage >&2; exit 1 ;;
    esac
done

source "$DEFAULT_BASE_DIR/script/load_config.sh"
load_config "$CONFIG_FILE"

resolve_path() {
    local path="$1"
    if [[ "$path" = /* ]]; then
        printf '%s\n' "$path"
    else
        printf '%s/%s\n' "$PROJECT_BASE_DIR" "$path"
    fi
}

require_value() {
    local name="$1"
    local value="${!name:-}"
    if [[ -z "$value" ]]; then
        echo "[ERROR] 配置缺少字段: $name" >&2
        exit 1
    fi
}

require_file() {
    local label="$1"
    local path="$2"
    if [[ ! -f "$path" ]]; then
        echo "[ERROR] $label 不存在: $path" >&2
        exit 1
    fi
}

require_executable() {
    local label="$1"
    local path="$2"
    if [[ ! -x "$path" ]]; then
        echo "[ERROR] $label 不存在或不可执行: $path" >&2
        exit 1
    fi
}

run_cmd() {
    echo "[CMD] $*" | tee -a "$LOG_FILE"
    "$@" 2>&1 | tee -a "$LOG_FILE"
}

append_summary() {
    printf '%s\n' "$*" >> "$SUMMARY_FILE"
}

for required in \
    PROJECT_BASE_DIR PROJECT_INPUT_FASTA PROJECT_INPUT_META PROJECT_SAMPLE_ID_COLUMN PROJECT_REFERENCE_ID \
    PATHS_PYTHON_DIR PATHS_TEMP_PREPARE_DIR PATHS_TEMP_FASTHAN_DIR PATHS_TEMP_MCAN_DIR \
    PATHS_PREPARE_REPORT_DIR PATHS_FASTHAN_TABLE_DIR PATHS_MCAN_TABLE_DIR PATHS_MCAN_REPORT_DIR \
    PATHS_SUMMARY_REPORT_DIR PATHS_LOG_DIR TOOLS_PYTHON_BIN TOOLS_FASTHAN_BIN TOOLS_MCAN_BIN \
    RUNTIME_SOFTWARE RUNTIME_THREADS RUNTIME_OVERWRITE RUNTIME_FASTHAN_ALGORITHMS RUNTIME_VISUAL_COLUMNS \
    RUNTIME_FASTHAN_COLOR_PALETTE RUNTIME_MCAN_COUNTRY_COLUMN RUNTIME_MCAN_STATE_COLUMN \
    RUNTIME_MCAN_CITY_COLUMN RUNTIME_MCAN_DEFAULT_DATE; do
    require_value "$required"
done

INPUT_FASTA="$(resolve_path "$PROJECT_INPUT_FASTA")"
INPUT_META="$(resolve_path "$PROJECT_INPUT_META")"
PYTHON_DIR="$(resolve_path "$PATHS_PYTHON_DIR")"
TEMP_PREPARE_DIR="$(resolve_path "$PATHS_TEMP_PREPARE_DIR")"
TEMP_FASTHAN_DIR="$(resolve_path "$PATHS_TEMP_FASTHAN_DIR")"
TEMP_MCAN_DIR="$(resolve_path "$PATHS_TEMP_MCAN_DIR")"
PREPARE_REPORT_DIR="$(resolve_path "$PATHS_PREPARE_REPORT_DIR")"
FASTHAN_TABLE_DIR="$(resolve_path "$PATHS_FASTHAN_TABLE_DIR")"
MCAN_TABLE_DIR="$(resolve_path "$PATHS_MCAN_TABLE_DIR")"
MCAN_REPORT_DIR="$(resolve_path "$PATHS_MCAN_REPORT_DIR")"
SUMMARY_REPORT_DIR="$(resolve_path "$PATHS_SUMMARY_REPORT_DIR")"
LOG_DIR="$(resolve_path "$PATHS_LOG_DIR")"
PYTHON_BIN="$(resolve_path "$TOOLS_PYTHON_BIN")"
FASTHAN_BIN="$(resolve_path "$TOOLS_FASTHAN_BIN")"
MCAN_BIN="$(resolve_path "$TOOLS_MCAN_BIN")"

mkdir -p "$TEMP_PREPARE_DIR" "$TEMP_FASTHAN_DIR" "$TEMP_MCAN_DIR" \
    "$PREPARE_REPORT_DIR" "$FASTHAN_TABLE_DIR" "$MCAN_TABLE_DIR" "$MCAN_REPORT_DIR" \
    "$SUMMARY_REPORT_DIR" "$LOG_DIR"

LOG_FILE="$LOG_DIR/pipeline.log"
SUMMARY_FILE="$SUMMARY_REPORT_DIR/pipeline_summary.txt"
if [[ "$RUNTIME_OVERWRITE" == "true" ]]; then
    : > "$LOG_FILE"
    : > "$SUMMARY_FILE"
else
    if [[ -s "$SUMMARY_FILE" ]]; then
        echo "[ERROR] summary 已存在且 overwrite=false: $SUMMARY_FILE" >&2
        exit 1
    fi
fi

require_file "输入 FASTA" "$INPUT_FASTA"
require_file "输入 metadata" "$INPUT_META"
require_executable "Python" "$PYTHON_BIN"

case "$RUNTIME_SOFTWARE" in
    fasthan|mcan|both) ;;
    *) echo "[ERROR] runtime.software 只能是 fasthan/mcan/both: $RUNTIME_SOFTWARE" >&2; exit 1 ;;
esac

if [[ "$RUNTIME_SOFTWARE" == "fasthan" || "$RUNTIME_SOFTWARE" == "both" ]]; then
    require_executable "fastHaN" "$FASTHAN_BIN"
fi
if [[ "$RUNTIME_SOFTWARE" == "mcan" || "$RUNTIME_SOFTWARE" == "both" ]]; then
    require_executable "McAN" "$MCAN_BIN"
fi

FILTERED_FASTA="$TEMP_PREPARE_DIR/filtered_sequences.fasta"
FASTHAN_PHYLIP="$TEMP_PREPARE_DIR/fasthan_input.phy"
QC_REPORT="$PREPARE_REPORT_DIR/input_qc.tsv"

append_summary "# network pipeline summary"
append_summary "config: $CONFIG_FILE"
append_summary "input_fasta: $INPUT_FASTA"
append_summary "input_meta: $INPUT_META"
append_summary "software: $RUNTIME_SOFTWARE"
append_summary "reference_id: $PROJECT_REFERENCE_ID"
append_summary "fasthan_color_palette: $RUNTIME_FASTHAN_COLOR_PALETTE"
append_summary ""

run_cmd "$PYTHON_BIN" "$PYTHON_DIR/1-1-prepare_input.py" \
    --fasta "$INPUT_FASTA" \
    --metadata "$INPUT_META" \
    --sample-id-column "$PROJECT_SAMPLE_ID_COLUMN" \
    --filtered-fasta "$FILTERED_FASTA" \
    --phylip "$FASTHAN_PHYLIP" \
    --qc-report "$QC_REPORT"

run_fasthan() {
    local algorithms_csv="$1"
    local algorithm out_prefix out_dir algorithm_upper args_var args_template args_string
    local -a algorithms fasthan_extra_args
    IFS=',' read -r -a algorithms <<< "$algorithms_csv"
    for algorithm in "${algorithms[@]}"; do
        algorithm="$(printf '%s' "$algorithm" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')"
        [[ -z "$algorithm" ]] && continue
        case "$algorithm" in
            original_tcs|modified_tcs|mjn|msn) ;;
            *)
                echo "[ERROR] 不支持的 fastHaN 算法: $algorithm" >&2
                exit 1
                ;;
        esac

        algorithm_upper="$(printf '%s' "$algorithm" | tr '[:lower:]' '[:upper:]')"
        args_var="RUNTIME_FASTHAN_${algorithm_upper}_ARGS"
        if [[ ! ${!args_var+x} ]]; then
            echo "[ERROR] 配置缺少 ${args_var}，请在 runtime 中为 $algorithm 设置额外参数" >&2
            exit 1
        fi
        args_template="${!args_var}"
        args_string="${args_template//\{threads\}/$RUNTIME_THREADS}"
        if [[ " $args_string " =~ [[:space:]]-i[[:space:]] || " $args_string " =~ [[:space:]]-o[[:space:]] ]]; then
            echo "[ERROR] $args_var 只允许填写算法额外参数，不要包含 -i 或 -o" >&2
            exit 1
        fi
        read -r -a fasthan_extra_args <<< "$args_string"

        out_dir="$FASTHAN_TABLE_DIR/$algorithm"
        mkdir -p "$out_dir"
        out_prefix="$out_dir/network"
        run_cmd "$FASTHAN_BIN" "$algorithm" -i "$FASTHAN_PHYLIP" "${fasthan_extra_args[@]}" -o "$out_prefix"

        if [[ -n "$RUNTIME_VISUAL_COLUMNS" ]]; then
            run_cmd "$PYTHON_BIN" "$PYTHON_DIR/1-2-generate_fasthan_visual_config.py" \
                --json-file "$out_prefix.json" \
                --metadata "$INPUT_META" \
                --sample-id-column "$PROJECT_SAMPLE_ID_COLUMN" \
                --columns "$RUNTIME_VISUAL_COLUMNS" \
                --output-dir "$out_dir" \
                --prefix "$algorithm" \
                --color-palette "$RUNTIME_FASTHAN_COLOR_PALETTE"
        fi
        append_summary "fastHaN_$algorithm: $out_dir"
        append_summary "fastHaN_${algorithm}_args: $args_string"
    done
}

run_mcan() {
    local mutation_file="$TEMP_MCAN_DIR/mcan.mutation"
    local meta_file="$TEMP_MCAN_DIR/mcan.meta"
    local sitemask_file="$TEMP_MCAN_DIR/siteMask"
    local convert_report="$MCAN_REPORT_DIR/mcan_convert_report.tsv"

    run_cmd "$PYTHON_BIN" "$PYTHON_DIR/1-3-convert_fasta_to_mcan.py" \
        --fasta "$FILTERED_FASTA" \
        --metadata "$INPUT_META" \
        --sample-id-column "$PROJECT_SAMPLE_ID_COLUMN" \
        --reference-id "$PROJECT_REFERENCE_ID" \
        --mutation-output "$mutation_file" \
        --meta-output "$meta_file" \
        --sitemask-output "$sitemask_file" \
        --report-output "$convert_report" \
        --country-column "$RUNTIME_MCAN_COUNTRY_COLUMN" \
        --state-column "$RUNTIME_MCAN_STATE_COLUMN" \
        --city-column "$RUNTIME_MCAN_CITY_COLUMN" \
        --default-date "$RUNTIME_MCAN_DEFAULT_DATE"

    run_cmd "$MCAN_BIN" \
        --mutation "$mutation_file" \
        --meta "$meta_file" \
        --sitemask "$sitemask_file" \
        --outDir "$MCAN_TABLE_DIR" \
        --oJSON \
        --oGraphML \
        --oTSV \
        --nthread "$RUNTIME_THREADS"
    append_summary "McAN: $MCAN_TABLE_DIR"
}

if [[ "$RUNTIME_SOFTWARE" == "fasthan" || "$RUNTIME_SOFTWARE" == "both" ]]; then
    run_fasthan "$RUNTIME_FASTHAN_ALGORITHMS"
fi

if [[ "$RUNTIME_SOFTWARE" == "mcan" || "$RUNTIME_SOFTWARE" == "both" ]]; then
    run_mcan
fi

append_summary ""
append_summary "qc_report: $QC_REPORT"
append_summary "log: $LOG_FILE"
echo "[DONE] network pipeline finished. Summary: $SUMMARY_FILE"
