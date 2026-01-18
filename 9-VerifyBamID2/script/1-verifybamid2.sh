#!/usr/bin/env bash
set -euo pipefail

#! conda install -c bioconda verifybamid2

base_dir="/mnt/f/Onedrive/文档（科研）/脚本/Download/10-mtDNA/9-VerifyBamID2/"
config_path="${base_dir}/conf/verifybamid2.yml"
config_loader="${base_dir}/python/load_script_config.py"

if [[ ! -f "$config_loader" ]]; then
  echo "Config loader not found: $config_loader" >&2
  exit 1
fi

args=("$@")
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      config_path="$2"
      shift 2
      ;;
    *)
      shift 1
      ;;
  esac
done

eval "$("python3" "$config_loader" --config "$config_path")"

# 保存原始参数，处理空参数情况
if [[ $# -gt 0 ]]; then
    args=("$@")
else
    args=()
fi


while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      config_path="$2"
      shift 2
      ;;
    --list)
      list_path="$2"
      shift 2
      ;;
    --output-dir)
      output_dir="$2"
      shift 2
      ;;
    --log-dir)
      log_dir="$2"
      shift 2
      ;;
    --tmp-dir)
      tmp_dir="$2"
      shift 2
      ;;
    --python)
      python_bin="$2"
      shift 2
      ;;
    --jobs)
      jobs="$2"
      shift 2
      ;;
    --force)
      force="true"
      shift 1
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

mkdir -p "$log_dir" "$tmp_dir" "$output_dir"
log_file="${log_dir}/1-verifybamid2.log"

{
  echo "$(date -Is)\tSTART"
  echo "CONFIG=${config_path}"
  echo "LIST=${list_path}"
  echo "OUTPUT_DIR=${output_dir}"
  echo "LOG_DIR=${log_dir}"
  echo "TMP_DIR=${tmp_dir}"
  echo "JOBS=${jobs}"

  if [[ ! -f "$config_path" ]]; then
    echo "Config not found: $config_path" >&2
    echo "$(date -Is)\tSTATUS=FAIL\tREASON=config_not_found"
    exit 1
  fi

  if [[ ! -f "$list_path" ]]; then
    echo "List not found: $list_path" >&2
    echo "$(date -Is)\tSTATUS=FAIL\tREASON=list_not_found"
    exit 1
  fi

  task_script="${base_dir}/python/verifybamid2_task.py"
  if [[ ! -f "$task_script" ]]; then
    echo "Task script not found: $task_script" >&2
    echo "$(date -Is)\tSTATUS=FAIL\tREASON=task_script_not_found"
    exit 1
  fi

  force_flag=""
  if [[ "$force" == "true" ]]; then
    force_flag="--force"
  fi

  if command -v parallel >/dev/null 2>&1; then
    parallel --colsep '\t' --jobs "$jobs" --joblog "${log_dir}/1-verifybamid2.parallel.log" \
      "$python_bin" "$task_script" \
      --sample {1} --bam {2} \
      --config "$config_path" \
      --output-dir "$output_dir" \
      --log-dir "$log_dir" \
      --tmp-dir "$tmp_dir" \
      $force_flag \
      :::: "$list_path"
  else
    awk -F'\t' 'NF>=2 {print $1"\t"$2}' "$list_path" | \
      xargs -P "$jobs" -L 1 -I {} bash -c \
      'sample="$(printf "%s" "{}" | cut -f1)"; bam="$(printf "%s" "{}" | cut -f2)"; \
      '"$python_bin"' '"$task_script"' --sample "$sample" --bam "$bam" \
      --config '"$config_path"' --output-dir '"$output_dir"' \
      --log-dir '"$log_dir"' --tmp-dir '"$tmp_dir"' '"$force_flag"''
  fi

  echo "$(date -Is)\tSTATUS=SUCCESS"
} >> "$log_file" 2>&1
