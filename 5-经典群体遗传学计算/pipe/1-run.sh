#!/usr/bin/env bash
# 脚本名称: 1-run.sh

# 描述:
# 本脚本用于自动化执行经典群体遗传学相关的分析流程。通过读取配置文件，自动设置参数、日志记录、样本表过滤及主分析流程。支持断点续跑、日志追踪、可选的样本表过滤、复杂变异跳过、显著性分析、随机种子设置及Tajima分析参数配置。

# 用法:
#   ./1-run.sh /path/to/conf/1-run.conf

# 主要流程:
#   1. 检查并加载配置文件，初始化日志。
#   2. 若检测到已成功完成的日志且未强制运行，则跳过执行。
#   3. 可选：根据配置过滤样本表。
#   4. 调用主分析脚本population_genetics.py，传递所有参数。
#   5. 成功完成后记录SUCCESS。

# 注意事项:
#   - 需提前配置好conda环境及相关依赖。
#   - 配置文件需包含所有必要参数。
#   - 日志文件自动追加，便于追踪和断点续跑。

set -euo pipefail

CONF_PATH="${1:-}"
if [[ -z "$CONF_PATH" ]]; then
  echo "用法: 1-run.sh /path/to/conf/1-run.conf" >&2
  exit 1
fi

source "$CONF_PATH"

LOG_FILE="${LOG_DIR}/1-run.sh.log"
mkdir -p "$LOG_DIR"

exec > >(tee -a "$LOG_FILE") 2>&1

if [[ -f "$LOG_FILE" && "${FORCE}" != "true" ]]; then
  if grep -q "SUCCESS" "$LOG_FILE"; then
    echo "检测到已成功完成的日志记录，跳过执行"
    exit 0
  fi
fi

trap 'echo "FAILED $(date) code=$?"' ERR

{
  echo "START $(date)"
  echo "CONF=$CONF_PATH"
  echo "VCF_PATH=$VCF_PATH"
  echo "SAMPLE_TABLE_PATH=$SAMPLE_TABLE_PATH"
  echo "FILTERED_SAMPLE_TABLE_PATH=$FILTERED_SAMPLE_TABLE_PATH"
  echo "OUTPUT_DIR=$OUTPUT_DIR"
  echo "GROUP_COLS=$GROUP_COLS"
  echo "ID_COL=$ID_COL"
  echo "CHROM_NAME=$CHROM_NAME"
  echo "SKIP_COMPLEX_VARIANTS=$SKIP_COMPLEX_VARIANTS"
  echo "USE_FILTERED_TABLE=$USE_FILTERED_TABLE"
  echo "ENABLE_SIGNIFICANCE=$ENABLE_SIGNIFICANCE"
  echo "PERMUTATION_N=$PERMUTATION_N"
  echo "BOOTSTRAP_N=$BOOTSTRAP_N"
  echo "RANDOM_SEED=$RANDOM_SEED"
  echo "CONDA_ENV=$CONDA_ENV"
}

FILTER_SCRIPT="$PROJECT_DIR/python/filter_sample_table.py"
MAIN_SCRIPT="$PROJECT_DIR/python/population_genetics.py"

RUN_SAMPLE_TABLE="$SAMPLE_TABLE_PATH"
if [[ "$USE_FILTERED_TABLE" == "true" ]]; then
  echo "################ 运行样本表过滤步骤 ################"
  conda run -n "$CONDA_ENV" python "$FILTER_SCRIPT" \
    --input "$SAMPLE_TABLE_PATH" \
    --output "$FILTERED_SAMPLE_TABLE_PATH" \
    --group-cols "$GROUP_COLS"
  RUN_SAMPLE_TABLE="$FILTERED_SAMPLE_TABLE_PATH"
fi

SKIP_FLAG=""
if [[ "$SKIP_COMPLEX_VARIANTS" == "true" ]]; then
  SKIP_FLAG="--skip-complex-variants"
fi

SIGNIF_FLAG=""
if [[ "${ENABLE_SIGNIFICANCE:-false}" == "true" ]]; then
  SIGNIF_FLAG="--enable-significance"
fi

SEED_FLAG=""
if [[ -n "${RANDOM_SEED:-}" ]]; then
  SEED_FLAG="--random-seed ${RANDOM_SEED}"
fi
echo "###########################################"
echo "############### 运行主分析 ################"
echo "###########################################"
echo "[请注意检查日志输出以确认参数设置正确]"
echo "[请注意检查每个群体数量是否均≥2否则报错]"
echo "###########################################"
conda run -n "$CONDA_ENV" python "$MAIN_SCRIPT" \
  --vcf "$VCF_PATH" \
  --sample-table "$RUN_SAMPLE_TABLE" \
  --group-cols "$GROUP_COLS" \
  --id-col "$ID_COL" \
  --output-dir "$OUTPUT_DIR" \
  --chrom "$CHROM_NAME" \
  --permutation-n "${PERMUTATION_N:-1000}" \
  --bootstrap-n "${BOOTSTRAP_N:-1000}" \
  $SEED_FLAG \
  $SIGNIF_FLAG \
  $SKIP_FLAG

echo "SUCCESS $(date)"
