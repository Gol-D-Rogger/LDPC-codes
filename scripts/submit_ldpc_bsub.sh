#!/usr/bin/env bash
# 批量提交 LDPC 仿真（LSF bsub）
# 需求:
#   - SNR 可配（列表或区间）
#   - matrix 路径可配，性能输出路径可配
#   - config 文件可配 (xxx.cnfg)
#   - 自动读取矩阵编号（从文件名提取数字），投递每个矩阵
#
# 用法示例:
#   bash scripts/submit_ldpc_bsub.sh \
#     --exec gen4_ldpc_sim/ssd_fc \
#     --config gen4_ldpc_sim/config/sdec.cnfg \
#     --matrix-dir gen4_ldpc_sim/matrix \
#     --out-dir runs/149x129 \
#     --snr "3.0 3.2 3.4 3.6" \
#     --queue regr_q
#
#   # 或用区间生成 SNR（start:step:end）
#   bash scripts/submit_ldpc_bsub.sh ... --snr-seq 3.0:0.1:3.4
#
#   # 仅预览（不提交）
#   bash scripts/submit_ldpc_bsub.sh ... --dry-run

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKSPACE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# 路径转换为绝对路径（相对仓库根）
to_abs_path() {
  local path="$1"
  if [[ -z "$path" ]]; then
    echo ""
  elif [[ "$path" = /* ]]; then
    echo "$path"
  else
    echo "$WORKSPACE_ROOT/$path"
  fi
}

# 默认值
EXEC="gen4_ldpc_sim/ssd_fc_dq"
CONFIG="gen4_ldpc_sim/config/sdec_dq.cnfg"
MATRIX_DIR="gen4_ldpc_sim/matrix"
OUT_DIR="runs/ldpc_batch"
QUEUE="${LSF_QUEUE:-regr_q}"
JOB_PREFIX="ldpc"
CWD=""          # 默认不指定，由命令内自己 cd
SNR_LIST=""
SNR_SEQ=""
DRY_RUN=0

# 生成矩阵相关（可选）
GEN_N="5"
GEN_K="20"
GEN_COUNT="1"
GEN_QUEUE=""
GEN_MODE="lsf"   # lsf|local
GEN_BIN="GenLDPC/ldpc_gen"
BUILD=0
GENERATED=0
declare -a NEW_MATRIX_IDS=()

# 重命名相关（可选）
RENAME_ENABLE=1
# 兼容 test_ldpc_sim.sh 的标准重命名

# 是否显式指定了 --matrix-dir，用于在生成矩阵时决定默认扫描目录
MATRIX_DIR_SET=0

usage() {
  echo "Usage: $0 [options]"
  echo "  --exec PATH            模拟器可执行文件 (默认: ${EXEC})"
  echo "  --config PATH          配置文件 .cnfg (默认: ${CONFIG})"
  echo "  --matrix-dir DIR       矩阵目录 (默认: ${MATRIX_DIR})"
  echo "  --out-dir DIR          性能输出目录 (默认: ${OUT_DIR})"
  echo "  --queue NAME           LSF 队列 (默认: ${QUEUE})"
  echo "  --job-prefix STR       任务名前缀 (默认: ${JOB_PREFIX})"
  echo "  --cwd DIR              bsub -cwd（如不设，则命令内部 cd 到 exec 所在目录）"
  echo "  --snr \"3.0 3.2 ...\"    SNR 列表（空格分隔，二选一）"
  echo "  --snr-seq A:B:C        SNR 区间 (起:步:止)，如 3.0:0.1:3.4（二选一）"
  echo "  --build                自动编译 GenLDPC 与 gen4_ldpc_sim（可选）"
  echo "  --gen-bin PATH         生成器 ldpc_gen 路径 (默认: ${GEN_BIN})"
  echo "  --gen-n N              生成矩阵：N（可选）"
  echo "  --gen-k K              生成矩阵：K（可选）"
  echo "  --gen-count NUM        生成矩阵：提交数量（phase=1..NUM）（可选）"
  echo "  --gen-queue NAME       生成矩阵使用的队列（缺省沿用 --queue）"
  echo "  --gen-mode MODE        生成模式：lsf|local（默认 lsf）"
  echo "  --gen-out-base DIR     生成输出根目录（传递给 ldpc_gen 的第4参，同时导出 GENLDPC_OUT_DIR）"
  echo "  --rename-enable        生成后启用标准重命名（默认开启）"
  echo "  --no-rename            跳过生成后的重命名步骤"
  echo "  --sim-matrix-dir DIR   仿真读取矩阵目录（传给 ssd_fc 第6参，并导出 LDPC_MATRIX_DIR）"
  echo "  --dry-run              只打印命令，不提交"
  echo "  --help                 显示帮助"
}

# 解析参数
while [[ $# -gt 0 ]]; do
  case "$1" in
    --exec) EXEC="$2"; shift 2;;
    --config) CONFIG="$2"; shift 2;;
    --matrix-dir) MATRIX_DIR="$2"; MATRIX_DIR_SET=1; shift 2;;
    --out-dir) OUT_DIR="$2"; shift 2;;
    --queue) QUEUE="$2"; shift 2;;
    --job-prefix) JOB_PREFIX="$2"; shift 2;;
    --cwd) CWD="$2"; shift 2;;
    --snr) SNR_LIST="$2"; shift 2;;
    --snr-seq) SNR_SEQ="$2"; shift 2;;
    --build) BUILD=1; shift;;
    --gen-bin) GEN_BIN="$2"; shift 2;;
    --gen-n) GEN_N="$2"; shift 2;;
    --gen-k) GEN_K="$2"; shift 2;;
    --gen-count) GEN_COUNT="$2"; shift 2;;
    --gen-queue) GEN_QUEUE="$2"; shift 2;;
    --gen-mode) GEN_MODE="$2"; shift 2;;
    --gen-out-base) GEN_OUT_BASE="$2"; shift 2;;
    --rename-enable) RENAME_ENABLE=1; shift;;
    --no-rename) RENAME_ENABLE=0; shift;;
    --sim-matrix-dir) SIM_MATRIX_DIR="$2"; shift 2;;
    --dry-run) DRY_RUN=1; shift;;
    --help|-h) usage; exit 0;;
    *) echo "Unknown option: $1"; usage; exit 1;;
  esac
done

# 路径归一化
EXEC_ABS="$(to_abs_path "$EXEC")"
CONFIG_ABS="$(to_abs_path "$CONFIG")"
GEN_BIN_ABS="$(to_abs_path "$GEN_BIN")"
if [[ -n "${GEN_OUT_BASE:-}" ]]; then
  GEN_OUT_BASE_ABS="$(to_abs_path "$GEN_OUT_BASE")"
else
  GEN_OUT_BASE_ABS="$(to_abs_path "GenLDPC/output")"
fi
EXEC="$EXEC_ABS"
CONFIG="$CONFIG_ABS"
GEN_BIN="$GEN_BIN_ABS"
GEN_OUT_BASE="$GEN_OUT_BASE_ABS"
if [[ -n "${SIM_MATRIX_DIR:-}" ]]; then
  SIM_MATRIX_DIR="$(to_abs_path "$SIM_MATRIX_DIR")"
fi
if [[ -n "${MATRIX_DIR:-}" ]]; then
  MATRIX_DIR="$(to_abs_path "$MATRIX_DIR")"
fi
OUT_DIR="$(to_abs_path "$OUT_DIR")"

# 检查依赖
#command -v bsub >/dev/null 2>&1 || { echo "bsub not found in PATH"; exit 1; }
[[ -x "$EXEC" ]] || { echo "Executable not found or not executable: $EXEC"; exit 1; }
[[ -f "$CONFIG" ]] || { echo "Config file not found: $CONFIG"; exit 1; }
# 自动编译（可选）
if ((BUILD)); then
  if ((DRY_RUN)); then
    echo "[DRY-RUN] make -C \"$WORKSPACE_ROOT/GenLDPC\" && make -C \"$WORKSPACE_ROOT/gen4_ldpc_sim\""
  else
    make -C "$WORKSPACE_ROOT/GenLDPC"
    make -C "$WORKSPACE_ROOT/gen4_ldpc_sim"
  fi
fi
# 若提供 --gen-n/--gen-k/--gen-count，则先生成矩阵（用作后续扫描来源）
if [[ -n "$GEN_N" || -n "$GEN_K" || -n "$GEN_COUNT" ]]; then
  if [[ -z "$GEN_N" || -z "$GEN_K" || -z "$GEN_COUNT" ]]; then
    echo "生成矩阵需同时指定 --gen-n --gen-k --gen-count"; exit 1
  fi

  mk=$((GEN_N - GEN_K))
  # 设定默认矩阵目录为 GenLDPC 输出（若未显式设置 --matrix-dir）
  if [[ $MATRIX_DIR_SET -eq 0 ]]; then
    MATRIX_DIR="$(to_abs_path "GenLDPC/output/${mk}x${GEN_N}/matrix")"
    MATRIX_DIR_SET=1
  else
    MATRIX_DIR="$(to_abs_path "$MATRIX_DIR")"
  fi

  # 源目录（用于重命名前后差集）
  GEN_SRC_H_DIR="$(to_abs_path "GenLDPC/output/${mk}x${GEN_N}/matrix")"
  GEN_SRC_MASK_DIR="$(to_abs_path "GenLDPC/output/${mk}x${GEN_N}/mask_matrix")"
  mkdir -p "$GEN_SRC_H_DIR" "$GEN_SRC_MASK_DIR"
  mapfile -t FILES_BEFORE_NEWFMT < <(find "$GEN_SRC_H_DIR" -maxdepth 1 -type f -name "*_QC_H_*.txt" -print 2>/dev/null | sort)

  if [[ "$GEN_MODE" == "lsf" ]]; then
    GEN_QUEUE_USE="${GEN_QUEUE:-$QUEUE}"
    GEN_JOB_NAME="${JOB_PREFIX}_gen_${GEN_N}x${GEN_K}"
    GEN_LOG_DIR="${OUT_DIR}/logs_gen"
    mkdir -p "$GEN_LOG_DIR"
    # 导出 GENLDPC_OUT_DIR，并传递为第4参数以保持最大兼容
    gen_inner_cmd="export GENLDPC_OUT_DIR=\"$GEN_OUT_BASE\"; \"$GEN_BIN\" $GEN_N $GEN_K \$LSB_JOBINDEX \"$GEN_OUT_BASE\""
    bsub_gen=(bsub -q "$GEN_QUEUE_USE" -J "${GEN_JOB_NAME}[1-${GEN_COUNT}]" \
                   -o "${GEN_LOG_DIR}/gen_%I.out" -e "${GEN_LOG_DIR}/gen_%I.err" -- bash -lc "$gen_inner_cmd")
    if ((DRY_RUN)); then
      printf "[DRY-RUN] "; printf "%q " "${bsub_gen[@]}"; echo
    else
      echo "Submitting generator job array: ${GEN_JOB_NAME}[1-${GEN_COUNT}]"
      "${bsub_gen[@]}"
      if command -v bwait >/dev/null 2>&1; then
        echo "Waiting for generator jobs to finish (bwait)..."
        bwait -w "ended(${GEN_JOB_NAME})"
      else
        echo "Waiting for generator jobs to finish (polling bjobs)..."
        while true; do
          if command -v bjobs >/dev/null 2>&1; then
            bjobs_output=$(bjobs -J "${GEN_JOB_NAME}" 2>/dev/null || true)
            rem=$(echo "$bjobs_output" | awk 'NR>1{c++} END{print c+0}')
            if [[ -z "$bjobs_output" || -z "$rem" || "$rem" -eq 0 ]]; then
              break
            fi
          else
            echo "[WARN] bjobs not found; skip polling." >&2
            break
          fi
          sleep 5
        done
      fi
    fi
  else
    for ((ph=1; ph<=GEN_COUNT; ph++)); do
      if ((DRY_RUN)); then
        echo "[DRY-RUN] GENLDPC_OUT_DIR=\"$GEN_OUT_BASE\" \"$GEN_BIN\" $GEN_N $GEN_K $ph \"$GEN_OUT_BASE\""
      else
        GENLDPC_OUT_DIR="$GEN_OUT_BASE" "$GEN_BIN" "$GEN_N" "$GEN_K" "$ph" "$GEN_OUT_BASE"
      fi
    done
  fi

  # 重命名（沿用 test_ldpc_sim.sh 逻辑）
  if ((RENAME_ENABLE)); then
    rename_script="$WORKSPACE_ROOT/scripts/rename_ldpc_matrices.sh"
    if ((DRY_RUN)); then
      echo "[DRY-RUN] bash \"$rename_script\" $mk $GEN_N --apply"
    else
      if bash "$rename_script" "$mk" "$GEN_N" --apply; then
        echo "矩阵重命名完成"
      else
        echo "矩阵重命名失败" >&2
        exit 1
      fi
    fi
  fi

  # 统计新增矩阵（重命名后格式化）
  mapfile -t FILES_AFTER_NEWFMT < <(find "$GEN_SRC_H_DIR" -maxdepth 1 -type f -name "*_QC_H_*.txt" -print 2>/dev/null | sort)
  declare -A seen_before=()
  for f in "${FILES_BEFORE_NEWFMT[@]}"; do seen_before["$f"]=1; done
  NEW_MATRIX_FILES=()
  for f in "${FILES_AFTER_NEWFMT[@]}"; do
    [[ -z "${seen_before[$f]+x}" ]] && NEW_MATRIX_FILES+=("$f")
  done

  NEW_MATRIX_IDS=()
  for f in "${NEW_MATRIX_FILES[@]}"; do
    base="$(basename "$f")"
    if [[ "$base" =~ _([0-9]+)\.txt$ ]]; then
      NEW_MATRIX_IDS+=("${BASH_REMATCH[1]}")
    fi
  done

  if [[ ${#NEW_MATRIX_IDS[@]} -eq 0 ]]; then
    mapfile -t FILES_AFTER_OLD < <(find "$GEN_SRC_H_DIR" -maxdepth 1 -type f -name "*_QC_H_*_*.txt" -print 2>/dev/null | sort)
    for f in "${FILES_AFTER_OLD[@]}"; do
      base="$(basename "$f")"
      if [[ "$base" =~ _([0-9]+)_([0-9]+)\.txt$ ]]; then
        NEW_MATRIX_IDS+=("${BASH_REMATCH[1]}_${BASH_REMATCH[2]}")
      fi
    done
  fi

  if [[ ${#NEW_MATRIX_IDS[@]} -gt 0 ]]; then
    GENERATED=1
    echo "新增矩阵: ${#NEW_MATRIX_IDS[@]} 个 -> ${NEW_MATRIX_IDS[*]}"
  else
    echo "未检测到新增矩阵文件"
  fi
fi

[[ -d "$MATRIX_DIR" ]] || { echo "Matrix dir not found: $MATRIX_DIR"; exit 1; }

# 生成 SNR 列表
declare -a SNR_ARR=()
if [[ -n "$SNR_LIST" ]]; then
  for s in $SNR_LIST; do SNR_ARR+=("$s"); done
elif [[ -n "$SNR_SEQ" ]]; then
  IFS=':' read -r A B C <<< "$SNR_SEQ"
  mapfile -t SNR_ARR < <(awk -v a="$A" -v b="$B" -v c="$C" 'BEGIN{for(x=a; x<=c+1e-9; x+=b) printf("%.6f\n", x)}')
else
  echo "必须提供 --snr 或 --snr-seq"; exit 1
fi

# 准备矩阵编号
declare -a MATRIX_IDS=()
if ((GENERATED)) && [[ ${#NEW_MATRIX_IDS[@]} -gt 0 ]]; then
  MATRIX_IDS=("${NEW_MATRIX_IDS[@]}")
else
  declare -A SEEN=()
  while IFS= read -r -d '' f; do
    base="$(basename -- "$f")"
    id=""
    if [[ "$base" =~ _([0-9]+)(\.[^.]+)?$ ]]; then
      id="${BASH_REMATCH[1]}"
    elif [[ "$base" =~ \.([0-9]+)$ ]]; then
      id="${BASH_REMATCH[1]}"
    fi
    if [[ -n "$id" && -z "${SEEN[$id]+x}" ]]; then
      MATRIX_IDS+=("$id"); SEEN["$id"]=1
    fi
  done < <(find "$MATRIX_DIR" -maxdepth 1 -type f -name '*.txt' -print0 | sort -z)

  if [[ ${#MATRIX_IDS[@]} -eq 0 ]]; then
    count=$(find "$MATRIX_DIR" -maxdepth 1 -type f -name '*.txt' | wc -l | tr -d ' ')
    if [[ "$count" -eq 0 ]]; then
      echo "矩阵目录为空：$MATRIX_DIR"; exit 1
    fi
    for ((i=1;i<=count;i++)); do MATRIX_IDS+=("$i"); done
  fi
  mapfile -t MATRIX_IDS < <(printf "%s\n" "${MATRIX_IDS[@]}" | sort -n)
fi

MATRIX_END="${MATRIX_IDS[-1]}"

echo "Detected matrix IDs: ${MATRIX_IDS[*]}"
echo "MATRIX_END = $MATRIX_END"
echo "SNRs: ${SNR_ARR[*]}"
echo "Queue: $QUEUE"
echo "Exec:  $EXEC"
echo "Config:$CONFIG"
echo "Matrix:$MATRIX_DIR"
echo "Out:   $OUT_DIR"
echo ""

mkdir -p "$OUT_DIR"

# 生成 bsub 命令并提交
for id in "${MATRIX_IDS[@]}"; do
  for snr in "${SNR_ARR[@]}"; do
    mat_dir="${OUT_DIR}/matrix${id}"
    mkdir -p "$mat_dir"
    log_o="${mat_dir}/snr${snr}.out"
    log_e="${mat_dir}/snr${snr}.err"

    # 绝对路径的 config，避免 -cwd 切换影响
    CONFIG_ABS="$CONFIG"
    [[ "$CONFIG_ABS" = /* ]] || CONFIG_ABS="$(cd "$(dirname "$CONFIG")" && pwd)/$(basename "$CONFIG")"

    # 命令：LDPC 模式 + 指定 config + AWGN + SNR + matrix_id
    # 在命令内部 cd 到 exec 所在目录，以满足相对路径的 ./matrix/ 访问
    exec_dir="$(cd "$(dirname "$EXEC")" && pwd)"
    exec_bin="$(basename "$EXEC")"
    # 导出矩阵目录并作为第6参数传递，便于 ssd_fc 定位矩阵
    if [[ -n "${SIM_MATRIX_DIR:-}" ]]; then
      inner_cmd="cd \"$exec_dir\" && export LDPC_MATRIX_DIR=\"$SIM_MATRIX_DIR\"; ./\"$exec_bin\" LDPC \"$CONFIG_ABS\" AWGN $snr $id \"$SIM_MATRIX_DIR\""
    else
      # 若启用重命名，则默认用重命名后的目录；否则不传第6参，让程序用默认
      if [[ $MATRIX_DIR_SET -eq 1 ]]; then
        inner_cmd="cd \"$exec_dir\" && export LDPC_MATRIX_DIR=\"$MATRIX_DIR\"; ./\"$exec_bin\" LDPC \"$CONFIG_ABS\" AWGN $snr $id \"$MATRIX_DIR\""
      else
        inner_cmd="cd \"$exec_dir\" && ./\"$exec_bin\" LDPC \"$CONFIG_ABS\" AWGN $snr $id"
      fi
    fi

    if [[ -n "$CWD" ]]; then
      BSUB_CWD=(-cwd "$CWD")
    else
      BSUB_CWD=()
    fi

    bsub_cmd=(bsub -q "$QUEUE" -J "${JOB_PREFIX}_M${id}_S${snr}" -o "$log_o" -e "$log_e" "${BSUB_CWD[@]}" -- bash -lc "$inner_cmd")

    if ((DRY_RUN)); then
      printf "[DRY-RUN] "; printf "%q " "${bsub_cmd[@]}"; echo
    else
      printf "Submitting %s (SNR=%s) ...\n" "${JOB_PREFIX}_M${id}" "$snr"
      "${bsub_cmd[@]}"
    fi
  done
done

echo "Done. Logs under: $OUT_DIR"
