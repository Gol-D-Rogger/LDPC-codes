#!/usr/bin/env bash
# LDPC 仿真测试脚本（无bsub版本）
# 功能：生成矩阵 → 运行仿真（串行）
# 
# 用法示例:
#   bash scripts/test_ldpc_sim.sh \
#     --gen-n 149 --gen-k 129 --gen-count 1 \
#     --snr "5.0" \
#     --out-dir runs/test_149x129

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKSPACE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ==================== 默认配置 ====================
# 路径配置
GEN_BIN="GenLDPC/ldpc_gen"
SIM_BIN="gen4_ldpc_sim/ssd_fc_dq"
CONFIG_FILE="gen4_ldpc_sim/config/sdec_dq.cnfg"

# 生成矩阵参数
GEN_N="149"           # 例如 149
GEN_K="129"           # 例如 129
GEN_COUNT="1"      # 单次生成的 target 数量（filenum，供后续扩展）
GEN_PHASE="1"      # 矩阵 phase 标识，用于区分不同子任务
GEN_OUT_BASE="GenLDPC/output"    # 输出根目录（可选，默认用环境变量）

# 仿真参数
SNR_LIST="4.9 5.0"        # SNR列表，例如 "5.0 5.5 6.0"
SNR_SEQ=""         # SNR区间，例如 "5.0:0.5:6.0"
CHANNEL="AWGN"     # 信道类型
OUT_DIR="runs/test_sim"  # 输出目录

# 控制选项
BUILD=1            # 是否自动编译（默认开启）
DRY_RUN=0          # 仅打印命令不执行
VERBOSE=1          # 详细输出

# ==================== 函数定义 ====================
log_info() {
    echo "[INFO] $*"
}

log_warn() {
    echo "[WARN] $*"
}

log_error() {
    echo "[ERROR] $*" >&2
}

usage() {
    cat << EOF
用法: $0 [选项]

生成矩阵:
  --gen-n N              矩阵参数 N（必需）
  --gen-k K              矩阵参数 K（必需）
  --gen-count NUM        单次生成的 target 数量（默认: ${GEN_COUNT}）
  --phase ID             phase 标识（默认: ${GEN_PHASE}）
  --gen-out-base DIR     生成输出根目录（可选）

仿真参数:
  --snr "A B C"          SNR列表（空格分隔）
  --snr-seq A:B:C        SNR区间（起:步:止）
  --channel TYPE         信道类型（默认: ${CHANNEL}）
  --out-dir DIR          输出目录（默认: ${OUT_DIR}）

路径配置:
  --gen-bin PATH         生成器路径（默认: ${GEN_BIN}）
  --sim-bin PATH         仿真器路径（默认: ${SIM_BIN}）
  --config PATH          配置文件（默认: ${CONFIG_FILE}）

控制选项:
  --no-build             不自动编译
  --dry-run              仅打印命令
  --quiet                减少输出
  --help                 显示帮助

示例:
  # 生成1个矩阵，SNR=5.0
  $0 --gen-n 149 --gen-k 129 --phase 1 --snr "5.0"
  
  # 指定 phase=3，SNR=5.0~6.0，步长0.5
  $0 --gen-n 149 --gen-k 129 --phase 3 --snr-seq "5.0:0.5:6.0"
EOF
}

# ==================== 参数解析 ====================
while [[ $# -gt 0 ]]; do
    case "$1" in
        --gen-n) GEN_N="$2"; shift 2;;
        --gen-k) GEN_K="$2"; shift 2;;
        --gen-count) GEN_COUNT="$2"; shift 2;;
        --phase|--phase-id) GEN_PHASE="$2"; shift 2;;
        --gen-out-base) GEN_OUT_BASE="$2"; shift 2;;
        --gen-bin) GEN_BIN="$2"; shift 2;;
        --sim-bin) SIM_BIN="$2"; shift 2;;
        --config) CONFIG_FILE="$2"; shift 2;;
        --snr) SNR_LIST="$2"; shift 2;;
        --snr-seq) SNR_SEQ="$2"; shift 2;;
        --channel) CHANNEL="$2"; shift 2;;
        --out-dir) OUT_DIR="$2"; shift 2;;
        --no-build) BUILD=0; shift;;
        --dry-run) DRY_RUN=1; shift;;
        --quiet) VERBOSE=0; shift;;
        --help|-h) usage; exit 0;;
        *) log_error "未知选项: $1"; usage; exit 1;;
    esac
done

# ==================== 参数验证 ====================
if [[ -z "$GEN_N" || -z "$GEN_K" ]]; then
    log_error "必须指定 --gen-n 和 --gen-k"
    usage
    exit 1
fi

if [[ -z "$SNR_LIST" && -z "$SNR_SEQ" ]]; then
    log_error "必须指定 --snr 或 --snr-seq"
    usage
    exit 1
fi

# 计算 m = n - k
m=$((GEN_N - GEN_K))
if [[ $m -le 0 ]]; then
    log_error "无效的参数: N=$GEN_N, K=$GEN_K (N必须大于K)"
    exit 1
fi

# ==================== 生成 SNR 列表 ====================
declare -a SNR_ARR=()
if [[ -n "$SNR_LIST" ]]; then
    for s in $SNR_LIST; do 
        SNR_ARR+=("$s")
    done
elif [[ -n "$SNR_SEQ" ]]; then
    IFS=':' read -r START STEP END <<< "$SNR_SEQ"
    mapfile -t SNR_ARR < <(awk -v a="$START" -v b="$STEP" -v c="$END" \
        'BEGIN{for(x=a; x<=c+1e-9; x+=b) printf("%.6f\n", x)}')
fi

# ==================== 路径设置 ====================
cd "$WORKSPACE_ROOT"

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

GEN_BIN_ABS="$(to_abs_path "$GEN_BIN")"
SIM_BIN_ABS="$(to_abs_path "$SIM_BIN")"
CONFIG_ABS="$(to_abs_path "$CONFIG_FILE")"
GEN_OUT_BASE_ABS="$(to_abs_path "$GEN_OUT_BASE")"

GEN_BIN_DISPLAY="${GEN_BIN_ABS#$WORKSPACE_ROOT/}"
GEN_BIN_DISPLAY="${GEN_BIN_DISPLAY:-$GEN_BIN_ABS}"
SIM_BIN_DISPLAY="${SIM_BIN_ABS#$WORKSPACE_ROOT/}"
SIM_BIN_DISPLAY="${SIM_BIN_DISPLAY:-$SIM_BIN_ABS}"
CONFIG_DISPLAY="${CONFIG_ABS#$WORKSPACE_ROOT/}"
CONFIG_DISPLAY="${CONFIG_DISPLAY:-$CONFIG_ABS}"

# 构造矩阵目录: <GEN_OUT_BASE>/<M>x<N>
MATRIX_BASE_DIR="${GEN_OUT_BASE_ABS}/${m}x${GEN_N}"
MATRIX_H_DIR="${MATRIX_BASE_DIR}/matrix"
MATRIX_MASK_DIR="${MATRIX_BASE_DIR}/mask_matrix"

MATRIX_BASE_DISPLAY="${MATRIX_BASE_DIR#$WORKSPACE_ROOT/}"
MATRIX_BASE_DISPLAY="${MATRIX_BASE_DISPLAY:-$MATRIX_BASE_DIR}"

# ==================== 显示配置 ====================
log_info "=========================================="
log_info "LDPC 仿真测试配置"
log_info "=========================================="
log_info "工作目录: $WORKSPACE_ROOT"
log_info "矩阵参数: M=$m, N=$GEN_N, K=$GEN_K"
log_info "生成数量(filenum): $GEN_COUNT"
log_info "Phase 标识: $GEN_PHASE"
log_info "矩阵目录: $MATRIX_BASE_DISPLAY"
log_info "SNR列表: ${SNR_ARR[*]}"
log_info "信道类型: $CHANNEL"
log_info "仿真器: $SIM_BIN_DISPLAY"
log_info "配置文件: $CONFIG_DISPLAY"
log_info "输出目录: $OUT_DIR"
log_info "=========================================="

# ==================== 编译 ====================
if ((BUILD)); then
    log_info "开始编译..."
    
    if ((DRY_RUN)); then
        echo "[DRY-RUN] make -C GenLDPC"
        echo "[DRY-RUN] make -C gen4_ldpc_sim ssd_fc_dq"
    else
        log_info "编译 GenLDPC..."
        make -C GenLDPC || { log_error "GenLDPC 编译失败"; exit 1; }
        
        log_info "编译 gen4_ldpc_sim (DQ版本)..."
        make -C gen4_ldpc_sim ssd_fc_dq || { log_error "gen4_ldpc_sim 编译失败"; exit 1; }
    fi
    
    log_info "编译完成"
fi

# 验证可执行文件
if [[ ! -x "$GEN_BIN_ABS" ]]; then
    log_error "生成器不存在或不可执行: $GEN_BIN_ABS"
    exit 1
fi

if [[ ! -x "$SIM_BIN_ABS" ]]; then
    log_error "仿真器不存在或不可执行: $SIM_BIN_ABS"
    exit 1
fi

if [[ ! -f "$CONFIG_ABS" ]]; then
    log_error "配置文件不存在: $CONFIG_ABS"
    exit 1
fi

# ==================== 生成矩阵 ====================
log_info "=========================================="
log_info "开始生成矩阵..."
log_info "=========================================="

# 创建输出目录
mkdir -p "$MATRIX_H_DIR" "$MATRIX_MASK_DIR"

# 记录生成前的文件（用于后续提取新文件）
if [[ ! $DRY_RUN -eq 1 ]]; then
mapfile -t FILES_BEFORE_NEWFMT < <(find "$MATRIX_H_DIR" -maxdepth 1 -type f -name "*_QC_H_*.txt" 2>/dev/null | sort)
fi

# 生成矩阵（单次调用，对应指定 phase）
log_info "生成矩阵（phase = $GEN_PHASE）..."

if [[ -n "$GEN_OUT_BASE" ]]; then
    gen_cmd="GENLDPC_OUT_DIR=\"$GEN_OUT_BASE_ABS\" \"$GEN_BIN_ABS\" $GEN_N $GEN_K $GEN_PHASE \"$GEN_OUT_BASE_ABS\""
else
    gen_cmd="\"$GEN_BIN_ABS\" $GEN_N $GEN_K $GEN_PHASE"
fi

if ((DRY_RUN)); then
    echo "[DRY-RUN] $gen_cmd"
else
    if ((VERBOSE)); then
        eval "$gen_cmd"
    else
        eval "$gen_cmd" > /dev/null 2>&1
    fi
    
    if [[ $? -eq 0 ]]; then
        log_info "✓ Phase $GEN_PHASE 生成成功"
    else
        log_error "✗ Phase $GEN_PHASE 生成失败"
        exit 1
    fi
fi

# 对新生成的矩阵执行标准重命名（仅在非 DRY_RUN）
if ((DRY_RUN)); then
    log_info "[DRY-RUN] bash \"$WORKSPACE_ROOT/scripts/rename_ldpc_matrices.sh\" $m $GEN_N --apply"
else
    log_info "执行矩阵重命名（标准编号）..."
    if bash "$WORKSPACE_ROOT/scripts/rename_ldpc_matrices.sh" "$m" "$GEN_N" --apply; then
        log_info "矩阵重命名完成"
    else
        log_error "矩阵重命名失败"
        exit 1
    fi
fi

# ==================== 提取新生成的矩阵ID ====================
declare -a MATRIX_IDS=()

if [[ ! $DRY_RUN -eq 1 ]]; then
    mapfile -t FILES_AFTER_NEWFMT < <(find "$MATRIX_H_DIR" -maxdepth 1 -type f -name "*_QC_H_*.txt" 2>/dev/null | sort)
    
    # 找出新文件
    declare -A seen_before=()
    for f in "${FILES_BEFORE_NEWFMT[@]}"; do
        seen_before["$f"]=1
    done
    
    for f in "${FILES_AFTER_NEWFMT[@]}"; do
        if [[ -z "${seen_before[$f]+x}" ]]; then
            # 新文件，提取ID
            basename_f="$(basename "$f")"
            if [[ "$basename_f" =~ _([0-9]+)\.txt$ ]]; then
                id="${BASH_REMATCH[1]}"
                MATRIX_IDS+=("$id")
                log_info "检测到新矩阵: $basename_f (ID: $id)"
            fi
        fi
    done
    
    if [[ ${#MATRIX_IDS[@]} -eq 0 ]]; then
        log_warn "未检测到新增的新格式矩阵，尝试扫描旧格式..."
        mapfile -t FILES_AFTER_OLD < <(find "$MATRIX_H_DIR" -maxdepth 1 -type f -name "*_QC_H_*_*.txt" 2>/dev/null | sort)
        for f in "${FILES_AFTER_OLD[@]}"; do
            basename_f="$(basename "$f")"
            if [[ "$basename_f" =~ _([0-9]+)_([0-9]+)\.txt$ ]]; then
                id1="${BASH_REMATCH[1]}"
                id2="${BASH_REMATCH[2]}"
                combined_id="${id1}_${id2}"
                MATRIX_IDS+=("$combined_id")
                log_warn "检测到未重命名的矩阵: $basename_f (旧格式ID: $combined_id)"
            fi
        done
    fi
    
    if [[ ${#MATRIX_IDS[@]} -eq 0 ]]; then
        log_error "未检测到新生成的矩阵文件"
        log_info "尝试列出 $MATRIX_H_DIR 目录:"
        ls -lh "$MATRIX_H_DIR"
        exit 1
    fi
else
    # DRY_RUN 模式，模拟ID
    MATRIX_IDS+=("phase_${GEN_PHASE}")
fi

log_info "检测到 ${#MATRIX_IDS[@]} 个矩阵: ${MATRIX_IDS[*]}"

# ==================== 运行仿真 ====================
log_info "=========================================="
log_info "开始运行仿真（串行）..."
log_info "=========================================="

mkdir -p "$OUT_DIR"

total_jobs=$((${#MATRIX_IDS[@]} * ${#SNR_ARR[@]}))
current_job=0

cd "$WORKSPACE_ROOT" || { log_error "无法切换到工作目录: $WORKSPACE_ROOT"; exit 1; }


for mat_id in "${MATRIX_IDS[@]}"; do
    for snr in "${SNR_ARR[@]}"; do
        current_job=$((current_job + 1))
        
        log_info "[$current_job/$total_jobs] 运行: Matrix=$mat_id, SNR=$snr"
        
        # 输出文件
        log_dir="$OUT_DIR/matrix_${mat_id}"
        mkdir -p "$log_dir"
        log_file="${log_dir}/snr_${snr}.log"
        
        # 构造仿真命令
        # 参数: LDPC <config> <channel> <snr> <matrix_id> <matrix_dir>
        # matrix_dir 应该指向 GenLDPC/output (不包含 MxN 子目录)
        matrix_dir_arg="$GEN_OUT_BASE_ABS"
        config_abs_path="$CONFIG_ABS"

        
        # 提取实际的matrix_id（第二个数字）
        if [[ "$mat_id" =~ _([0-9]+)$ ]]; then
            actual_matrix_id="${BASH_REMATCH[1]}"
        else
            actual_matrix_id="$mat_id"
        fi
        
        sim_cmd="\"$SIM_BIN_ABS\" LDPC \"$config_abs_path\" $CHANNEL $snr $actual_matrix_id \"$matrix_dir_arg\""
        
        if ((DRY_RUN)); then
            echo "[DRY-RUN] cd $WORKSPACE_ROOT && $sim_cmd 2>&1 | tee \"$log_file\""
        else
            if ((VERBOSE)); then
                log_info "工作目录: $(pwd)"
                log_info "执行: $sim_cmd"
            fi
            
            # 执行并记录
            (
                cd "$WORKSPACE_ROOT" || exit 1
                eval "$sim_cmd" 2>&1 | tee "$log_file"
            )
            exit_code=${PIPESTATUS[0]}
            
            if [[ $exit_code -eq 0 ]]; then
                log_info "✓ 完成: Matrix=$mat_id, SNR=$snr"
            else
                log_error "✗ 失败: Matrix=$mat_id, SNR=$snr (退出码: $exit_code)"
                log_error "查看日志: $log_file"
                # 【可选】显示错误上下文
                if ((VERBOSE)); then
                    log_error "错误日志末尾:"
                    tail -20 "$log_file" | sed 's/^/  /'
                fi
                exit 1
            fi
        fi
    done
done

# ==================== 完成 ====================
log_info "=========================================="
log_info "仿真完成！"
log_info "=========================================="
log_info "输出目录: $OUT_DIR"
log_info "矩阵目录: $MATRIX_BASE_DIR"
log_info ""
log_info "查看结果:"
log_info "  ls -lh $OUT_DIR/*/snr_*.log"
log_info ""

if ((DRY_RUN)); then
    log_info "这是 DRY-RUN 模式，实际未执行任何命令"
fi
