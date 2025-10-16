#!/usr/bin/env bash
# 独立测试：矩阵重命名 → 仿真
# 目的：验证完整流程（从现有矩阵文件开始）

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKSPACE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "$WORKSPACE_ROOT"

# ==================== 颜色输出 ====================
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

log_info() { echo -e "${BLUE}[INFO]${NC} $*"; }
log_success() { echo -e "${GREEN}[✓]${NC} $*"; }
log_warn() { echo -e "${YELLOW}[⚠]${NC} $*"; }
log_error() { echo -e "${RED}[✗]${NC} $*"; }
log_step() { echo -e "\n${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"; echo -e "${GREEN}$*${NC}"; echo -e "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n"; }

# ==================== 配置 ====================
M=20
N=149
SRC_DIR="GenLDPC/${M}x${N}"        # 生成器默认输出位置
DEST_DIR="GenLDPC/output/${M}x${N}"  # 仿真器期望的位置
TEST_OUTPUT_DIR="runs/rename_test"

SNR=5.0
CONFIG_FILE="gen4_ldpc_sim/config/sdec_dq.cnfg"
SIM_BIN="gen4_ldpc_sim/ssd_fc_dq"
AUTO_CONFIRM="${AUTO_CONFIRM:-0}"
RESET_DEST="${RESET_DEST:-1}"

# ==================== 步骤 0: 环境检查 ====================
log_step "步骤 0: 环境检查"

if [[ ! -d "$SRC_DIR" ]]; then
    log_error "源目录不存在: $SRC_DIR"
    exit 1
fi

if [[ ! -x "$SIM_BIN" ]]; then
    log_error "仿真器不存在或不可执行: $SIM_BIN"
    log_info "尝试编译..."
    make -C gen4_ldpc_sim ssd_fc_dq || { log_error "编译失败"; exit 1; }
fi

if [[ ! -f "$CONFIG_FILE" ]]; then
    log_error "配置文件不存在: $CONFIG_FILE"
    exit 1
fi

log_success "环境检查通过"

# ==================== 步骤 1: 复制文件到目标目录 ====================
log_step "步骤 1: 复制矩阵文件到仿真器期望目录"

log_info "源目录: $SRC_DIR"
log_info "目标目录: $DEST_DIR"

# 创建目标目录
mkdir -p "$DEST_DIR/matrix"
mkdir -p "$DEST_DIR/mask_matrix"
mkdir -p "$DEST_DIR/cycle_record"

if [[ "$RESET_DEST" == "1" ]]; then
    log_info "清理目标目录旧文件..."
    find "$DEST_DIR/matrix" -type f -name "*.txt" -delete 2>/dev/null || true
    find "$DEST_DIR/mask_matrix" -type f -name "*.txt" -delete 2>/dev/null || true
    find "$DEST_DIR/cycle_record" -type f -name "*.txt" -delete 2>/dev/null || true
fi

# 显示源文件
log_info "源文件列表:"
echo "━━━ H 矩阵 ━━━"
find "$SRC_DIR/matrix" -type f -name "*.txt" 2>/dev/null | sort | while read f; do
    echo "  $(basename "$f")"
done

echo ""
echo "━━━ Mask 矩阵 ━━━"
find "$SRC_DIR/mask_matrix" -type f -name "*.txt" 2>/dev/null | sort | while read f; do
    echo "  $(basename "$f")"
done

echo ""
echo "━━━ Cycle 记录 ━━━"
find "$SRC_DIR/cycle_record" -type f -name "*.txt" 2>/dev/null | sort | while read f; do
    echo "  $(basename "$f")"
done

# 复制文件
log_info "开始复制..."
cp -v "$SRC_DIR/matrix/"* "$DEST_DIR/matrix/" 2>/dev/null || log_warn "无 H 矩阵文件"
cp -v "$SRC_DIR/mask_matrix/"* "$DEST_DIR/mask_matrix/" 2>/dev/null || log_warn "无 Mask 文件"
cp -v "$SRC_DIR/cycle_record/"* "$DEST_DIR/cycle_record/" 2>/dev/null || log_warn "无 Cycle 文件"

log_success "文件复制完成"

# ==================== 步骤 2: 检查文件格式 ====================
log_step "步骤 2: 检查文件命名格式"

log_info "当前文件格式:"
find "$DEST_DIR/matrix" -type f -name "*.txt" | head -3 | while read f; do
    basename_f="$(basename "$f")"
    echo "  $basename_f"
    
    # 分析格式
    if [[ "$basename_f" =~ _([0-9]+)_([0-9]+)\.txt$ ]]; then
        file_num="${BASH_REMATCH[1]}"
        phase="${BASH_REMATCH[2]}"
        echo "    └─ 格式: 旧格式 (file_num=$file_num, phase=$phase)"
        echo "    └─ 需要重命名为: ${basename_f%_*_*.txt}_${phase}.txt"
    elif [[ "$basename_f" =~ _([0-9]+)\.txt$ ]]; then
        id="${BASH_REMATCH[1]}"
        echo "    └─ 格式: 新格式 (id=$id)"
        echo "    └─ 无需重命名"
    else
        echo "    └─ 格式: 未知"
    fi
    echo ""
done

# ==================== 步骤 3: 执行重命名（预览） ====================
log_step "步骤 3: 重命名预览"

log_info "执行重命名脚本（预览模式）..."
set +e
bash scripts/rename_ldpc_matrices.sh "$M" "$N"
preview_status=$?
set -e
if [[ $preview_status -eq 0 ]]; then
    log_success "预览完成"
else
    log_error "预览失败 (退出码: $preview_status)"
    exit 1
fi

# 等待用户确认
echo ""
log_warn "请检查上述预览结果"
if [[ "$AUTO_CONFIRM" == "1" ]]; then
    log_info "自动确认重命名"
else
    read -p "继续执行重命名? [y/N] " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        log_error "用户取消"
        exit 1
    fi
fi

# ==================== 步骤 4: 执行重命名（应用） ====================
log_step "步骤 4: 执行重命名"

log_info "应用重命名..."
set +e
bash scripts/rename_ldpc_matrices.sh "$M" "$N" --apply
apply_status=$?
set -e
if [[ $apply_status -eq 0 ]]; then
    log_success "重命名完成"
else
    log_error "重命名失败 (退出码: $apply_status)"
    exit 1
fi

# 验证结果
log_info "重命名后的文件:"
find "$DEST_DIR/matrix" -type f -name "*.txt" | sort | while read f; do
    echo "  $(basename "$f")"
done

# ==================== 步骤 5: 提取矩阵ID ====================
log_step "步骤 5: 提取可用的矩阵ID"

declare -a MATRIX_IDS=()
while IFS= read -r f; do
    basename_f="$(basename "$f")"
    if [[ "$basename_f" =~ _([0-9]+)\.txt$ ]]; then
        id="${BASH_REMATCH[1]}"
        mask_pattern="$DEST_DIR/mask_matrix/*_mask_${id}.txt"
        cycle_pattern="$DEST_DIR/cycle_record/test_*_${id}.txt"
        if compgen -G "$mask_pattern" > /dev/null; then
            MATRIX_IDS+=("$id")
        else
            log_warn "跳过矩阵 $basename_f: 缺少对应的 mask 文件 (ID=$id)"
        fi
    fi
done < <(find "$DEST_DIR/matrix" -maxdepth 1 -type f -name "*_QC_H_*.txt" | sort -V)

if [[ ${#MATRIX_IDS[@]} -eq 0 ]]; then
    log_error "未找到任何可用的矩阵文件"
    log_info "目录内容:"
    ls -lh "$DEST_DIR/matrix/"
    exit 1
fi

log_info "检测到 ${#MATRIX_IDS[@]} 个矩阵: ${MATRIX_IDS[*]}"

# ==================== 步骤 6: 运行仿真 ====================
log_step "步骤 6: 运行仿真测试"

mkdir -p "$TEST_OUTPUT_DIR"

CONFIG_ABS="$WORKSPACE_ROOT/$CONFIG_FILE"
MATRIX_DIR_ARG="$WORKSPACE_ROOT/${DEST_DIR%/${M}x${N}}"

log_info "配置:"
log_info "  仿真器: $SIM_BIN"
log_info "  配置文件: $CONFIG_ABS"
log_info "  矩阵目录: $MATRIX_DIR_ARG"
log_info "  信道/SNR: AWGN / $SNR dB"
log_info "  矩阵ID: ${MATRIX_IDS[*]}"

total=${#MATRIX_IDS[@]}
success=0
failed=0

for mat_id in "${MATRIX_IDS[@]}"; do
    log_info "[$((success+failed+1))/$total] 测试矩阵 ID=$mat_id..."
    
    log_file="$TEST_OUTPUT_DIR/matrix_${mat_id}_snr${SNR}.log"
    
    # 构造命令
    cmd="$WORKSPACE_ROOT/$SIM_BIN LDPC \"$CONFIG_ABS\" AWGN $SNR $mat_id \"$MATRIX_DIR_ARG\""
    
    log_info "执行: $cmd"
    
    # 运行仿真
    if eval "$cmd" > "$log_file" 2>&1; then
        log_success "矩阵 $mat_id 仿真成功"
        ((success++))
        
        # 显示关键结果
        if grep -q "BER" "$log_file"; then
            echo "  结果摘要:"
            grep -E "(BER|FER|Average)" "$log_file" | head -5 | sed 's/^/    /'
        fi
    else
        log_error "矩阵 $mat_id 仿真失败"
        ((failed++))
        
        # 显示错误信息
        echo "  错误日志 (最后20行):"
        tail -20 "$log_file" | sed 's/^/    /'
    fi
    
    echo ""
done

# ==================== 步骤 7: 总结 ====================
log_step "测试完成"

log_info "━━━ 测试统计 ━━━"
log_info "  总计: $total"
log_success "  成功: $success"
if [[ $failed -gt 0 ]]; then
    log_error "  失败: $failed"
else
    log_info "  失败: $failed"
fi

log_info ""
log_info "━━━ 输出位置 ━━━"
log_info "  矩阵目录: $DEST_DIR"
log_info "  日志目录: $TEST_OUTPUT_DIR"

if [[ $failed -eq 0 ]]; then
    log_success "🎉 所有测试通过！"
    exit 0
else
    log_error "❌ 有 $failed 个测试失败"
    exit 1
fi