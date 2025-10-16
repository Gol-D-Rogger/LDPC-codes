#!/usr/bin/env bash
# 快速测试脚本 - 使用小参数
# 功能：生成矩阵 -> 编译DQ版本 -> 运行仿真

set -euo pipefail

# ==================== 配置参数 ====================
# 矩阵生成参数（小参数快速测试）
GEN_N=100           # 列数（确保 N-K >= 6 以匹配 g=6）
GEN_K=80           # 行数（K < N）
GEN_COUNT=1        # 生成1个矩阵
GEN_PHASE=1        # phase编号

# 仿真参数
SNR="4.9 5.0"          # 单个SNR点
MATRIX_ID=1        # 使用第1个矩阵

# 路径配置
WORKSPACE_ROOT="/workspaces/LDPC-codes"
GEN_BIN="${WORKSPACE_ROOT}/GenLDPC/ldpc_gen"
SIM_BIN="${WORKSPACE_ROOT}/gen4_ldpc_sim/ssd_fc_dq"
CONFIG="${WORKSPACE_ROOT}/runs/quick_test/sdec_quick.cnfg"  # 动态生成的配置
OUT_DIR="${WORKSPACE_ROOT}/runs/quick_test"

# 自动计算
M=$((GEN_N - GEN_K))  # M = N - K = 5
MATRIX_BASE_DIR="${WORKSPACE_ROOT}/GenLDPC/output/${M}x${GEN_N}"
MATRIX_DIR="${MATRIX_BASE_DIR}/matrix"
MASK_DIR="${MATRIX_BASE_DIR}/mask_matrix"

# ==================== 函数定义 ====================
log_info() {
    echo -e "\n\033[1;36m[INFO]\033[0m $*"
}

log_success() {
    echo -e "\033[1;32m[SUCCESS]\033[0m $*"
}

log_error() {
    echo -e "\033[1;31m[ERROR]\033[0m $*"
}

# ==================== 步骤1: 编译 ====================
log_info "步骤1: 编译 GenLDPC 和 gen4_ldpc_sim (DQ版本)"

cd "${WORKSPACE_ROOT}"

# 确保输出目录存在
mkdir -p "${OUT_DIR}"

if ! make -C GenLDPC 2>&1 | tail -5; then
    log_error "GenLDPC 编译失败"
    exit 1
fi
log_success "GenLDPC 编译完成"

if ! make -C gen4_ldpc_sim ssd_fc_dq 2>&1 | tail -5; then
    log_error "gen4_ldpc_sim (DQ) 编译失败"
    exit 1
fi
log_success "gen4_ldpc_sim (DQ) 编译完成"

# ==================== 步骤2: 生成矩阵 ====================
log_info "步骤2: 生成 LDPC 矩阵 (N=${GEN_N}, K=${GEN_K}, Phase=${GEN_PHASE})"

# 创建输出目录
mkdir -p "${MATRIX_DIR}" "${MASK_DIR}"

# 记录生成前的文件
BEFORE_H=$(ls -1 "${MATRIX_DIR}"/*_*_*.txt 2>/dev/null | wc -l)
BEFORE_M=$(ls -1 "${MASK_DIR}"/*_*_*.txt 2>/dev/null | wc -l)

log_info "生成前: H矩阵=${BEFORE_H}个, Mask矩阵=${BEFORE_M}个"

# 执行生成
cd "${WORKSPACE_ROOT}/GenLDPC"
log_info "执行: ./ldpc_gen ${GEN_N} ${GEN_K} ${GEN_PHASE}"

# 测试模式下强制保存（即使 PHI 不可逆也会写出文件），便于连通性验证
export LDPC_GEN_FORCE_SAVE=1
if ! ./ldpc_gen ${GEN_N} ${GEN_K} ${GEN_PHASE} 2>&1 | tee "${OUT_DIR}/gen_matrix.log"; then
    log_error "矩阵生成失败"
    exit 1
fi

# 记录生成后的文件
AFTER_H=$(ls -1 "${MATRIX_DIR}"/*_*_*.txt 2>/dev/null | wc -l)
AFTER_M=$(ls -1 "${MASK_DIR}"/*_*_*.txt 2>/dev/null | wc -l)

NEW_H=$((AFTER_H - BEFORE_H))
NEW_M=$((AFTER_M - BEFORE_M))

log_success "矩阵生成完成: 新增 H矩阵=${NEW_H}个, Mask矩阵=${NEW_M}个"

# 列出新生成的文件
log_info "生成的矩阵文件:"
ls -lh "${MATRIX_DIR}"/*_*_*.txt 2>/dev/null | tail -3
log_info "生成的掩码文件:"
ls -lh "${MASK_DIR}"/*_*_*.txt 2>/dev/null | tail -3

# ==================== 步骤3: 检查文件名格式 ====================
log_info "步骤3: 检查文件名格式"

# 找到最新的H矩阵文件
H_FILE=$(ls -t "${MATRIX_DIR}"/*_*_*.txt 2>/dev/null | head -1)
M_FILE=$(ls -t "${MASK_DIR}"/*_*_*.txt 2>/dev/null | head -1)

if [[ -z "$H_FILE" ]]; then
    log_error "未找到生成的H矩阵文件"
    exit 1
fi

log_info "H矩阵文件: $(basename "$H_FILE")"
log_info "Mask文件: $(basename "$M_FILE")"

H_BASENAME=$(basename "$H_FILE")
if [[ "$H_BASENAME" =~ _([0-9]+)_([0-9]+)\.txt$ ]]; then
    FILE_NUM="${BASH_REMATCH[1]}"
    PHASE_NUM="${BASH_REMATCH[2]}"
    log_info "文件编号: file_num=${FILE_NUM}, phase=${PHASE_NUM}"
else
    log_error "无法从文件名提取编号: $H_BASENAME"
    exit 1
fi

# ==================== 生成仿真配置（与小矩阵匹配） ====================
log_info "生成 DQ 配置文件: ${CONFIG}"
mkdir -p "${OUT_DIR}"
cat >"${CONFIG}" <<EOF
12                              // meta data size (B)
256                             // pad bit size (bit)
512                             // LBA size(B)
8                               // LBA # per DSP CW
${M}                              // row number of base matrix (N-K)
${GEN_N}                          // column number of base matrix (N)
256                             // circulant number
6                               // Dense matrix size or rows of E matrix in base matrix
6                               // column weight of base matrix
1000                            // maximum simulation number
10                              // maximum error number
FC_FDEC                         // LDPC decoder(FC_SKIP/FC_FDEC/FC_RDEC/FC_MIX/FC_FDEC2)
24                              // Fast decoder maximum iteration number 24
0                               // Iteration to turn on FDEC column skip feature (0:always OFF)
32                              // Retry decoder maximum iteration number
0.625                           // Retry decoder alpha
1                               // Quantization Mode
8                               // APP/Q bit width
6                               // R bit width (Min need 1-bit less)
3                               // Fraction bits
7                               // LLR bit width
3                               // LLR fraction bits
MANUAL                          // LLR mode (IDEAL/MANUAL/VENDOR0/VENDOR1)
1                               // Number of NAND soft data bits
0 0.15 -0.15 0.3 -0.3 0.5 -0.5  // VERF values
1.875 -1.875                    // Hard decoding LLR for retry decoder
1060                            // init_synd_wt threshold to skip FDEC
104C11DB7                       // CRC 32 polynomial values
104C11DB7                       // 32-bit LSFR for randomizer
EOF

# ==================== 步骤4: 运行仿真 ====================
log_info "步骤4: 运行 LDPC 仿真 (SNR=${SNR}, Matrix phase=${PHASE_NUM})"

cd "${WORKSPACE_ROOT}/gen4_ldpc_sim"

# 创建输出目录
SIM_OUT_DIR="${OUT_DIR}/matrix_${FILE_NUM}_${PHASE_NUM}"
mkdir -p "${SIM_OUT_DIR}"

# 设置矩阵目录环境变量
export LDPC_MATRIX_DIR="${WORKSPACE_ROOT}/GenLDPC/output"

# 构造仿真命令
SIM_CMD="./ssd_fc_dq LDPC ${CONFIG} AWGN ${SNR} ${PHASE_NUM} ${LDPC_MATRIX_DIR}"

log_info "仿真命令: ${SIM_CMD}"
log_info "输出日志: ${SIM_OUT_DIR}/sim.log"

# 运行仿真并记录输出
if ${SIM_CMD} 2>&1 | tee "${SIM_OUT_DIR}/sim.log"; then
    log_success "仿真运行完成"
else
    EXIT_CODE=$?
    log_error "仿真运行失败 (退出码: ${EXIT_CODE})"
    log_info "检查日志: ${SIM_OUT_DIR}/sim.log"
    exit 1
fi

# ==================== 步骤5: 结果总结 ====================
log_info "步骤5: 测试总结"

echo ""
echo "=========================================="
echo "          快速测试完成"
echo "=========================================="
echo "矩阵参数: N=${GEN_N}, K=${GEN_K}, M=${M}"
echo "矩阵文件: ${FILE_NUM}_${PHASE_NUM}"
echo "SNR:      ${SNR} dB"
echo "输出目录: ${OUT_DIR}"
echo "=========================================="
echo ""

log_info "查看完整日志:"
echo "  矩阵生成: ${OUT_DIR}/gen_matrix.log"
echo "  仿真输出: ${SIM_OUT_DIR}/sim.log"

log_success "全部测试完成！"
