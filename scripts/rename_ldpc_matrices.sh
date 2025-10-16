#!/usr/bin/env bash
# LDPC 矩阵智能重命名脚本
# 功能：将 _X_Y.txt 格式重命名为连续编号 _N.txt
# 特性：自动检测已有最大编号，新文件从下一个编号开始
#
# 用法：
#   预览：bash rename_ldpc_matrices.sh <M> <N>
#   执行：bash rename_ldpc_matrices.sh <M> <N> --apply
#
# 示例：
#   bash scripts/rename_ldpc_matrices.sh 20 149 --apply

set -uo pipefail

# 颜色定义
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

log_info() { echo -e "${BLUE}[INFO]${NC} $*"; }
log_ok() { echo -e "${GREEN}[OK]${NC} $*"; }
log_warn() { echo -e "${YELLOW}[WARN]${NC} $*"; }
log_error() { echo -e "${RED}[ERROR]${NC} $*" >&2; }

# ==================== 参数解析 ====================
if [[ $# -lt 2 ]]; then
    cat << 'EOF'
用法: ./rename_ldpc_matrices.sh <M> <N> [--apply]

参数:
  M        矩阵行数 (M = N - K)
  N        矩阵列数
  --apply  执行重命名（默认仅预览）

功能:
  - 扫描 GenLDPC/output/<M>x<N> 目录
  - 将 *_X_Y.txt 格式重命名为 *_N.txt
  - 自动检测已有的 *_N.txt 文件，从最大编号+1开始
  - 按 (file_num, phase) 排序后连续编号

示例:
  # 预览
  ./rename_ldpc_matrices.sh 20 149

  # 执行
  ./rename_ldpc_matrices.sh 20 149 --apply

  # 目录已有 _1.txt, _2.txt
  # 新文件 _1_3.txt, _1_4.txt 会被重命名为 _3.txt, _4.txt
EOF
    exit 1
fi

M="$1"
N="$2"
APPLY=0

shift 2
while [[ $# -gt 0 ]]; do
    case "$1" in
        --apply) APPLY=1; shift;;
        *) log_error "未知选项: $1"; exit 1;;
    esac
done
# ==================== 路径设置 ====================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKSPACE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
BASE_DIR="${WORKSPACE_ROOT}/GenLDPC/output/${M}x${N}"
MATRIX_DIR="${BASE_DIR}/matrix"
MASK_DIR="${BASE_DIR}/mask_matrix"
CYCLE_DIR="${BASE_DIR}/cycle_record"

log_info "=========================================="
log_info "LDPC 矩阵智能重命名"
log_info "=========================================="
log_info "矩阵尺寸: ${M}x${N}"
log_info "基础目录: ${BASE_DIR}"
log_info "模式: $(((APPLY)) && echo '执行重命名' || echo '预览模式')"
log_info "=========================================="
echo ""

if [[ ! -d "$BASE_DIR" ]]; then
    log_error "目录不存在: $BASE_DIR"
    exit 1
fi

# ==================== 智能重命名函数 ====================
rename_files_smart() {
    local dir="$1"
    local pattern_old="$2"  # 旧格式匹配
    local pattern_new="$3"  # 新格式匹配
    local desc="$4"
    
    if [[ ! -d "$dir" ]]; then
        log_warn "$desc 目录不存在: $dir"
        return
    fi
    
    # 查找已有的新格式文件（例如 *_3.txt）
    local max_existing=0
    while IFS= read -r -d '' f; do
        local basename_f=$(basename "$f")
        # 匹配新格式：...QC_H_<num>.txt
        if [[ "$basename_f" =~ _QC_H_([0-9]+)\.txt$ ]]; then
            local num="${BASH_REMATCH[1]}"
            ((num > max_existing)) && max_existing=$num
        elif [[ "$basename_f" =~ _mask_([0-9]+)\.txt$ ]]; then
            local num="${BASH_REMATCH[1]}"
            ((num > max_existing)) && max_existing=$num
        elif [[ "$basename_f" =~ test_[0-9]+x[0-9]+_([0-9]+)\.txt$ ]]; then
            local num="${BASH_REMATCH[1]}"
            ((num > max_existing)) && max_existing=$num
        fi
    done < <(find "$dir" -maxdepth 1 -type f -name "$pattern_new" -print0 2>/dev/null || true)
    
    # 查找需要重命名的旧格式文件（例如 *_1_2.txt）
    declare -a old_files=()
    declare -a sort_keys=()
    
    while IFS= read -r -d '' f; do
        local basename_f=$(basename "$f")
        
        # 匹配旧格式：...QC_H_<file_num>_<phase>.txt
        if [[ "$basename_f" =~ ^(.*)_QC_H_([0-9]+)_([0-9]+)(\.txt.*)$ ]]; then
            local prefix="${BASH_REMATCH[1]}"
            local file_num="${BASH_REMATCH[2]}"
            local phase="${BASH_REMATCH[3]}"
            local suffix="${BASH_REMATCH[4]}"
            # 排序键：file_num*10000 + phase，确保正确排序
            local sort_key=$((file_num * 10000 + phase))
            old_files+=("$f")
            sort_keys+=("$sort_key")
        elif [[ "$basename_f" =~ ^(.*)_mask_([0-9]+)_([0-9]+)(\.txt.*)$ ]]; then
            local prefix="${BASH_REMATCH[1]}"
            local file_num="${BASH_REMATCH[2]}"
            local phase="${BASH_REMATCH[3]}"
            local suffix="${BASH_REMATCH[4]}"
            local sort_key=$((file_num * 10000 + phase))
            old_files+=("$f")
            sort_keys+=("$sort_key")
        elif [[ "$basename_f" =~ ^test_([0-9]+)x([0-9]+)_([0-9]+)_([0-9]+)(\.txt.*)$ ]]; then
            local m="${BASH_REMATCH[1]}"
            local n="${BASH_REMATCH[2]}"
            local file_num="${BASH_REMATCH[3]}"
            local phase="${BASH_REMATCH[4]}"
            local suffix="${BASH_REMATCH[5]}"
            local sort_key=$((file_num * 10000 + phase))
            old_files+=("$f")
            sort_keys+=("$sort_key")
        fi
    done < <(find "$dir" -maxdepth 1 -type f -name "$pattern_old" -print0 2>/dev/null || true)
    
    if [[ ${#old_files[@]} -eq 0 ]]; then
        log_info "━━━ $desc ━━━"
        log_info "  未找到需要重命名的文件"
        log_info "  已有最大编号: $max_existing"
        echo ""
        return
    fi
    
    # 排序：按 sort_key 排序文件
    declare -a sorted_indices=()
    while IFS= read -r idx; do
        sorted_indices+=("$idx")
    done < <(for i in "${!sort_keys[@]}"; do
        echo "${sort_keys[$i]} $i"
    done | sort -n | awk '{print $2}')
    
    log_info "━━━ $desc (${#old_files[@]} 个待重命名) ━━━"
    log_info "  已有最大编号: $max_existing"
    log_info "  新文件将从 $((max_existing + 1)) 开始编号"
    echo ""
    
    local new_id=$((max_existing + 1))
    local renamed=0
    local skipped=0
    
    for idx in "${sorted_indices[@]}"; do
        local old_file="${old_files[$idx]}"
        local dir_path=$(dirname "$old_file")
        local basename_f=$(basename "$old_file")
        
        # 构造新文件名
        local new_basename=""
        if [[ "$basename_f" =~ ^(.*)_QC_H_[0-9]+_[0-9]+(\.txt.*)$ ]]; then
            new_basename="${BASH_REMATCH[1]}_QC_H_${new_id}${BASH_REMATCH[2]}"
        elif [[ "$basename_f" =~ ^(.*)_mask_[0-9]+_[0-9]+(\.txt.*)$ ]]; then
            new_basename="${BASH_REMATCH[1]}_mask_${new_id}${BASH_REMATCH[2]}"
        elif [[ "$basename_f" =~ ^test_([0-9]+)x([0-9]+)_[0-9]+_[0-9]+(\.txt.*)$ ]]; then
            new_basename="test_${BASH_REMATCH[1]}x${BASH_REMATCH[2]}_${new_id}${BASH_REMATCH[3]}"
        else
            log_warn "  跳过未识别格式: $basename_f"
            ((skipped++))
            continue
        fi
        
        local new_file="${dir_path}/${new_basename}"
        
        # 检查目标文件是否已存在
        if [[ -e "$new_file" ]]; then
            log_error "  [$new_id] 目标文件已存在，跳过: $new_basename"
            ((skipped++))
        else
            if ((APPLY)); then
                if mv "$old_file" "$new_file"; then
                    log_ok "  [$new_id] $basename_f → $new_basename"
                    ((renamed++))
                else
                    log_error "  [$new_id] 重命名失败，跳过: $basename_f"
                    ((skipped++))
                fi
            else
                echo -e "  ${YELLOW}[预览]${NC} [$new_id] $basename_f → $new_basename"
            fi
        fi
        
        ((new_id++))
    done
    
    echo ""
    if ((APPLY)); then
        log_info "  统计: ✓ 重命名 $renamed 个, ⊘ 跳过 $skipped 个"
    else
        log_info "  统计: 预览 $renamed 个, 跳过 $skipped 个"
    fi
    echo ""
}

# ==================== 执行重命名 ====================
# 1. H 矩阵
rename_files_smart "$MATRIX_DIR" "*_QC_H_*_*.txt*" "*_QC_H_*.txt" "H 矩阵"

# 2. Mask 矩阵
rename_files_smart "$MASK_DIR" "*_mask_*_*.txt*" "*_mask_*.txt" "Mask 矩阵"

# 3. Cycle 记录
rename_files_smart "$CYCLE_DIR" "test_*_*_*.txt*" "test_*_*.txt" "Cycle 记录"

# ==================== 总结 ====================
log_info "=========================================="
if ((APPLY)); then
    log_ok "重命名完成！"
    echo ""
    log_info "验证结果："
    echo ""
    
    if [[ -d "$MATRIX_DIR" ]]; then
        log_info "H 矩阵文件："
        find "$MATRIX_DIR" -maxdepth 1 -type f -name "*.txt" | sort -V | while read f; do
            echo "  → $(basename "$f")"
        done
        echo ""
    fi
else
    log_info "这是预览模式，未执行任何重命名"
    log_info "要执行重命名，请使用: $0 $M $N --apply"
fi
log_info "=========================================="

# 正常退出
exit 0