#!/usr/bin/env bash
# 用法：
#   预览： ./seqrename.sh  fileA fileB fileC ...
#   执行： ./seqrename.sh --apply  fileA fileB fileC ...
# 说明：
#   会把  ..._X_Y.ext  改成  ..._1.ext, ..._2.ext ...
#   顺序就是你传入参数的顺序（选中的顺序）

set -euo pipefail

apply=0
if [[ "${1:-}" == "--apply" ]]; then
  apply=1
  shift
fi

if [[ "$#" -eq 0 ]]; then
  echo "用法: $0 [--apply] <选中的文件...>"
  exit 1
fi

n=1
for f in "$@"; do
  # 只处理存在的普通文件
  if [[ ! -f "$f" ]]; then
    echo "跳过(不是文件): $f"
    continue
  fi

  dir=$(dirname -- "$f")
  base=$(basename -- "$f")

  # 允许无扩展名，允许任意扩展名；匹配最后的 _数字_数字
  if [[ "$base" =~ ^(.*)_([0-9]+)_([0-9]+)(\.[^.]+)?$ ]]; then
    prefix="${BASH_REMATCH[1]}"
    ext="${BASH_REMATCH[4]:-}"
    new="${prefix}_${n}${ext}"

    # 避免覆盖
    if [[ -e "$dir/$new" ]]; then
      echo "目标已存在，跳过：$dir/$new"
    else
      if ((apply)); then
        mv -- "$f" "$dir/$new"
        echo "重命名：$f  ->  $dir/$new"
      else
        echo "[预览] mv -- '$f' '$dir/$new'"
      fi
      n=$((n+1))
    fi
  else
    echo "未匹配到末尾“_数字_数字”，跳过：$f"
  fi
done

