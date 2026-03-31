#!/bin/zsh
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/../../../../" && pwd)
SRC_DIR="$REPO_ROOT/IBEX/src"
DST_DIR="$SCRIPT_DIR/include"

mkdir -p "$DST_DIR"

for header in \
  alloc.h \
  finite_lib.h \
  intio.h \
  ldpc_matrix_inverse.h \
  mod2convert.h \
  mod2dense.h \
  mod2sparse.h \
  rand.h \
  transceiver.h \
  vec_op.h
do
  cp "$SRC_DIR/$header" "$DST_DIR/$header"
done

echo "Synced headers into $DST_DIR"
