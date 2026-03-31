# IBEX RTL LUT Quick Alignment Guide

This note is a fast semantic alignment guide for readers who already understand QC-LDPC matrix basics and only need to align RTL LUT meaning.

Scope:
- Focus on runtime RTL CN alignment semantics.
- Keep notation close to existing docs/code, with first-use aliases.
- Exclude validation-method discussion.

---

## 0. 30s Summary

- `WrapBase (base_r)` is the virtual point one `delta` step before the first active CPM in row `r`.
- Row shift chain is generated left-to-right:
  `S_r(n) = (WrapBase[r] + (n + 1) * delta[r]) mod Z`.
- `WrapBaseDelta (y_r)` is Definition A:
  minimal `k` such that `(last_element + k * delta) mod Z == 0`.
- RTL CN engine aligns CN to cycle-0 first, then rotates by one `delta` per active column, then applies row wrap:
  `wrap_shift = WrapBase + WrapBaseDelta * delta`.
- The common doc formula `(y_r + 1) * delta + base_r` is not in conflict with code:
  code already consumed one `delta` rotation in the per-column loop.

---

## 1. Notation Map (Doc <-> Intuitive Alias <-> Code)

| Concept | Doc name | Alias in this note | Code-side name |
|---|---|---|---|
| Row base seed | `WrapBase` | `base_r` | `rtl_view.wrap_base[i]` |
| Row delta step | `delta_r` | `delta_r` | `h_matrix.delta[i]` |
| Wrap delta-count | `WrapBaseDelta` | `y_r` | `rtl_view.wrap_base_delta[i]` |
| First active shift | `first_element` | `first_shift` | `rtl_view.first_shift[i]` |
| Cycle-0 aligned CN | `S'r0` | `cycle-0 aligned CN` | `rotate_cn_row(cn_row, first_shift, bits)` output |
| Row wrap rotation | `decoder_rotation` | `wrap_shift` | `wrap_base + wrap_base_delta * delta` |

---

## 2. Two Views That Often Get Mixed

### 2.1 Physical matrix order

Physical column order is still:
- left: payload columns
- right: parity columns

### 2.2 Anchor-centered signed-shift order

When formulas use:
`{j*delta, (j-1)*delta, ..., 0, -delta, ..., -x_r*delta}`
they describe offsets around an anchor (`shift=0`), not physical left/right placement.

Practical reading:
- negative side corresponds to one side of the anchor (often payload side in the chosen indexing),
- positive side corresponds to the other side.

So confusion usually comes from mixing:
- "matrix physical order"
- and "signed offset around anchor."

---

## 3. Sr Formula: Why It Becomes "Shift Then XOR"

Row syndrome:
`Sr = sum_k(Hr'_k * C'_k)`.

Each non-zero CPM multiply equals cyclic shift of subvector:
`Hr'_k * C'_k = C'_k >> shift_k` (mod `Z` ring).

Hence:
`Sr = sum_k(C'_k >> shift_k)` over GF(2) (bitwise XOR).

If `shift_k < 0`, then:
`C >> (-a) == C << a` (same cyclic ring).

This is why Sr formulas look like mixed `>>` and `<<` terms.

---

## 4. Cycle Alignment Logic (What RTL Actually Needs)

RTL wants one fixed datapath per cycle, so it rotates the row CN state such that "current active CPM" is always aligned to position 0.

For row `r`:

1. Initial alignment:
   `S'r0 = Sr << (j * delta_r)` (cycle-0 aligned CN)
2. Per-cycle step:
   `S'r(t+1) = S'rt >> delta_r`
3. Last active position:
   `S'r(j+x_r) = Sr << (-x_r * delta_r)`

Meaning:
- no per-column custom rotation logic;
- only one fixed rotate-by-delta operation each cycle.

---

## 5. Why Wrap Is Needed After Column Loop

After processing all active columns in one row-iteration:
- CN has already accumulated multiple in-loop rotations.
- Without wrap compensation, next iteration would start from wrong phase.

So row-end wrap applies:
`wrap_shift = base_r + y_r * delta_r` (mod `Z`).

In code:
- column loop already executes one `>> delta` after each active column;
- therefore the separate wrap expression does not carry the extra `+1`.

This matches the doc-level identity once total in-loop rotation is counted.

---

## 6. Minimal Worked Example (Compact)

Assume:
- `Z=16`, `delta=2`, `j=2`, `x_r=2`
- active shift offsets around anchor: `{+4,+2,0,-2,-4}`

Then:
- `Sr = C0>>4 + C1>>2 + C2 + C3<<2 + C4<<4`
- `S'r0 = Sr<<4` (align `C0`)
- `S'r1 = S'r0>>2` (align `C1`)
- `S'r2 = S'r1>>2` (align `C2`)
- ...

The pattern is always:
- first align to cycle-0
- then shift by `delta` each cycle
- then apply row wrap to return to next-iteration start phase.

---

## 7. Direct Code Touchpoints

- RTL view build and row parameters:
  `DVCtrans/IBEXsrc/MP_Framework/ibex_rtl_engine.cpp`
- Initial CN align (`first_shift`), per-column delta rotate, row-end wrap:
  `DVCtrans/IBEXsrc/MP_Framework/ibex_rtl_engine.cpp`
- Source matrix/RTL config generation context:
  `IBEX/src/ldpc_codec.cpp`
- Canonical semantics reference:
  `IBEX/doc/IBEX_Matrix_LUT.md`

---

## 8. Quick Checklist (When Reading Waveforms or Dumps)

- If CN phase looks off by a row-wise constant:
  check whether you compare raw `Sr` vs cycle-aligned `S'r0`.
- If wrap appears missing one `delta`:
  include the last in-loop `>> delta`.
- If payload/parity side interpretation feels inverted:
  confirm whether current formula uses physical order or anchor-signed order.

