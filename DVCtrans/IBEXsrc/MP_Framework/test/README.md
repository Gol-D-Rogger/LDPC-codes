# MP_Framework test harness

This directory provides a standalone sanitizer test framework for
`DVCtrans/IBEXsrc/MP_Framework`.

## Goal

- Build and run `MP_Framework/ldpc_codec.cpp` and `ldpc_c_model.c` under
  AddressSanitizer/UndefinedBehaviorSanitizer.
- Sweep representative IBEX matrix scenarios and exercise:
  - matrix load / crop
  - delta-preserving shift rebuild
  - IBEX encoder
  - clean-channel detector path
  - `BF_IBEX` decode
- Reproduce DV-facing memory bugs with a local, repeatable harness.

## Layout

- `include/`
  - copied headers from `IBEX/src`
- `sync_headers.sh`
  - refreshes copied headers from `IBEX/src`
- `Makefile`
  - builds the sanitizer binaries and run targets
- `mp_overflow_harness.cpp`
  - direct `ldpc_packet` scenario driver for matrix/config/encode/decode sweeps
- `mp_cmodel_harness.cpp`
  - DPI-style driver for `ldpc_c_model.c`
- `dv_official_lut_harness.cpp`
  - standalone driver for `DVCtrans/IBEXsrc/DVsrc/ldpc.h::f_create_include_files()`
- `include/svdpi.h`
  - minimal open-array shim so `ldpc_c_model.c` can be tested outside SV

## Build

From this directory:

```sh
make asan
```

Build the DPI-style harness:

```sh
make cmodel-asan
```

Build the official non-external-matrix LUT generator:

```sh
make dv-official-lut
```

## Run

Config-only sweep:

```sh
make run-config
```

Encode sweep:

```sh
make run-encode
```

Layer-decode sweep:

```sh
make run-layer
```

BF_IBEX sweep:

```sh
make run-bf
```

Full sweep:

```sh
make run-all
```

DPI/C-model regression for the exact `4096B + 320B/321B` parity cases:

```sh
make run-cmodel
```

Generate official DV LUT outputs from `DVsrc/ldpc.h`:

```sh
make run-dv-official-lut
```

This writes:

```sh
./output/ldpc_matrix.vh
./output/ldpc_matrix.h
```

Single custom byte-count case on the direct harness:

```sh
./build/mp_overflow_asan --mode bf_ibex --single 6 64 --user-bytes 4096 --parity-bytes 321
```

## Notes

- The harness defaults to matrix root `../../../../IBEX/ibex_matrix_flat_13rate`.
- `mp_overflow_harness.cpp` currently sweeps:
  - `M = 5..17`
  - `K = 64..67`
  - `parity_bytes = (M == 5) ? 320 : (M*64 - 1)`
- `mp_cmodel_harness.cpp` currently covers the DV-style exact byte cases:
  - `user_bytes = 4096`
  - `parity_bytes = 320`
  - `parity_bytes = 321`
- `dv_official_lut_harness.cpp` does not use external matrix files.
  It directly exercises the official built-in matrix family generator from
  `DVCtrans/IBEXsrc/DVsrc/ldpc.h`.
- Reproduced and fixed bug:
  - `ldpc_packet::ldpc_decoder()` had an early unconditional parity copy for
    non-`BF_IBEX` / non-`LAYER_G2` modes.
  - When `blk_len - info_len < hm_m` (shortened parity stream), that pre-copy
    read past `det_blk`.
  - The fix removes that early parity copy and lets the later
    shortening-aware parity import path handle all cases.
- The build uses:
  - `-fsanitize=address,undefined`
  - `-fno-omit-frame-pointer`
  - `-O1 -g`
