# VREF Sweep Scripts

This folder contains two utilities:

- `submit_vref_sweep.py`: generate VREF-biased `.cnfg`/TOML files and submit parallel AutoFER jobs with `bsub`.
- `plot_vref_sweep.py`: read AutoFER `*_progress.xlsx` files and plot RBER-FER curves.
- `plot_vref_sweep.m`: MATLAB version of the plotting script; saves both `.fig` and `.png`.

## 1. Submit VREF Sweep

Example:

```bash
python3 scripts/vref_sweep/submit_vref_sweep.py \
  --template-cnfg IBEX/config/Ibex_hd_row11.cnfg \
  --template-toml <your_autofer_template.toml> \
  --exe <path_to>/ssd_fc \
  --auto-fer-dir scripts/auto_fer_eval \
  --out-dir output/vref_sweep \
  --queue regr_q
```

Behavior:

- Verifies `--exe` basename is `ssd_fc`.
- `--auto-fer-dir` points to the folder containing `auto_fer_eval.py` and helper scripts.
  You can also pass `--auto-fer-eval <path>/auto_fer_eval.py` directly.
- Patches line 34 / `VERF values` in the generated `.cnfg`.
- Generates one TOML per VREF case.
- Submits one outer `bsub` job per case:

```bash
bsub -q <queue> -J vref_<case> -o <out>/logs/<case>.autofer.log \
  python3 scripts/auto_fer_eval/auto_fer_eval.py adaptive --config <out>/toml/<case>.toml
```

Use `--dry-run` to print commands without submitting. Use `--no-submit` to only generate files.

## 2. VREF Combination Rule

Base VREF line:

```text
0 0.15 -0.15 0.3 -0.3 0.5 -0.5
```

For each common step (`0.025`, `0.05` by default), every soft pair independently chooses:

- `inward`
- `keep`
- `outward`

The all-keep case is excluded for each step. Baseline is included once. Total cases by default:

```text
1 baseline + 2 steps * (3^3 - 1) = 53
```

Generated metadata:

- `cases.csv`: case names, VREF values, generated paths, and expected progress xlsx path.
- `submit_all.sh`: exact `bsub` commands.
- `submitted_jobs.csv`: outer submission outputs.

## 3. Plot VREF Sweep

Example:

```bash
python3 scripts/vref_sweep/plot_vref_sweep.py \
  --input-root output/vref_sweep \
  --out-dir output/vref_sweep/plots \
  --min-fail-cw 10
```

Behavior:

- Reads `*_progress.xlsx` files under `--input-root`.
- Uses `RAW_BER` as x-axis and `LDPC_FER` as y-axis.
- `--min-fail-cw N` keeps only rows with `FAIL_CW >= N`.
- Uses `semilogy`.
- Uses `cases.csv` so legend entries are shown as `VREF=...`.
- The overlay plot uses a broad colormap so curves are easier to distinguish.
- All overlay curves use the same line width; only baseline is drawn with a dashed line.

Outputs:

- `plots/single/<case>.png`
- `plots/all_vref_overlay.png`

## 4. MATLAB Plot Version

MATLAB command:

```matlab
addpath('scripts/vref_sweep');
plot_vref_sweep('--input-root', 'output/vref_sweep', ...
                '--out-dir', 'output/vref_sweep/plots_matlab', ...
                '--min-fail-cw', '10');
```

Command line:

```bash
matlab -batch "addpath('scripts/vref_sweep'); plot_vref_sweep('--input-root','output/vref_sweep','--out-dir','output/vref_sweep/plots_matlab','--min-fail-cw','10');"
```

Outputs:

- `plots_matlab/all_vref_overlay.fig`
- `plots_matlab/all_vref_overlay.png`

The MATLAB version only generates the combined overlay plot. It uses the same visual rule as the Python overlay: all curves have the same line width, baseline is dashed, and colors are spread by colormap.
