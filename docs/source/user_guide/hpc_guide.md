# Running DiskMINT on an HPC Cluster

This guide explains how to run DiskMINT model grids on a SLURM-based HPC cluster using the example workflow in `examples/example_diskmint_hpc/`. The same pattern scales from a handful of test models to several hundred production runs.

A complete, ready-to-run example is in:
```
DiskMINT/examples/example_diskmint_hpc/hpc_diskmint_radmc3d_models/
```

---

## Overview

The workflow has three components:

1. **Generate** — create a CSV describing the model grid (on your local machine)
2. **Submit** — submit a SLURM array job; each array element runs one model
3. **Monitor** — track job status and check logs

Each model runs as an independent SLURM job, so all models in a grid run simultaneously. A 60-model grid takes approximately the same wall time as a single model (~24 h).

---

## Directory Structure

```
project/
├── generate_grid.py            # Step 1: generate parameter grid
├── run_model.py                # Called by each SLURM array element
├── run_single_model.slurm      # SLURM job script
├── submit_all_grids.sh         # Master submission script
├── check_status.sh             # Job monitoring helper
├── model_parameters_Xmsolar.csv  # Base parameter file per stellar mass
├── input/                      # Dust opacity + stellar spectrum files
├── grids/                      # Generated parameter grid CSVs
└── logs/                       # Per-job SLURM output/error logs
```

---

## Step 1 — Prepare Input Files

Before running on the cluster, gather the files that every model will share:

- **Stellar spectrum** (`BTSettl_Xmsolar_*.inp`) — BT-Settl spectra for 0.1–2.0 M☉ are in `input/`
- **Dust opacity tables** (`dustkappa_*.inp`, `aave.inp`, `fracs.inp`, `fracs_numb.inp`, `ndsd.inp`) — the included optool tables are in `input/`
- **Base parameter CSV** (`model_parameters_Xmsolar.csv`) — one file per stellar mass; grid-variable parameters (`mdiskd`, `ratio_g2d_global`, `Rtap`) are overridden at runtime

---

## Step 2 — Generate the Parameter Grid (Local)

Run `generate_grid.py` on your local machine before uploading to the cluster. It produces a CSV where each row is one model.

```bash
# Default grid: 8 models (2 Mdisk/Mstar × 4 g2d, Rtap = 100 au)
python generate_grid.py \
    --output grids/grid_1p0ms_default.csv \
    --grid default \
    --mstar 1.0 \
    --model_dir /cluster/path/to/models/

# Extended grid: 80 models (5 masses × 4 g2d × 4 Rtap)
python generate_grid.py \
    --output grids/grid_1p0ms_extended.csv \
    --grid extended \
    --mstar 1.0 \
    --model_dir /cluster/path/to/models/

# Test grid: 2 models for quick verification
python generate_grid.py \
    --output grids/test_grid.csv \
    --grid test \
    --mstar 1.0
```

The CSV columns are `ratio_mdisk, ratio_g2d, rtap, model_mark_short, model_directory`. The `model_directory` column gives the absolute path where that model will write its `data/` and `output/` subdirectories.

### Custom grids

To define your own parameter ranges, edit `generate_grid.py` and add a new `grid_type` block:

```python
elif grid_type == 'my_survey':
    ratio_mdisk_values = [1e-4, 3e-4, 1e-3, 3e-3]
    ratio_g2d_values   = [30., 100., 300.]
    rtap_values        = [50., 150., 300.]
```

Or create a CSV manually with those three columns and use it directly with the SLURM script.

---

## Step 3 — Configure the SLURM Script

Open `run_single_model.slurm` and adjust the resource directives for your cluster:

```bash
#SBATCH --time=24:00:00          # Wall time per model (increase for VHSE-heavy runs)
#SBATCH --cpus-per-task=94       # Match to your node size; set nthreads= same value in run_model.py
#SBATCH --mem=128G               # Typical memory requirement; reduce nr/ntheta if too much
#SBATCH --partition=standard     # Your cluster partition name
#SBATCH --account=your_account   # Your allocation account
```

Also update the module load line (or conda activation) at the top of the script:

```bash
module load anaconda/2024.06     # Replace with your cluster's module name
# or:
# source ~/.bashrc && conda activate diskmint_env
```

---

## Step 4 — Submit Jobs

### Single grid

```bash
N=$(( $(wc -l < grids/grid_1p0ms_default.csv) - 1 ))
sbatch --array=0-$((N-1)) \
    run_single_model.slurm \
    grids/grid_1p0ms_default.csv \
    1.0 \
    20260405_0 \
    hpc_grid1p0msolar \
    model_parameters_1p0Msolar.csv
```

Arguments to `run_single_model.slurm`:
1. Grid CSV path
2. Stellar mass (M☉)
3. Model date identifier (used in output names)
4. Model mark / label
5. Base parameter CSV file

### Multiple grids with the master script

Edit the bottom section of `submit_all_grids.sh` to list your grids:

```bash
submit_grid 0.5 "default" "grid0p5msolar" "model_parameters_0p5Msolar.csv"
submit_grid 1.0 "extended" "grid1p0msolar" "model_parameters_1p0Msolar.csv"
submit_grid 2.0 "default" "grid2p0msolar" "model_parameters_2p0Msolar.csv"
```

Then run:

```bash
bash submit_all_grids.sh
```

Each call to `submit_grid` launches a separate SLURM array. All arrays run in parallel across stellar masses; within each array, all models run in parallel. Job IDs are logged to `job_submissions.log`.

---

## Step 5 — Monitor Jobs

```bash
# Overall status
squeue -u $USER

# Compact summary with elapsed time
sacct -j <JOB_ID> --format=JobID,JobName,State,Elapsed,CPUTime

# Per-model log (running)
tail -f logs/diskmint_<JOB_ID>_<ARRAY_INDEX>.out

# Quick check across all logs
bash check_status.sh
```

---

## Output Layout

Models are written to the `model_directory` column of the grid CSV:

```
model_directory/
├── output/
│   └── diskmintmodel{DATE}_t{ID}_{MARK}_mgas{Mgas}ms_mdust{Mdust}ms_gtd{G2D}_rc{Rtap}_alphav{alpha}/
│       ├── COinitgrid-*.dat
│       └── COendgrid-*.chem
└── data/              ← emptied after each job completes (saves disk space)
```

Model names follow the convention:
`diskmintmodel{DATE}_t{ID}_{MARK}_mgas{Mgas}ms_mdust{Mdust}ms_gtd{G2D}_rc{Rtap}_alphav{alpha}`

---

## Resource Guidelines

| Grid size | Typical resources per model |
|---|---|
| Test (2–4 models) | 4 cores, 32 GB, 2 h |
| Standard (10–60 models) | 32–48 cores, 64 GB, 24 h |
| Production (100+ models) | 48–96 cores, 128 GB, 24–48 h |

To limit how many jobs run at once (e.g., on shared partitions):

```bash
sbatch --array=0-59%10 run_single_model.slurm ...  # at most 10 simultaneous
```

---

## Troubleshooting

**Job fails immediately**
```bash
cat logs/diskmint_<JOB_ID>_<ARRAY_INDEX>.err
# Common causes: missing input files, wrong conda env name, bad module load line
```

**Out of memory**
```bash
# Increase in SLURM script:
#SBATCH --mem=256G
# Or reduce in parameter CSV: nphot = 1e6, nr = 20, ntheta = 60
```

**VHSE not converging**

Reduce `n_vhse_loop` to 20 for a first pass. 
Typically the VHSE takes about 10 iterations to converge, if it does not converge beyond 20 iterations.
You may need to check the grid set up.
For physically extreme models (very large or very small disks), the iteration may need `thetaup` or `rout` adjustments — see {doc}`build_your_own_model`.


**Jobs pending for a long time**
```bash
scontrol show partition <partition>   # check queue limits
# Or submit smaller batches: --array=0-19 then 20-39 etc.
```
