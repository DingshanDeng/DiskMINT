# Scalable DiskMINT-RADMC3D Workflow

## Overview

This workflow enables efficient, parallel execution of DiskMINT models across multiple parameter grids using SLURM array jobs.

**Key improvements over original approach:**
- ✅ Parallel execution (all models run simultaneously)
- ✅ Easy to extend to different stellar masses
- ✅ Modular grid definitions
- ✅ Better error tracking and monitoring
- ✅ Scalable from 12 to 1000+ models

## Directory Structure

```
project/
├── run_model.py              # Main model runner (replaces 0-model_btsettl_2p0Msolar.py)
├── generate_grid.py          # Parameter grid generator
├── run_single_model.slurm    # SLURM script for individual model
├── submit_all_grids.sh       # Master submission script
├── check_status.sh           # Monitoring script
├── input/                    # Input files (dustkappa, etc.)
├── data/                     # Working directory for RADMC3D
├── output/                   # Final model outputs
├── grids/                    # Generated parameter grids (CSV)
└── logs/                     # SLURM output/error logs
```

## Quick Start

### 1. Setup

```bash
# Make scripts executable
chmod +x *.sh *.py

# Create necessary directories
mkdir -p input data output grids logs

# Copy your input files to input/
cp /path/to/dustkappa* input/
cp /path/to/*.inp input/
cp /path/to/model_parameters_template.csv ./
```

### 2. Generate Parameter Grids

```bash
# Default grid (12 models: 3 masses × 4 g2d)
python generate_grid.py --output grids/grid_2p0msolar.csv --grid default

# Extended grid (60 models: 5 masses × 4 g2d × 3 Rtap)
python generate_grid.py --output grids/grid_2p0msolar_ext.csv --grid extended

# Test grid (4 models for quick testing)
python generate_grid.py --output grids/test_grid.csv --grid test
```

### 3. Submit Jobs

**Option A: Submit single grid**
```bash
sbatch --array=0-11 run_single_model.slurm \
    grids/grid_2p0msolar.csv 2.0 20250824_3 grid2p0msolar_btsettl
```

**Option B: Use master script (recommended)**
```bash
# Edit submit_all_grids.sh to configure your grids
# Then submit:
bash submit_all_grids.sh
```

### 4. Monitor Progress

```bash
# Check job status
bash check_status.sh

# Or use SLURM commands directly
squeue -u $USER
sacct -j <JOB_ID> --format=JobID,JobName,State,Elapsed,CPUTime

# Check specific model log
tail -f logs/diskmint_<JOB_ID>_<ARRAY_INDEX>.out
```

## Extending to Multiple Stellar Masses

Edit `submit_all_grids.sh` to add more grids:

```bash
# Add these lines in the submission section:
submit_grid 0.5 "default" "grid0p5msolar_btsettl"
submit_grid 1.0 "default" "grid1p0msolar_btsettl"
submit_grid 1.5 "default" "grid1p5msolar_btsettl"
submit_grid 2.0 "extended" "grid2p0msolar_btsettl_extended"
```

This will launch 4 separate job arrays, each processing its own parameter grid in parallel.

## Customizing Parameter Grids

### Method 1: Use grid generator with custom ranges

Edit `generate_grid.py` to add a new grid type:

```python
elif grid_type == 'my_custom_grid':
    ratio_mdisk_values = [1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3]
    ratio_g2d_values = [10., 50., 100., 200., 300.]
    rtap_values = [50., 100., 200., 300., 400.]
```

### Method 2: Create CSV manually

Create a CSV file with columns: `ratio_mdisk,ratio_g2d,rtap`

```csv
ratio_mdisk,ratio_g2d,rtap
1e-5,10.0,300.0
1e-5,30.0,300.0
1e-4,10.0,300.0
...
```

## Resource Management

### Adjusting job resources

Edit `run_single_model.slurm`:

```bash
#SBATCH --time=48:00:00          # Increase for larger models
#SBATCH --cpus-per-task=48       # Match to your node capabilities
#SBATCH --mem=128G               # Adjust based on model size
```

### Limiting concurrent jobs

```bash
# Submit with percentage limit (e.g., max 10 jobs at once)
sbatch --array=0-99%10 run_single_model.slurm ...

# Or submit in chunks
sbatch --array=0-24 ...   # First 25 models
sbatch --array=25-49 ...  # Next 25 models
```

## Output Organization

Models are saved as:
```
output/
├── diskmintmodel20250824_3_t0_grid2p0msolar_btsettl_..._parameters.dat
├── diskmintmodel20250824_3_t0_grid2p0msolar_btsettl_.../
│   ├── COinitgrid-*.dat
│   ├── COendgrid-*.chem
│   └── ...
└── ...
```

Model naming: `diskmintmodel{DATE}_t{ID}_{MARK}_mgas{Mgas}_mdust{Mdust}_gtd{G2D}_rc{Rtap}_alphav{alpha}`

## Troubleshooting

### Job fails immediately
```bash
# Check error log
cat logs/diskmint_<JOB_ID>_<ARRAY_INDEX>.err

# Common issues:
# - Missing input files → check input/ directory
# - Wrong conda environment → verify conda activate works
# - Parameter file not found → check path in run_model.py
```

### Out of memory
```bash
# Increase memory in SLURM script
#SBATCH --mem=128G  # or higher

# Or reduce nphot in run_model.py
para.edit_parameter("nphot", new_value=1e6)  # instead of 1e7
```

### Jobs pending
```bash
# Check queue status
squeue -u $USER -t PENDING

# Check partition limits
scontrol show partition <partition_name>
```

## Comparison: Old vs New Workflow

| Aspect | Old (Sequential) | New (Array Jobs) |
|--------|-----------------|------------------|
| Parallelization | None | Full (all models simultaneously) |
| Total wall time (12 models) | ~12 × 24h = 288h | ~24h (max single model) |
| Extensibility | Edit nested loops | Add line to submit script |
| Error handling | One failure kills all | Independent jobs |
| Monitoring | Single log file | Per-model logs |
| Scalability | Poor (100s of models) | Excellent (1000s of models) |

## Best Practices

1. **Start with test grid**: Verify setup with 2-4 models before full runs
2. **Use job arrays**: SLURM handles scheduling/dependencies automatically  
3. **Monitor early**: Check first few jobs complete successfully
4. **Save grids**: Keep CSV files for reproducibility
5. **Log everything**: job_submissions.log tracks all submitted grids

## Integration with Chemistry Workflow

This follows the same pattern as your chemistry scaling:
- Grid generation → CSV files
- SLURM array jobs → parallel execution
- Master script → easy multi-grid submission
- Status monitoring → track progress

You can now easily run 100s of DiskMINT models across different stellar masses, disk parameters, and stellar atmosphere models in parallel.
