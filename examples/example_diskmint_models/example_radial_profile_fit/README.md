# Radial Profile Fitting Example

This example shows how to install an externally derived radial profile including
1. a dust surface density table
2. a radially-varying gas-to-dust ratio 
into a DiskMINT model, instead of relying on the built-in power-law surface density profile.

Two scripts are provided:

- **`0-model_radial_profile_fit_simple.py`** — a flat, single-run script with
  no command-line arguments. Start here.
- **`1-model_radial_profile_fit_advanced.py`** — an argparse-driven script
  with two named modes and CLI knobs for VHSE iterations, dust settling, and
  skipping chemistry. Use this once you're comfortable with the workflow and
  want to compare a constant vs. radially-varying gas-to-dust ratio.

Both scripts read the same `radial_profile_fit_parameters.csv` and share the
same `input/` and `model_inputs/` folders.

## Directory Structure

```text
example_radial_profile_fit/
├── 0-model_radial_profile_fit_simple.py
├── 1-model_radial_profile_fit_advanced.py
├── compare_radial_profiles.py
├── radial_profile_fit_parameters.csv
├── README.md
├── input/
│   ├── BTSettl_1p0Msolar_1pc_um_ergpcm2hzs.inp
│   ├── dustkappa_optool_20250122_*.inp
│   ├── aave.inp
│   ├── fracs.inp
│   ├── fracs_numb.inp
│   └── ndsd.inp
├── model_inputs/
│   ├── sigma_reference_template.dat
│   └── ratio_g2d_reference_template.dat
└── target_profiles/
    ├── example_continuum_profile.csv
    └── example_c18o_profile.csv
```

## Quick Start (Simple Script)

```bash
python -u 0-model_radial_profile_fit_simple.py 2>&1 | tee -a output_radial_profile_fit_simple.log
```

This reads `radial_profile_fit_parameters.csv`, installs
`model_inputs/sigma_reference_template.dat` as `sigmad_ref`, keeps
`ratio_g2d_global` constant, and runs with dust settling enabled
(`bool_dust_settling = True`) by default. Results are saved under
`output/example_radial_profile_fit_simple/`. There are no flags to set —
edit the variables near the top of the script directly if you want to change
run options (VHSE iterations, dust settling, chemistry on/off, etc.).

## Advanced / Two-Mode Script

```bash
# Baseline model with a constant gas-to-dust ratio
python -u 1-model_radial_profile_fit_advanced.py --mode model_A_constant_g2d 2>&1 | tee -a output_radial_profile_fit.log

# Radial fitting template: dust surface density AND gas-to-dust ratio both from reference tables
python -u 1-model_radial_profile_fit_advanced.py --mode model_B_radial_profile_fit 2>&1 | tee -a output_radial_profile_fit.log

# Quicker density/thermal test before running chemistry
python -u 1-model_radial_profile_fit_advanced.py --mode model_B_radial_profile_fit --skip-chemistry 2>&1 | tee -a output_radial_profile_fit.log
```

| Flag | Description |
|---|---|
| `--mode` | `model_A_constant_g2d` (constant `ratio_g2d_global`) or `model_B_radial_profile_fit` (adds a radially-varying `ratio_g2d` from `ratio_g2d_reference_template.dat`) |
| `--n-vhse-loop` | Maximum VHSE iterations (default 10) |
| `--dust-settling` / `--no-dust-settling` | Stable DiskMINT dust settling (default: enabled) |
| `--skip-chemistry` | Skip the chemical network for a faster thermal/density-only run |

## Shared Inputs

`sigma_reference_template.dat` and `ratio_g2d_reference_template.dat` are
plain two-column ASCII tables, `[radius_cm, value]`, with radius
monotonically increasing:

```text
radius_cm  sigma_dust_g_cm-2      # sigma_reference_template.dat
radius_cm  gas_to_dust_ratio      # ratio_g2d_reference_template.dat
```

`radial_profile_fit_parameters.csv` is shared by both scripts. Its
`mdiskd`/`ratio_g2d_global`/`pl_sufdens`/`pl_tapoff`/`Rtap` defaults already
match what the advanced script's `model_A` mode sets explicitly, and what the
simple script relies on implicitly (it does not re-set them).

## Comparing Model Output to Target Profiles

```bash
python compare_radial_profiles.py
```

This plots the bundled `model_inputs/` and `target_profiles/` tables as a
quick sanity check — it does **not** plot live RADMC-3D or chemistry output.
Output: `profile_comparison/radial_profile_fit_inputs.png`. Extend this
script with your own continuum-image or C18O line-radiative-transfer
products for a real model-vs-observation comparison.

## Notes and Gotchas

- Both scripts share `data/` as their RADMC-3D scratch directory. Run one
  script at a time — whichever runs second overwrites the intermediate files
  left by the first. Final saved outputs are safe: they're namespaced under
  `output/<chemical_save_name>/`, and the three possible names
  (`example_radial_profile_fit_simple`,
  `example_radial_profile_fit_model_A_constant_g2d`,
  `example_radial_profile_fit_model_B_radial_g2d`) never collide.
- Dust settling is enabled by default in both scripts
  (`bool_dust_settling = True`). Use `--no-dust-settling` with the advanced
  script to disable it; the simple script has no flags, so edit
  `bool_dust_settling` near the top of the file directly.
- Full write-up and pytest validation instructions:
  `docs/source/examples/example_radial_profile_fit.md`.
