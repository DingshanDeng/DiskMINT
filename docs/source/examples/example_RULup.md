# Example: Disk Around RU Lup

This example reproduces the best-fit DiskMINT model for RU Lup, a classical T Tauri star in the Lupus star-forming region, as presented in [Deng et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023ApJ...954..165D). It demonstrates the full pipeline including VHSE, the chemical network, and line radiative transfer with LIME.

The example files are located at:

```
examples/example_diskmint_models/example_RULup/
```

## Directory Structure

```
example_RULup/
├── example_model_RULup.py # Main model script
├── example_RULup_parameters.csv # Parameter file
├── example_model_RULup_VHSE_bestfit_parameters.csv # Written-out best-fit parameters
├── example_model_RULup_VHSE_bestfit_parameters_sigmad_ref.dat # Reference surface density
├── data/ # RADMC-3D working directory
│ ├── RULupSpecModel_wuv.inp # Stellar + UV spectrum
│ ├── dustkappa_WD01_amax0p3cm_*.inp # Dust opacity files (20 size bins)
│ ├── aave.inp # Average grain sizes
│ ├── fracs.inp # Dust mass fractions
│ ├── fracs_numb.inp # Dust number fractions
│ ├── ndsd.inp # Normalized dust size distribution
│ ├── spectrum.out # RADMC-3D SED output
│ ├── ratio_g2d_grid.dat # Gas-to-dust ratio on the model grid
│ ├── colimemodel.c # LIME model file (C source)
│ ├── wrapper_LIME_example.py # Python wrapper for running LIME
│ ├── molecule_c18o.inp # C¹⁸O molecular data for LIME
│ ├── limegrid.dat # LIME grid file
│ ├── example_model_RULup_VHSE_bestfit_reducedRGH22/ # Chemical network output
│ └── images/ # LIME synthetic images output
└── LIME_example_RULup/ # LIME line radiative transfer example
```

## Model Description

RU Lup is a well-known classical T Tauri star with a massive, CO-rich disk. The DiskMINT best-fit model derives a self-consistent thermochemical disk structure and estimates the disk gas mass from C¹⁸O line emission.

| Parameter | Value |
|---|---|
| Stellar mass | 0.7 M☉ |
| Stellar effective temperature | 4060 K |
| Stellar radius | 2.5 R☉ |
| Mass accretion rate | 1 × 10⁻⁷ M☉/yr |
| Disk dust mass | 4.0 × 10⁻⁴ M☉ |
| Gas-to-dust mass ratio | 30 |
| Tapering-off radius | 63 au |
| Disk inclination | 18.8° |
| Disk outer grid edge | 500 au |
| Dust grain size (a_min to a_max) | 1 × 10⁻⁶ to 0.3 cm |
| Dust size distribution slope | 3.5 |
| Dust bulk density | 3.224 g cm⁻³ |
| Dust composition | WD01 (Weingartner & Draine 2001) |

For the full best-fit parameters and context, see [Deng et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023ApJ...954..165D).

## Running the Model

Run the main script from inside the `example_RULup/` directory: 
(It is recommended to run it inside a `screen` session so it can continue if your connection drops:)

```bash
# screen -S my_session
# conda activate diskmint_env
python -u example_model_RULup.py 2>&1 | tee -a output.log
```

The script performs the following steps in order:

1. **Reads parameters** from `example_RULup_parameters.csv`.
2. **Copies input files** (stellar spectrum, dust opacities) to the data/ working directory.
3. **Sets up the `Mint` model object** and configures the run options.
4.  **Executes the full model** via `exe.runmodel()`, which internally runs:
    - Vertical hydrostatic equilibrium (VHSE) iteration (up to 20 loops until convergence)
    - RADMC-3D thermal Monte Carlo and SED calculation
    - Chemical network (reducedRGH22)
5. Saves a copy of the model output to save_dir (set in the script).

### Key Model Options

```python
bool_MakeDustKappa  = True    # Compute dust size fractions from the opacity files
bool_SED            = True    # Compute the SED with RADMC-3D
bool_VHSE           = True    # Solve vertical hydrostatic equilibrium iteratively
bool_dust_settling  = False   # Dust settling is OFF — gas and dust are well-coupled
bool_chemistry      = True    # Run the chemical network
bool_savemodel      = True    # Save a copy of the model output to save_dir
n_vhse_loop         = 20      # Maximum VHSE iterations
```

Note that `bool_dust_settling = False` in this example. For RU Lup, the best-fit model assumes well-coupled gas and dust without settling.

### Parameter File

`example_RULup_parameters.csv` uses the standard DiskMINT `CSV` format (comma-separated, `#` for comments):

```
name, value, units, type, description
```

Key differences from a generic model setup:

- Larger grid: `nr = 170`, `ntheta = 200` to resolve the extended RU Lup disk
- WD01 dust: `dustopacname_1 = WD01_amax0p3cm` (Weingartner & Draine 2001 opacity model)
- Larger a_max: `amax_1 = 0.3` cm
- DSHARP bulk density: `rhobulk = 3.224` g cm⁻³
- Parameters can be overridden directly in the script using:

```python
para.edit_parameter("parameter_name", new_value=...)
```

### Line Radiative Transfer with LIME

After the DiskMINT model completes, the chemical output (CO abundances on a cylindrical grid) can be passed to LIME for non-LTE line radiative transfer. The data/ directory contains the necessary LIME input files:

- `colimemodel.c`: The LIME model file written in C. It reads the DiskMINT output grid (`limegrid.dat`) and interpolates the gas density, temperature, and C¹⁸O abundance at arbitrary positions requested by LIME.
- `wrapper_LIME_example.py`: A Python wrapper script that calls LIME with the correct settings (transitions, resolution, velocity range).
- `molecule_c18o.inp`: The C¹⁸O molecular data file (energy levels, Einstein coefficients) required by LIME.

To run LIME, ensure it is installed and then execute the wrapper:

```bash
python -u wrapper_LIME_example.py 2>&1 | tee -a output_lime.log
```

Synthetic channel maps will be saved to `data/images/`. These can then be convolved with an observational beam for comparison with ALMA data.

## Input Files Reference

| File | Description
|---|---|
| `RULupSpecModel_wuv.inp`	| RU Lup stellar spectrum including UV, in RADMC-3D format |
| `dustkappa_WD01_amax0p3cm_N.inp`	| WD01 dust opacity for grain size bin N (0–19), in RADMC-3D format |
| `aave.inp`	| Average grain size for each size bin |
| `fracs.inp`	| Dust mass fraction for each size bin |
| `fracs_numb.inp`	| Dust number fraction for each size bin |
| `ndsd.inp`	| Normalized dust size distribution |
| `colimemodel.c`	| LIME model source code for interpolating DiskMINT output |
| `molecule_c18o.inp`	| C¹⁸O molecular data for LIME |