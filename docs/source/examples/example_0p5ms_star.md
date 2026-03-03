# Example: Disk Around a 0.5 $M_{\odot}$ Star

This example demonstrates a complete DiskMINT workflow for modeling a protoplanetary disk around a 0.5 solar-mass star, using physical parameters inspired by HH30. It covers the full pipeline: parameter setup, model execution, and post-processing analysis using Jupyter notebooks.

The example files are located at:

```
examples/example_diskmint_models/example_diskmint_model_0p5ms_star/
```

## Directory Structure

```
example_diskmint_model_0p5ms_star/
├── 0-model_0p5ms_star_example.py                # Main model script
├── model_parameters_0p5ms_star.csv              # Parameter file
├── 1-check_density_distribution_thermal.ipynb   # Notebook: thermal & density check
├── 2-check_radiative_transfer_SED.ipynb         # Notebook: SED check
├── 3-check_density_distribution_chemical.ipynb  # Notebook: chemical network check
├── input/                                       # Pre-computed input files
│   ├── BTSettl_0p5Msolar_1pc_um_ergpcm2hzs.inp  # Stellar spectrum (0.5 M☉)
│   ├── BTSettl_*.inp                             # Stellar spectra for other masses
│   ├── dustkappa_optool_20250122_*.inp           # Dust opacity files (20 size bins)
│   ├── aave.inp                                  # Average grain sizes
│   ├── fracs.inp                                 # Dust mass fractions
│   ├── fracs_numb.inp                            # Dust number fractions
│   └── ndsd.inp                                  # Normalized dust size distribution
└── lineop                                        # Line opacity data (for chemistry)
```

## Model Description

The model represents a disk around a star with the following properties:

| Parameter | Value |
|---|---|
| Stellar mass | 0.5 M☉ |
| Stellar effective temperature | 3802 K |
| Stellar radius | 1.9 R☉ |
| Disk dust mass | 1.0 × 10⁻⁴ M☉ |
| Gas-to-dust mass ratio | 100 |
| Tapering-off radius | 100 au |
| Disk outer grid edge | 1000 au |
| Dust grain size (a_min to a_max) | 1 × 10⁻⁶ to 0.1 cm |
| Dust size distribution slope | 3.5 |
| Dust composition | DIANA (silicate + carbon, via optool) |

The stellar spectrum is taken from the BT-Settl model grid. The `input/` directory also provides BT-Settl spectra for other stellar masses (0.1, 0.3, 0.7, 1.0, and 2.0 M☉) that can be substituted by editing `fmodel_filename` in the parameter file.

## Step 0: Running the Model

The main script is `0-model_0p5ms_star_example.py`. It is recommended to run it inside a `screen` session so it can continue if your connection drops:

```bash
# screen -S my_session
# conda activate diskmint_env
python -u 0-model_0p5ms_star_example.py 2>&1 | tee -a output.log
```

The script performs the following steps in order:

1. **Reads parameters** from `model_parameters_0p5ms_star.csv`.
2. **Copies input files** (stellar spectrum, pre-existing dust opacities, grain size files) to the working `data/` directory.
3. **Generates dust opacity files** using `optool` for DIANA dust (silicate + carbon mixture), with a power-law grain size distribution across 20 logarithmically spaced size bins.
4. **Sets up the `Mint` model object** and configures the run options.
5. **Executes the full model** via `exe.runmodel()`, which internally runs:
   - Vertical hydrostatic equilibrium (VHSE) iteration (up to 20 loops until convergence)
   - Dust settling
   - RADMC-3D thermal Monte Carlo and SED calculation
   - Chemical network (reducedRGH22)
6. **Saves chemical output files** (`COinitgrid-*.dat`, `COendgrid-*.chem`) to `model_outputs/data_files/`.

### Key Model Options

These boolean flags near the top of the script control which physics modules are active:

```python
bool_MakeDustKappa  = True   # Recompute dust opacity files with optool
bool_SED            = True   # Compute the SED with RADMC-3D
bool_VHSE           = True   # Solve vertical hydrostatic equilibrium iteratively
bool_dust_settling  = True   # Apply dust settling
bool_chemistry      = True   # Run the chemical network
n_vhse_loop         = 20     # Maximum VHSE iterations
```

### Running a Parameter Grid

The script is structured with a `for` loop so you can sweep over multiple parameter combinations (e.g., different dust masses, grain sizes, or gas-to-dust ratios) in a single run. Each iteration produces a uniquely named output subdirectory under `output/`. The model name encodes the key parameters for easy identification:

```
diskmint_similar_to_HH30_example_20260223_t0_mdust1p0e-04ms_gtd100_pla3p50_amax0p10
```

## Parameter File

`model_parameters_0p5ms_star.csv` is a human-readable CSV file that defines all model parameters. Lines starting with `#` are treated as comments. The columns are:

```
name, value, units, type, description
```

Parameters are organized into four sections:

- **RADMC-3D setup**: Grid resolution (`nr`, `ntheta`), photon count (`nphot`), number of CPU threads (`nthreads`), grid boundaries (`rin`, `rout`, `thetaup`)
- **Stellar parameters**: `mstar`, `rstar`, `tstar`, `mdotacc`, and the stellar spectrum filename (`fmodel_filename`)
- **Disk structure**: Surface density power law (`pl_sufdens`), tapering (`pl_tapoff`, `Rtap`), scale height (`hp100`, `plh`, `scaleheight_index`), dust settling parameters (`visc_alpha`, `vel_frag`), and inclination
- **Dust opacity**: Grain size range (`amin_1`, `amax_1`), number of species (`nr_dust_1`), size slope (`pla_dustsize`), bulk density (`rhobulk`), opacity file names (`dustopacname_1`, `dustopacname_2`)
- **Chemistry**: Save directory and model name (`chemical_save_dir`, `chemical_save_name`), LIME grid resolution (`nr_cyl_LIME`, `nz_cyl_LIME`)

You can edit this file to change any parameter before running the model, or override individual parameters directly in the script using:

```python
para.edit_parameter("parameter_name", new_value=...)
```

## Post-Processing Notebooks

After the model run completes, three Jupyter notebooks are provided for checking the outputs. In each notebook, update `work_dir` to point to your working directory and `i_model_select` to choose which model to inspect.

### Notebook 1: Density and Thermal Structure

**`1-check_density_distribution_thermal.ipynb`**

Reads the RADMC-3D output and visualizes:

- **Surface density profiles** (Σ) for the total dust, individual grain size bins, and the gas
- **2D gas number density maps** (log scale, shown with both logarithmic and linear r-axis)
- **2D dust and gas temperature maps**, including the smallest and largest grain size bins
- **Mass conservation check**: prints the integrated dust and gas masses and the resulting gas-to-dust ratio to verify the model setup is consistent with input parameters

### Notebook 2: Spectral Energy Distribution

**`2-check_radiative_transfer_SED.ipynb`**

Reads the RADMC-3D `spectrum.out` file and plots the SED as νF_ν vs. wavelength. It overlays:

- The input BT-Settl stellar spectrum
- The disk+star SED computed by RADMC-3D

The plot covers 0.1–5000 μm and is scaled to a fiducial distance of 150 pc. Observational data can be added in the designated cell (marked `# TBA`).

### Notebook 3: Chemical Network Output

**`3-check_density_distribution_chemical.ipynb`**

Reads the chemical network output files (`COinitgrid-*.dat` and `COendgrid-reducedRGH22-*.chem`) and produces a panel figure showing:

- **Total hydrogen number density** (n_H) in 2D
- **Gas temperature** (T_gas) in 2D
- **C¹⁸O abundance** after the chemical network
- **C¹⁸O J=2–1 emitting layer** (line luminosity per grid cell)

The C¹⁸O integrated line luminosities for the J=2–1 and J=3–2 transitions are also printed to the terminal.

## Input Files Reference

| File | Description |
|---|---|
| `BTSettl_0p5Msolar_1pc_um_ergpcm2hzs.inp` | BT-Settl stellar spectrum at 1 pc, in RADMC-3D format (wavelength in μm, flux in erg/cm²/Hz/s) |
| `BTSettl_*Msolar_*.inp` | BT-Settl spectra for 0.1, 0.3, 0.7, 1.0, and 2.0 M☉ stars |
| `dustkappa_optool_20250122_N.inp` | Dust opacity for grain size bin N (0–19), in RADMC-3D format |
| `dustkappa_optool_20250122_combined.inp` | Combined opacity file for all grain size bins |
| `aave.inp` | Average grain size for each size bin |
| `fracs.inp` | Dust mass fraction for each size bin |
| `fracs_numb.inp` | Dust number fraction for each size bin |
| `ndsd.inp` | Normalized dust size distribution |
| `lineop` | Line opacity table used by the chemical network |
