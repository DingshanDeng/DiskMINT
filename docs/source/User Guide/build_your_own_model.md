# Build Your DiskMINT Model

This guide walks you through adapting DiskMINT for a new target, from setting up your parameter file to interpreting the outputs. It assumes you have already run the RU Lup example from {doc}`../Quick Start/quickstart`.

---

## 1. Setting Up Your Parameter CSV

Every DiskMINT model is driven by a parameter CSV file. The file uses comma as the separator; lines beginning with `#` are ignored. The five columns are:

```
name, value, units, valuetype, description
```

The model scripts read this file with `para.read_parameters_from_csv(filename=...)`, and individual values can be overridden at runtime with `para.edit_parameter(name, new_value=...)`.

Template files for a range of stellar masses are provided in `examples/example_diskmint_hpc/hpc_diskmint_radmc3d_models/`.

### Parameter sections

**Section 1 — RADMC-3D Setup**

Controls the spatial grid and photon statistics for the Monte Carlo radiative transfer.

| Parameter | Typical value | Notes |
|-----------|---------------|-------|
| `nphot` | `1e7` | Photon packages. Increase for cleaner temperature structure; each ×10 costs ×10 runtime. |
| `nthreads` | 10–94 | CPU threads for RADMC-3D. Match to your hardware. |
| `nr` | 30 | Radial grid cells. |
| `ntheta` | 100 | Polar grid cells. |
| `rin` | 0.05 au | Inner edge — push inward to cover the dust sublimation radius. |
| `rout` | 1000 au | Outer edge — must enclose all dust and CO gas. |
| `thetaup` | 0.7 | Upper polar boundary (radians from the pole). |
| `scattering_mode` | 0 | 0 = no scattering (faster); 1 = isotropic; 2 = anisotropic. |

**Section 2 — Stellar Parameters**

| Parameter | Description |
|-----------|-------------|
| `fmodel_filename` | Stellar spectrum file in RADMC-3D format (`.inp`). Must be in `data_dir`. BT-Settl spectra for 0.1–2.0 M☉ are included in `examples/example_diskmint_hpc/input/`. |
| `mstar` | Stellar mass in solar masses. |
| `rstar` | Stellar radius in solar radii. |
| `mdotacc` | Mass accretion rate (solar masses per year). |
| `G0Hab_set` | UV flux at the stellar surface in Habing units. Controls photodissociation chemistry. |

**Section 3 — Disk Parameters**

*Masses*

| Parameter | Description |
|-----------|-------------|
| `mdiskd` | Disk dust mass in solar masses. |
| `ratio_g2d_global` | Global gas-to-dust mass ratio. Gas mass = `mdiskd × ratio_g2d_global`. |

*Surface density*

The surface density follows a tapered power law:
Σ(r) ∝ (r/Rtap)^`pl_sufdens` × exp[−(r/`Rtap`)^`pl_tapoff`]

| Parameter | Typical value | Description |
|-----------|---------------|-------------|
| `pl_sufdens` | −1.0 | Power law index. |
| `pl_tapoff` | 1.0 | Tapering index (viscous solution: 2 + `pl_sufdens`). |
| `Rtap` | 100 au | Characteristic (tapering) radius. |

*Scale height*

Two modes are supported (`scaleheight_index`):
- **Mode 1** (`scaleheight_index = 1`): `hp = hp100 × (r/100 au)^plh` — DIANA-style
- **Mode 2** (`scaleheight_index = 2`): `hp = hprcr × Rtap × (r/Rtap)^plh` — DALI-style

| Parameter | Description |
|-----------|-------------|
| `hp100` | Pressure scale height at 100 au (mode 1). |
| `hprcr` | hp/Rtap at Rtap (mode 2). |
| `plh` | Flaring index; typical 1.1–1.25. |

*Dust physics*

| Parameter | Description |
|-----------|-------------|
| `visc_alpha` | Turbulence parameter for dust settling (Dubrulle prescription). |
| `vel_frag` | Fragmentation threshold velocity (cm s⁻¹). |

**Section 4 — Dust Opacity Setup**

| Parameter | Description |
|-----------|-------------|
| `dustopacname_1` | Base name of the `dustkappa_*.inp` files (without the size index). |
| `dustopacname_2` | Optional base name for a second dust composition. The code supports two compositions, but the single-composition setup (`dustopacname_1` only) is the best-tested path. |
| `nr_dust_1` | Number of size bins. |
| `amin_1` / `amax_1` | Min/max grain size for this opacity (cm). |
| `amin_all` / `amax_all` | Overall min/max grain size across all dust populations. Keep these consistent with the per-composition ranges when generating size-distribution helper files. |
| `pla_dustsize` | Grain size distribution power law index (MRN: 3.5). |
| `rhobulk` | Bulk grain density (g cm⁻³). |

See [Section 2](#2-configuring-dust-opacities) for how to prepare the opacity files.

**Section 5 — Chemistry Setup**

| Parameter | Description |
|-----------|-------------|
| `chemical_save_dir` | Directory where chemistry outputs are written. |
| `chemical_save_name` | Model name string (also used for output filenames). |
| `nr_cyl_LIME` / `nz_cyl_LIME` | Grid size for LIME line RT (if used). |

---

## 2. Configuring Dust Opacities

DiskMINT requires pre-computed dust opacity tables in RADMC-3D `dustkappa_*.inp` format. A set of ready-to-use optool-computed tables is included in `examples/example_diskmint_hpc/input/`.

### Using the included opacity files

The example uses a 20-bin log-spaced size distribution (`nr_dust_1 = 20`, `amin_1 = 1e-6 cm`, `amax_1 = 0.1 cm`, `pla_dustsize = 3.5`) computed with optool. To use the same opacities:

1. Set `dustopacname_1 = optool_20250122` in your parameter CSV.
2. Copy all `dustkappa_optool_20250122_*.inp` files, plus `aave.inp`, `fracs.inp`, `fracs_numb.inp`, and `ndsd.inp`, into your model input directory (`data_dir` in the older scripts, or the directory passed as `file_dir="..."` when you build `model.Mint(...)`).

If you use two dust compositions, prepare one `dustkappa_{name}_{i}.inp` sequence per composition. In practice, most examples and the most stable workflow use just one composition.

### Computing your own opacities with optool

```bash
# Example
optool -c pyrmg70 0.87 -c org 0.13 -a 0.01 -na 1 -l 0.1 3000 200 -radmc
```

Run optool once per size bin (or write a loop), naming outputs `dustkappa_{name}_{i}.inp` for `i = 0` to `nr_dust_1 - 1` (smallest to largest).

Then set `bool_MakeDustKappa = True` in your run script to let DiskMINT compute `aave.inp`, `ndsd.inp`, `fracs.inp`, and `fracs_numb.inp` automatically from `amin_1`, `amax_1`, `amin_all`, `amax_all`, `pla_dustsize`, and `rhobulk`.

These are the current on-disk filenames used by the code. Older notes may refer to `frac.inp` or `frac_nb.inp`, but the maintained filenames are `fracs.inp` and `fracs_numb.inp`.

### UV spectrum (required for chemistry)

A stellar spectrum file in RADMC-3D format is required. BT-Settl spectra for standard stellar masses are included in `examples/example_diskmint_hpc/input/`. If UV observations are available for your target, fit a blackbody to estimate the UV component and splice it onto the photospheric spectrum. Scaling TW Hya's UV spectrum to your target luminosity is an acceptable first approximation.

---

## 3. Running the Model

### Script structure

A minimal model runner looks like this:

```python
import diskmint.model as model
import diskmint.execute as exe

para = model.Parameters()
para.read_parameters_from_csv(filename="my_parameters.csv", directory=".", extension="")

# Override parameters for this run
para.edit_parameter("nthreads", new_value=10)
para.edit_parameter("mdiskd",   new_value=1e-3)

# Configure model options
mint = model.Mint(para, file_dir="data/")
mint.bool_MakeDustKappa = True
mint.bool_VHSE          = True
mint.n_vhse_loop        = 20   # number of RT+hydrostatic iterations
mint.bool_chemistry     = True
mint.bool_savemodel     = True
mint.chem_code_dir      = "/path/to/DiskMINT/chemistry/reducedRGH22"

exe.runmodel(mint, test_alliteration=False)
```

For a complete working example, see `examples/example_diskmint_models/example_RULup/example_model_RULup.py`.

### What `runmodel()` produces

After a successful run, the `data_dir` contains standard RADMC-3D input/output files including density and temperature grids, plus DiskMINT-specific files:

| File pattern | Contents |
|---|---|
| `dust_density.inp` | Dust density on the (r, θ) grid |
| `dust_temperature.dat` | Dust temperature from RT |
| `gas_density.inp` | Gas density derived from dust + g2d ratio |
| `gas_temperature.inp` | Gas temperature (= dust temperature or decoupled near star) |
| `COinitgrid-*.dat` | Initial CO abundance grid (input to chemistry) |
| `COinitgrid-GSinit_*.dat` | Same with Gaussian vertical prior (pre-VHSE) |
| `COendgrid-*.chem` | Final $\mathrm{C^{18}O}$ abundance grid after chemistry network |

If `bool_savemodel = True`, a copy of the named model directory is written to `chemical_save_dir`.

### VHSE iteration

When `bool_VHSE = True`, the code alternates between:
1. Running RADMC-3D to compute the temperature structure
2. Solving the vertical hydrostatic equilibrium (VHS) to update the density structure

This loop runs `n_vhse_loop` times (10–20 is typical). Each RT pass takes roughly 5–30 minutes depending on `nphot` and `nthreads`. Convergence can be monitored by watching the scale-height profile stabilize between iterations.

---

## 4. Reading the Chemistry Network Outputs

### The `.chem` file

The final CO abundance grid is written to `COendgrid-{model_name}.chem`. It is a plain-text array with four columns:

| Column | Units | Description |
|--------|-------|-------------|
| `r` | cm | Cylindrical radius |
| `z` | cm | Height above midplane |
| `log10(n_H2 / n_H)` | — | log₁₀ H₂ abundance relative to H nuclei |
| `log10(n_$\mathrm{C^{18}O}$ / n_H)` | — | log₁₀ $\mathrm{C^{18}O}$ abundance relative to H nuclei |

Reading in Python:

```python
import numpy as np

data = np.loadtxt("COendgrid-mymodel.chem")
r_cm, z_cm, log_h2, log_c18o = data[:,0], data[:,1], data[:,2], data[:,3]

r_au = r_cm / 1.496e13
z_au = z_cm / 1.496e13
```

Here the r and z grids are the same used in the `COinitgrid-`, and the information of the number density of H nuclei can be read from the `COinitgrid-` file. 

### Using the outputs for line RT

The `.chem` grid can be passed directly to LIME or RADMC-3D for molecular line radiative transfer. The `nr_cyl_LIME` and `nz_cyl_LIME` parameters in the CSV control the cylindrical interpolation grid fed to LIME.

A complete LIME example is in `examples/example_diskmint_models/example_RULup/LIME_example_RULup/`.

For worked analysis using the chemistry outputs, see:
- {doc}`../examples/example_0p5ms_star` (0.5 M☉ example with chemistry)

---

## 5. Local vs. HPC

### Running locally

A single model with `bool_VHSE = True`, `n_vhse_loop = 10`, and `nphot = 1e7` takes roughly:

| Hardware | Wall time (approx.) |
|---|---|
| Laptop (4–8 cores) | 4–12 hours |
| Desktop / workstation (16–32 cores) | 1–4 hours |
| HPC node (94 cores) | <1 hour |

Set `nthreads` to the number of physical cores available. For exploratory work, reduce `nphot` to `1e6` or disable VHSE (`bool_VHSE = False`) for a fast sanity check.

### Running on an HPC cluster

For grid-scale runs (tens to hundreds of models), use the SLURM array job workflow in `examples/example_diskmint_hpc/`. Each model runs as an independent array element, so all models in a grid run simultaneously.

See the {doc}`hpc_guide` for step-by-step instructions including job scripts and a worked multi-grid example.
