# Setup Script Walkthrough

This page explains `0-model_0p5ms_star_example.py` section by section. Reading it alongside
the script gives you a clear mental model of how a DiskMINT run is structured, so you can
adapt the example to your own target.

The full script is at:
```
examples/example_diskmint_models/example_diskmint_model_0p5ms_star/0-model_0p5ms_star_example.py
```

---

## 1. Imports

```python
import os, sys, copy, glob, shutil
import numpy as np
import time

import diskmint.constants as const
import diskmint.model as model
import diskmint.disk_density as dd
import diskmint.execute as exe
```

The four `diskmint` imports cover everything you need:

| Module | Role |
|---|---|
| `diskmint.constants` | Physical constants (`const.ms`, `const.rs`, `const.mu`, …) |
| `diskmint.model` | `Parameters` class (reads/edits the CSV) and `Mint` model object |
| `diskmint.disk_density` | Density structure utilities (used internally by `Mint`) |
| `diskmint.execute` | `runmodel()` — top-level function that drives the full pipeline |

The script also loads `diskmint_utils`, a collection of analysis helpers that live in
`examples/example_utils/`. They are **not** part of the installed `diskmint` package, so they
must be added to `sys.path` manually:

```python
utils_dir = os.path.join(current_dir, "..", "..", "example_utils")
sys.path.append(utils_dir)
import diskmint_utils as utils
```

---

## 2. Directory Setup

```python
working_dir = current_dir
data_dir    = os.path.join(working_dir, 'data')
read_dir    = os.path.join(working_dir, 'input')
save_dir    = os.path.join(working_dir, 'output')
```

| Variable | Purpose |
|---|---|
| `working_dir` | Root of the example directory — where the script and CSV live |
| `data_dir` | RADMC-3D working directory — all model input/output files go here; can be a separate path on a large-storage disk |
| `read_dir` | Pre-computed input files (stellar spectra, dust opacity files) to copy in before running |
| `save_dir` | Directory where chemical output files are archived after each run |

`data_dir` is kept separate from `working_dir` because a single model run can produce
hundreds of megabytes of RADMC-3D output. Pointing `data_dir` to a scratch partition prevents
your home directory from filling up.

---

## 3. Model Boolean Flags

```python
bool_MakeDustKappa  = True   # recompute dust opacities with optool
bool_SED            = True   # compute the SED with RADMC-3D
bool_VHSE           = True   # solve vertical hydrostatic equilibrium iteratively
n_vhse_loop         = 20     # maximum VHSE iterations (convergence usually < 10)
bool_dust_settling  = True   # apply dust settling
bool_chemistry      = True   # run the chemical network
```

These flags let you skip expensive steps when iterating on a single component:

- Set `bool_MakeDustKappa = False` after the first run to reuse already-computed opacity
  files — dust opacity computation with `optool` is slow for many size bins.
- Set `bool_VHSE = False` to use a simple Gaussian vertical structure (useful for a quick
  sanity check before a full run).
- Set `bool_chemistry = False` to skip the Fortran chemical network and only get the
  thermal/SED output.

---

## 4. Reading the Parameter File

```python
para = model.Parameters()
para.read_parameters_from_csv(filename='model_parameters_0p5ms_star',
                               directory=working_dir,
                               extension='.csv')
```

`model.Parameters()` creates a parameter object pre-loaded with default values.
`read_parameters_from_csv` then overwrites those defaults with the values in your CSV.
Lines starting with `#` are treated as comments and ignored.

After loading, individual parameters can be overridden in the script without editing the CSV:

```python
para.edit_parameter("chemical_save_dir", new_value=save_dir)
```

This is useful when you want to sweep one parameter programmatically (see the model loop
below) while keeping the CSV as a stable baseline.

---

## 5. The Model Loop

```python
for mdisk_dust_t, pla_set_t, amax_set_t in zip([1.0e-4], [3.5], [0.1]):
    para.edit_parameter("mdiskd",          new_value=mdisk_dust_t)
    para.edit_parameter("ratio_g2d_global",new_value=ratio_g2d_t)
    para.edit_parameter("pla_dustsize",    new_value=pla_set_t)
    para.edit_parameter("amax_all",        new_value=amax_set_t)
    ...
```

The script wraps the entire run inside a `for` loop, even though this example only has one
model. This structure makes it trivial to run a small parameter grid: just extend each list
inside `zip(...)` with additional values. Each iteration produces a separate, uniquely-named
output directory under `output/`.

The three parameters varied here are:

| Parameter | Meaning |
|---|---|
| `mdiskd` | Total dust disk mass [g] — set relative to solar mass via `const.ms` |
| `ratio_g2d_global` | Global gas-to-dust mass ratio (M_gas / M_dust) |
| `pla_dustsize` | Power-law slope of the dust size distribution |
| `amax_all` | Maximum grain size [cm] |

---

## 6. Dust Opacity Setup

```python
mass_fraction_carbon    = 0.13
mass_fraction_silicate  = 1 - mass_fraction_carbon   # 0.87
rho_carbon    = 1.8   # g/cm³
rho_silicate  = 3.01  # g/cm³
porosity_t    = 0.0

rhobulk = (1 - porosity_t) / (mass_fraction_carbon / rho_carbon
                               + mass_fraction_silicate / rho_silicate)
```

This block computes the effective bulk density of a DIANA dust grain (87% silicate + 13%
amorphous carbon by mass) with no porosity. The formula is a harmonic mean weighted by mass
fractions — the standard mixing rule for composite grain densities.

```python
nr_dust_set = 20
para = utils.set_dust_parameters_to_nrdust1(
    para, amin_set=1e-6, amax_set=amax_set, nr_dust_set=nr_dust_set,
    pla_set=pla_set, rhobulk_set=rhobulk_set)
```

`set_dust_parameters_to_nrdust1` populates the dust size bins in `para`. With `nr_dust_set=20`,
DiskMINT splits the grain size distribution from `amin=1e-6 cm` to `amax` into 20
logarithmically-spaced bins, each with its own opacity file.

```python
utils.wrapper_optool_opac(
    data_dir=data_dir, dust_opac_name=para.dustopacname_1,
    amin=para.amin_1, amax=para.amax_1, pla=para.pla_dustsize,
    nr_dust=para.dust_spec_nr,
    mass_frac_1=mass_fraction_silicate, mass_frac_2=mass_fraction_carbon,
    porosity=porosity_t)
```

`wrapper_optool_opac` calls `optool` under the hood to compute Mie scattering and absorption
opacities for each size bin. The resulting `dustkappa_*.inp` files are written to `data_dir`
in RADMC-3D format. This step requires `optool` to be on your `PATH` (see
[Requirements](../../Quick Start/requirements)).

---

## 7. Model Naming Convention

```python
name_of_this_model = '%s_%s_t%i_mdust%.1ems_gtd%.0f_pla%.2f_amax%.2f' % (
    model_name_prefix, model_date, i_model,
    para.get_parameter_value('mdiskd'),
    para.get_parameter_value('ratio_g2d_global'),
    para.get_parameter_value('pla_dustsize'),
    para.get_parameter_value('amax_all'))
name_of_this_model = name_of_this_model.replace('.', 'p')
```

The name encodes all varied parameters. Example output:

```
diskmint_example_0p5ms_star_20260223_t0_mdust1p0e-04ms_gtd100_pla3p50_amax0p10
```

Dots are replaced by `p` so the name is safe as a directory or file name on all platforms.
This convention makes it easy to identify a model from its output directory name alone,
without reading the parameter file.

---

## 8. Creating the `Mint` Object

```python
mint = model.Mint(para, file_dir=data_dir)

mint.bool_MakeDustKappa  = bool_MakeDustKappa
mint.bool_SED            = bool_SED
mint.bool_VHSE           = bool_VHSE
mint.n_vhse_loop         = n_vhse_loop
mint.bool_dust_settling  = bool_dust_settling
mint.bool_chemistry      = bool_chemistry
mint.chem_code_dir       = chem_code_dir
```

`model.Mint` is DiskMINT's central object. It holds the parameter set (`para`) and all
run-time options. Assigning the boolean flags here (rather than passing them to a function)
keeps the interface explicit — you can inspect or change any flag between steps if needed.

`file_dir=data_dir` tells `Mint` where to read and write RADMC-3D files.

---

## 9. Running the Model

```python
exe.runmodel(mint, test_alliteration=False)
```

This single call drives the complete pipeline in sequence:

1. Write RADMC-3D input files (grid, density, stellar spectrum, dust opacities)
2. Run the VHSE iteration loop (`n_vhse_loop` maximum iterations, stops on convergence)
3. Apply dust settling (re-distributes grain size bins vertically)
4. Run RADMC-3D thermal Monte Carlo to get dust and gas temperatures
5. Compute the SED if `bool_SED=True`
6. Run the Fortran chemical network (`reducedRGH22`) to get C¹⁸O abundances

`test_alliteration=False` disables a diagnostic mode that prints each VHSE step in detail;
leave it `False` for normal runs.

Runtime is typically 1–3 hours on a workstation with 24 threads (dominated by RADMC-3D Monte
Carlo and the VHSE iterations).

---

## 10. Saving Chemical Output Files

```python
data_save_dir = os.path.join(working_dir, 'model_outputs', 'data_files')
os.makedirs(data_save_dir, exist_ok=True)

for file_t in ['COinitgrid-*.dat', 'COendgrid-*.chem']:
    os.system('mv ' + os.path.join(data_dir, file_t) + ' ' + data_save_dir)
```

After `runmodel()` finishes, the chemical output files are moved out of `data_dir` (the
RADMC-3D scratch space) into a permanent archive under `model_outputs/data_files/`. This
prevents them from being overwritten if you re-run the model with different parameters.

| File pattern | Contents |
|---|---|
| `COinitgrid-*.dat` | Initial chemical grid (H nuclei density, gas temperature, UV field) at the start of the chemical network |
| `COendgrid-reducedRGH22-*.chem` | Final chemical network output — C¹⁸O abundance per grid cell after the full reaction network |

These files are what the analysis notebooks read (see Notebooks 1 and 3).
