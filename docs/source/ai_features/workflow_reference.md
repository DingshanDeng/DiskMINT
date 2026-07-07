# Workflow Reference

Step-by-step description of the DiskMINT pipeline: what each stage does, which files it reads and writes, and how stages connect. Intended as a quick lookup for the AI assistant and as a human-readable pipeline overview.

For the full narrative walkthrough see {doc}`../quick_start/quickstart`. For parameter details see {doc}`parameters_reference`.

---

## Pipeline Overview

```
Parameter CSV
    │
    ▼
[Stage 1] problem_setup()           ← writes RADMC-3D input files to data_dir
    │
    ▼
[Stage 2] radmc3d mctherm           ← computes initial dust temperature
    │
    ▼
[Stage 3] VHSE iteration loop       ← alternates density update ↔ temperature update
    │       (n_vhse_loop iterations, exits on convergence)
    │
    ▼
[Stage 4] chemistry_setup()         ← interpolates density/temperature onto chemistry grid
    │
    ▼
[Stage 5] runchemistry()            ← calls disk_main binary (Fortran chemistry network)
    │
    ▼
[Stage 6] savemodel()               ← copies final files to chemical_save_dir/chemical_save_name/
```

All stages are orchestrated by `exe.runmodel(mint)` in `diskmint.execute`.

---

## Stage 1 — Problem Setup

**Function:** `disk_density.problem_setup(mint)`

**What it does:**
- Writes the RADMC-3D grid (`amr_grid.inp`)
- Computes the initial Gaussian dust density (`dust_density.inp`)
- Computes the initial gas density (`gas_density.inp`) from dust × `ratio_g2d_global`
- Writes the gas temperature file (= dust temperature, or decoupled near star if `R_temp_trans > 0`)
- Writes the stellar spectrum (`stars.inp`) and opacity index (`dustopac.inp`)
- If `bool_MakeDustKappa = True`: computes `aave.inp`, `ndsd.inp`, `fracs.inp`, `fracs_numb.inp` from grain size parameters

**Files written to `data_dir`:**
`amr_grid.inp`, `dust_density.inp`, `gas_density.inp`, `gas_temperature.inp`, `stars.inp`, `dustopac.inp`, `wavelength_micron.inp`

---

## Stage 2 — Initial Thermal Calculation

**Command:** `radmc3d mctherm` (called via `os.system`)

**What it does:** Monte Carlo radiative transfer to compute dust temperature from the stellar radiation field.

**Files written:** `dust_temperature.dat`

**Convergence check:** If `dust_temperature.dat` is not produced on the first attempt, DiskMINT retries up to 3 times.

---

## Stage 3 — VHSE Iteration Loop

**Condition:** only runs if `bool_VHSE = True`

**Each iteration (`i_iter`) does:**

1. If `bool_dust_settling = True`: call `set_dust_settling()` → writes `dust_density_settled.new`, renamed to `dust_density.inp`
2. Run `radmc3d mctherm` again → updates `dust_temperature.dat`
3. Solve vertical hydrostatic equilibrium → writes `dust_density.new`, renamed to `dust_density.inp`; also writes `ratio_g2d_grid.new` → `ratio_g2d_grid.dat`

**Per-iteration snapshots saved:**
- `dust_density_mdust_VHSE_{i}.dat`
- `dust_temperature_mdust_VHSE_{i}.dat`
- `gas_temperature_mdust_VHSE_{i}.inp`

**Convergence criterion:** max(|Δρ_gas / ρ_gas|) < 0.05 across the grid. The loop also exits if differences start increasing (divergence detected) or if `n_vhse_loop` is reached.

**If `bool_VHSE = False` but `bool_dust_settling = True`:** one dust settling pass is applied to the initial Gaussian model, followed by a single thermal calculation.

---

## Stage 4 — Chemistry Grid Setup

**Function:** `disk_density.chemistry_setup(mint, save_name=...)`

**What it does:** Interpolates the converged density and temperature structure onto the cylindrical chemistry grid (`nr_cyl_LIME × nz_cyl_LIME`) and writes the initial abundance file.

**Files written:**
- `COinitgrid-{save_name}.dat` — initial CO abundance grid (input to chemistry)
- `COinitgrid-GSinit_{save_name}.dat` — same with Gaussian vertical prior (written before VHSE; used for comparison)

---

## Stage 5 — Chemistry Network

**Function:** `disk_density.runchemistry(mint, chem_code_dir=..., G0Hab_set=...)`

**What it does:** Calls the compiled Fortran binary `disk_main` from `chem_code_dir` (resolved from `DISKMINT_BIN_DIR` or auto-detected). The reduced chemical network (Ruaud et al. 2022) evolves H₂, CO, $\mathrm{C^{18}O}$, and associated species on the cylindrical grid.

**Key input parameter:** `G0Hab_set` — UV flux at stellar surface in Habing units. Scales the photodissociation rates.

**Files written:**
- `COendgrid-{save_name}.chem` — final abundance grid (see {doc}`output_format_reference`)

**Chemistry binary location:** `$DISKMINT_BIN_DIR/disk_main` (or `chemistry/bin/disk_main` relative to the DiskMINT package).

---

## Stage 6 — Save Model

**Function:** `disk_density.savemodel(mint, save_dir=..., model_dir=...)`

**Condition:** only runs if `bool_savemodel = True`

**What it does:** Copies all files from `data_dir` to `chemical_save_dir/chemical_save_name/`. If chemistry was run, also moves `COinitgrid-*.dat` and `COendgrid-*.chem` there, and renames the chemistry output directory to `chemistry-{network_name}`.

**Final output location:**
```
{chemical_save_dir}/{chemical_save_name}/
├── amr_grid.inp
├── dust_density.inp
├── dust_temperature.dat
├── gas_density.inp
├── gas_temperature.inp
├── dust_density_mdust_VHSE_*.dat       (per-iteration snapshots)
├── dust_temperature_mdust_VHSE_*.dat
├── COinitgrid-*.dat
├── COendgrid-*.chem
└── chemistry-reducedRGH22/             (chemistry network outputs)
```

---

## Key Entry Points in the Code

| Function | Module | Purpose |
|---|---|---|
| `exe.runmodel(mint)` | `diskmint.execute` | Top-level pipeline orchestrator |
| `model.Parameters()` | `diskmint.model` | Parameter container; read from CSV with `.read_parameters_from_csv()` |
| `model.Mint(para, file_dir=...)` | `diskmint.model` | Model object; holds parameters + runtime flags |
| `dd.problem_setup(mint)` | `diskmint.disk_density` | Stage 1: writes RADMC-3D inputs |
| `dd.chemistry_setup(mint, ...)` | `diskmint.disk_density` | Stage 4: writes chemistry grid |
| `dd.runchemistry(mint, ...)` | `diskmint.disk_density` | Stage 5: calls Fortran binary |
| `dd.savemodel(mint, ...)` | `diskmint.disk_density` | Stage 6: copies outputs |

---

## Runtime Flags — Stable vs. Experimental

**When generating or suggesting a run script, use only the stable flags below. Never include experimental flags in generated code.**

### Stable flags (safe to use)

| Flag | Default | Purpose |
|---|---|---|
| `mint.bool_MakeDustKappa` | `False` | Compute dust opacities from grain parameters |
| `mint.bool_SED` | `False` | Compute SED |
| `mint.bool_VHSE` | `False` | Solve vertical hydrostatic equilibrium |
| `mint.n_vhse_loop` | `10` | Max VHSE iterations (exits early on convergence) |
| `mint.bool_dust_settling` | `False` | Enable dust settling |
| `mint.bool_chemistry` | `False` | Run Fortran chemistry network |
| `mint.bool_savemodel` | `False` | Copy outputs to `chemical_save_dir` |
| `mint.chem_code_dir` | `'reducedRGH22'` | Chemistry network subdirectory name |

### Experimental flags (do NOT include in generated scripts)

The following flags are under active development and may produce unexpected or incorrect results. **Never include them in any generated or suggested script** unless the user explicitly asks about them by name:

- `mint.bool_temp_decouple`
- `mint.bool_dust_fragmentation`
- `mint.bool_dust_radial_drifting`
- `mint.bool_dust_inner_rim`
- `mint.bool_same_rc_as_radmc3d`

If the user asks about these flags, explain that they are experimental features not yet ready for general use and advise leaving them at their defaults (`False`).

---

## Minimal Run Script

```python
import diskmint.model as model
import diskmint.execute as exe

para = model.Parameters()
para.read_parameters_from_csv(filename="my_parameters.csv", directory=".", extension="")
para.edit_parameter("nthreads", new_value=10)

mint = model.Mint(para, file_dir="data/")
mint.bool_MakeDustKappa = True
mint.bool_VHSE          = True
mint.n_vhse_loop        = 10
mint.bool_chemistry     = True
mint.bool_savemodel     = True

exe.runmodel(mint)
```
