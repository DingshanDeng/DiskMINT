# Output Format Reference

Directory layout, file formats, and column definitions for all DiskMINT output files. Includes Python snippets for loading each file type.

For the pipeline that produces these files see {doc}`workflow_reference`. For the User Guide description see {doc}`../user_guide/build_your_own_model`.

---

## Final Output Directory

When `bool_savemodel = True`, all results are written to:

```
{chemical_save_dir}/{chemical_save_name}/
```

Typical contents after a full run (`bool_VHSE = True`, `bool_chemistry = True`):

```
model_name/
├── amr_grid.inp                              RADMC-3D spherical grid
├── dust_density.inp                          Final dust density (RADMC-3D format)
├── dust_temperature.dat                      Final dust temperature (RADMC-3D format)
├── gas_density.inp                           Gas density (= dust × g2d ratio)
├── gas_temperature.inp                       Gas temperature
├── ratio_g2d_grid.dat                        Spatially varying g2d ratio
├── stars.inp                                 Stellar spectrum (RADMC-3D format)
├── dustopac.inp                              Opacity index file
├── dust_density_mdust_VHSE_{i}.dat          Density snapshot after iteration i
├── dust_temperature_mdust_VHSE_{i}.dat      Temperature snapshot after iteration i
├── gas_temperature_mdust_VHSE_{i}.inp       Gas temperature snapshot after iteration i
├── COinitgrid-{name}.dat                    Initial CO chemistry grid
├── COinitgrid-GSinit_{name}.dat             Same, with Gaussian vertical prior
├── COendgrid-{name}.chem                    Final $\mathrm{C^{18}O}$ abundance grid ← main science output
├── chemistry-reducedRGH22/                  Chemistry network binary outputs
└── {name}_parameters.csv                    Copy of the input parameter CSV
```

---

## Chemistry Output: `.chem` File

**Filename pattern:** `COendgrid-{chemical_save_name}.chem`

Plain-text file; one row per grid point on the cylindrical chemistry grid (`nr_cyl_LIME × nz_cyl_LIME` points).

| Column | Units | Description |
|--------|-------|-------------|
| `r` | cm | Cylindrical radius |
| `z` | cm | Height above midplane |
| `log10(n_H2 / n_H)` | — | log₁₀ H₂ abundance relative to H nuclei |
| `log10(n_{\mathrm{C^{18}O}} / n_H)` | — | log₁₀ $\mathrm{C^{18}O}$ abundance relative to H nuclei |

**Loading in Python:**

```python
import numpy as np

data    = np.loadtxt("COendgrid-mymodel.chem")
r_cm    = data[:, 0]
z_cm    = data[:, 1]
log_h2  = data[:, 2]
log_c18o = data[:, 3]

au = 1.496e13  # cm per au
r_au = r_cm / au
z_au = z_cm / au

# Abundance in linear units
abundance_c18o = 10**log_c18o
```

**Related files:**

| File | Description |
|------|-------------|
| `COinitgrid-{name}.dat` | Same 4-column format; initial (pre-network) CO abundance |
| `COinitgrid-GSinit_{name}.dat` | Same, using Gaussian vertical prior before VHSE |

### Chemical Output Analysis Workflow

The `.chem` raw columns are typically accessed through helper functions in `examples/example_utils/diskmint_utils.py`. This file lives in the DiskMINT git repository, not in the pip-installed package. Add its directory to `sys.path` before importing:

```python
import sys
sys.path.append("/your/path/to/DiskMINT/examples/example_utils")  # replace with your DISKMINT_REPO path
import diskmint_utils as utils

# Step 1 — load init + end grids together
df = utils.read_model(
    chem_save_name    = "mymodel_name",
    data_dir_init     = "output/mymodel_name/",
    data_dir_chemical = "output/mymodel_name/",
    have_CO2          = 'reducedRGH22'
)

# Step 2 — extract 2D arrays on the chemistry cylindrical grid
(rr_2D, zrr_2D), (rr_grid, zrr_grid), (tauv_star_2D, tauv_zup_2D), \
nH_2D, Tgas_2D, nH2_2D, nC18O_2D = utils.name_modelpara(df, nr_LIME=215, nz_LIME=200)
# rr_grid, zrr_grid  : 2D radius [au] and z/r
# nH_2D              : H nuclei number density [cm⁻³]  (from COinitgrid)
# Tgas_2D            : gas temperature [K]
# tauv_star_2D       : UV optical depth from the star
# nC18O_2D           : log₁₀(n_C18O / n_H)  (from COendgrid)

# Step 3 — compute C¹⁸O line luminosities and emitting layer
whole_data = utils.compute_emittinglayer(df, dv_set=0.0, nr=215, nz=200)
# prints C¹⁸O J=2-1 and J=3-2 luminosities in L☉

r_au_2D, z_au_2D, dlum21_2D, dlum32_2D, tauup_CO_2D = \
    utils.read_emittinglayer(wholedata=whole_data, nr=215, nz=200)
# dlum21_2D    : differential C¹⁸O J=2-1 luminosity per cell [L☉/cell]
# dlum32_2D    : differential C¹⁸O J=3-2 luminosity per cell [L☉/cell]
# tauup_CO_2D  : C¹⁸O vertical optical depth from above

# Step 4 — 4-panel diagnostic plot (nH, Tgas, τ_UV, nC18O + emitting layer)
import matplotlib.pyplot as plt
fig, axes = plt.subplots(1, 4, figsize=(16, 3))
utils.plot_density_emitting_layer(
    fig, axes,
    rr_grid, zrr_grid,
    nH_2D, Tgas_2D, tauv_star_2D, nC18O_2D, dlum21_2D
)
```

---

## RADMC-3D Density Files

**`dust_density.inp`** — dust density on the RADMC-3D (r, θ) spherical grid.

RADMC-3D binary/ASCII format. Load with the `modelgrid.readData()` utility:

```python
from diskmint import modelgrid

d = modelgrid.readData(dtemp=False, ddens=True, gtemp=False)
# d.rhodust  shape: (nr, ntheta, nphi, nr_dust_species)
# d.rhogas   not available directly; compute as sum(rhodust) × ratio_g2d
```

Or via radmc3dPy (see install notes — required for SED):

```python
import radmc3dPy.analyze as ra
d = ra.readData(dtemp=False, ddens=True)
rhodust = d.rhodust   # shape: (nr, ntheta, nphi, nr_dust_species)
```

---

## RADMC-3D Temperature Files

**`dust_temperature.dat`** — dust temperature on the (r, θ) grid.

```python
d = modelgrid.readData(dtemp=True, ddens=False, gtemp=False)
# d.dusttemp  shape: (nr, ntheta, nphi, nr_dust_species)
```

**`gas_temperature.inp`** — gas temperature (= dust temperature, or decoupled if `R_temp_trans > 0`).

```python
d = modelgrid.readData(dtemp=False, ddens=False, gtemp=True)
# d.gastemp  shape: (nr, ntheta, nphi, 1)
```

---

## Gas Density Reconstruction

Gas density is not stored directly — reconstruct it from dust density and the spatially-varying gas-to-dust ratio:

```python
import numpy as np
import diskmint.constants as const
import diskmint.modelgrid as modelgrid

d = modelgrid.readData(dtemp=True, ddens=True, gtemp=True)
rhodust     = d.rhodust                        # shape: (nr, ntheta, nphi, n_dust)
rhodust_all = rhodust.sum(axis=3)              # sum over dust species

ratio_g2d_grid = np.loadtxt("ratio_g2d_grid.dat").reshape(rhodust_all.shape)
ratio_g2d_grid[ratio_g2d_grid == np.inf] = np.nan

rhogas = rhodust_all * ratio_g2d_grid
rhogas[rhogas <= const.mu] = const.mu          # floor at mean molecular weight
```

---

## Coordinate System

RADMC-3D uses spherical coordinates (r, θ, φ). θ is co-latitude, measured from the z-axis (θ = 0 at the pole, θ = π/2 at the midplane). The opening angle z/r is:

```python
rc     = d.grid.x   # radii [cm]
thetac = d.grid.y   # co-latitude

rc_2D, theta_2D, _ = np.meshgrid(rc, thetac, d.grid.z, indexing='ij')
zr_2D = np.pi/2 - theta_2D   # opening angle z/r above midplane (dimensionless)
```

For 2D structure plots: x-axis = `rc / au` (radius in au), y-axis = `zr_2D` (z/r).

---

## Surface Density

Integrate the gas density vertically to get surface density Σ(r):

```python
au = 1.496e13  # cm per au
trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz

zz_2D = np.sin(zr_2D) * rc_2D             # z [cm]
# Integrate along θ axis (index 1); factor of 4π accounts for both hemispheres and azimuth
sigma_gas = 4.0 * np.pi * trapz(rhogas[:, :, 0], x=-zz_2D[:, :, 0], axis=1)
# sigma_gas  shape: (nr,)  units: g/cm²
```

---

## SED Output

**`spectrum.out`** — produced when `bool_SED = True`. Load with radmc3dPy (see install notes — required for SED):

```python
import radmc3dPy.analyze as ra
import numpy as np

s   = ra.readSpectrum("spectrum.out")
lam = s[:, 0]   # wavelength [μm]
fnu = s[:, 1]   # flux density [Jy] at 1 pc

distance_pc = 150.0
fnu_scaled = fnu / distance_pc**2           # scale to source distance [Jy]

nu   = 3e14 / lam                           # frequency [Hz]  (c in μm/s)
nufnu = nu * fnu_scaled                     # νFν [Jy·Hz = erg/s/cm²]
```

---

## VHSE Iteration Snapshots

The files `dust_density_mdust_VHSE_{i}.dat` and `dust_temperature_mdust_VHSE_{i}.dat` have the same format as their non-indexed counterparts. They are saved for diagnostic purposes; `i = 0` is the initial Gaussian, `i = N` is the last iteration before convergence exit.

To check convergence by eye, load successive snapshots and compare the vertical density profiles at a few radii.

---

## Parameter CSV Copy

A copy of the input parameter CSV is saved as `{chemical_save_name}_parameters.csv` inside the output directory. This ensures full reproducibility — parameters, grid, and outputs travel together.

---

## Model Naming Convention

The `chemical_save_name` (and therefore the output directory name) follows this pattern:

```
diskmintmodel{DATE}_t{ID}_{MARK}_mgas{Mgas}ms_mdust{Mdust}ms_gtd{G2D}_rc{Rtap}_alphav{alpha}
```

with `.` replaced by `p`. Example:
```
diskmintmodel20260405_t3_grid1p0msolar_mgas1p00e-02ms_mdust1p00e-04ms_gtd100_rc100p0_alphav5p0e-03
```

Parsing the name gives the key physical parameters without opening any file.
