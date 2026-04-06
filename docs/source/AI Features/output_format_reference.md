# Output Format Reference

Directory layout, file formats, and column definitions for all DiskMINT output files. Includes Python snippets for loading each file type.

For the pipeline that produces these files see {doc}`workflow_reference`. For the User Guide description see {doc}`../User Guide/build_your_own_model`.

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
├── COendgrid-{name}.chem                    Final C18O abundance grid ← main science output
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
| `log10(n_C18O / n_H)` | — | log₁₀ C18O abundance relative to H nuclei |

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

Or via radmc3dPy (if installed):

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
