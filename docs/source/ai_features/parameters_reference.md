# Parameters Reference

Compact lookup table for every configurable parameter in DiskMINT. Parameters are set in the model CSV file and read with `para.read_parameters_from_csv()`; individual values can be overridden at runtime with `para.edit_parameter(name, new_value=...)`.

**CSV column order:** `name, value, units, valuetype, description`

Lines beginning with `#` are ignored. See {doc}`../user_guide/build_your_own_model` for the full guide.

---

## Section 1 — RADMC-3D Setup

| Parameter | Type | Units | Typical value | Description |
|-----------|------|-------|---------------|-------------|
| `nphot` | float64 | — | `1e7` | Monte Carlo photon packages. More = smoother temperature; each decade costs ~10× runtime. |
| `nthreads` | float64 | — | 10 | CPU threads for RADMC-3D. Set to physical core count. |
| `scattering_mode` | float64 | — | 0 | 0 = none, 1 = isotropic, 2 = anisotropic. |
| `nr` | float64 | — | 30 | Radial grid cells. |
| `r_step_max` | float64 | au | 5.0 | Maximum radial step size. |
| `ntheta` | float64 | — | 100 | Polar (θ) grid cells. |
| `nphi` | float64 | — | 1 | Azimuthal grid cells (1 = axisymmetric). |
| `rin` | float64 | au | 0.05 | Inner grid edge. Push to dust sublimation radius for VHSE. |
| `rout` | float64 | au | 1000 | Outer grid edge. Must enclose all dust and CO gas. |
| `thetaup` | float64 | rad | 0.7 | Upper polar boundary (radians from pole). |

---

## Section 2 — Stellar Parameters

| Parameter | Type | Units | Description |
|-----------|------|-------|-------------|
| `fmodel_filename` | str | — | Stellar spectrum file in RADMC-3D `.inp` format. Must be in `data_dir`. |
| `tstar` | float64 | K | Stellar effective temperature (informational; spectrum file takes precedence). |
| `mstar` | float64 | M☉ | Stellar mass. |
| `rstar` | float64 | R☉ | Stellar radius. |
| `mdotacc` | float64 | M☉ yr⁻¹ | Mass accretion rate. |
| `G0Hab_set` | float64 | Habing | UV flux at the stellar surface. Controls photodissociation chemistry. Typical values: 2e9 (0.1 M☉) → 2e11 (2.0 M☉). |

---

## Section 3 — Disk Parameters

### 3.1 Masses

| Parameter | Type | Units | Description |
|-----------|------|-------|-------------|
| `mdiskd` | float64 | M☉ | Disk dust mass. |
| `ratio_g2d_global` | float64 | — | Global gas-to-dust mass ratio. Gas mass = `mdiskd × ratio_g2d_global`. |

### 3.2 Surface Density

Σ(r) ∝ (r/Rtap)^`pl_sufdens` × exp[−(r/Rtap)^`pl_tapoff`]

| Parameter | Type | Units | Typical | Description |
|-----------|------|-------|---------|-------------|
| `pl_sufdens` | float64 | — | −1.0 | Power law index. Viscous disk: −1 to −0.5. |
| `pl_tapoff` | float64 | — | 1.0 | Tapering exponent. Viscous solution: 2 + `pl_sufdens`. |
| `Rtap` | float64 | au | 100 | Characteristic (tapering) radius. |

### 3.3 Scale Height

Two modes controlled by `scaleheight_index`:
- **1** (DIANA): `hp = hp100 × (r / 100 au)^plh`
- **2** (DALI): `hp = hprcr × Rtap × (r / Rtap)^plh`

| Parameter | Type | Units | Typical | Description |
|-----------|------|-------|---------|-------------|
| `scaleheight_index` | float64 | — | 1 | 1 = DIANA-style, 2 = DALI-style. |
| `hp100` | float64 | au | 10 | Scale height at 100 au (mode 1). |
| `hprcr` | float64 | — | 0.1 | hp/Rtap at Rtap (mode 2). |
| `plh` | float64 | — | 1.15 | Flaring index. Typical range: 1.1–1.25. |

### 3.4 Dust Physics

| Parameter | Type | Units | Description |
|-----------|------|-------|-------------|
| `visc_alpha` | float64 | — | Turbulence parameter for Dubrulle dust settling. |
| `vel_frag` | float64 | cm s⁻¹ | Fragmentation threshold velocity. |
| `a_drift` | float64 | cm | Grain size above which radial drift is applied. |
| `radius_drift` | float64 | au | Outer radius for radial drift treatment. |
| `mdiskd_2` | float64 | M☉ | Secondary dust component mass (set to 0 to disable). |
| `R_temp_trans` | float64 | au | Radius inside which gas temperature is decoupled from dust. |
| `fact_Tgas_2_Tdust` | float64 | — | Factor by which gas temperature exceeds dust temperature inside `R_temp_trans`. |

### 3.5 Geometry

| Parameter | Type | Units | Description |
|-----------|------|-------|-------------|
| `incl` | float64 | degrees | Disk inclination (RADMC-3D convention). |
| `phi` | float64 | degrees | Disk position angle (RADMC-3D convention). |

---

## Section 4 — Dust Opacity Setup

| Parameter | Type | Units | Description |
|-----------|------|-------|-------------|
| `dustopacname_1` | str | — | Base name for `dustkappa_{name}_{i}.inp` files. |
| `dustopacname_2` | str | — | Base name for a second dust composition (set same as `_1` if unused). |
| `nr_dust_1` | float64 | — | Number of size bins for opacity 1. |
| `nr_dust_2` | float64 | — | Number of size bins for opacity 2 (0 = unused). |
| `dust_spec_nr` | float64 | — | Total dust species (= `nr_dust_1 + nr_dust_2`). |
| `amin_1` / `amin_2` | float64 | cm | Minimum grain size. |
| `amax_1` / `amax_2` | float64 | cm | Maximum grain size. |
| `amin_all` / `amax_all` | float64 | cm | Overall min/max grain size across all compositions. |
| `pla_dustsize` | float64 | — | Grain size distribution power law index (MRN default: 3.5). |
| `rhobulk` | float64 | g cm⁻³ | Bulk grain density. |

---

## Section 5 — Chemistry Setup

| Parameter | Type | Units | Description |
|-----------|------|-------|-------------|
| `chemical_save_dir` | str | — | Directory where model outputs are written. |
| `chemical_save_name` | str | — | Model name (used as output subdirectory and filename prefix). |
| `nr_cyl_LIME` | float64 | — | Radial grid size for LIME interpolation. |
| `nz_cyl_LIME` | float64 | — | Vertical grid size for LIME interpolation. |

---

## Runtime-Only Parameters (set in Python, not CSV)

These are attributes of the `Mint` object set in the run script, not in the CSV:

| Attribute | Default | Description |
|-----------|---------|-------------|
| `bool_MakeDustKappa` | `True` | If True, compute `aave`, `ndsd`, `fracs`, `fracs_numb` from CSV grain parameters. If False, read pre-computed `.inp` files from `data_dir`. |
| `bool_SED` | `False` | Compute the SED with `radmc3d sed`. |
| `bool_VHSE` | `True` | Run the VHSE iteration loop. |
| `n_vhse_loop` | 10 | Maximum number of VHSE iterations. Convergence may exit early. |
| `bool_dust_settling` | `True` | Apply Dubrulle dust settling within the VHSE loop. |
| `bool_chemistry` | `True` | Run the reduced chemical network after VHSE. |
| `bool_savemodel` | `True` | Copy final model files to `chemical_save_dir/chemical_save_name/`. |
| `chem_code_dir` | auto | Path to the compiled chemistry binary directory. Set via `DISKMINT_BIN_DIR` or automatically resolved. |
