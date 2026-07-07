# DiskMINT-GARDEN

DiskMINT-GARDEN (Grid of Astrochemical Radiative Disk EmissioN) is the machine-learning inference API for estimating disk properties from observed millimeter continuum and $\mathrm{C^{18}O}$ fluxes. Unlike a full DiskMINT model run, it does not run RADMC-3D or the chemistry network. It uses XGBoost regressors trained on the DiskMINT-GARDEN model grid.

**Available since DiskMINT v1.7.0**

The current API provides inverse inference for:

- Gas mass
- Dust mass
- Gas-to-dust ratio
- Characteristic radius $R_c$

## Installation

Install the optional ML dependencies before using the API:

```bash
pip install "diskmint[garden]"
```

## Single-Target Inference

Run inference for one target with a single call:

```python
import diskmint.garden.infer as infer

result = infer.from_observations(
    flux_mm=120.0,       # millimeter continuum flux density [mJy]
    flux_c18o=850.0,     # C18O integrated line flux [mJy km/s]
    distance=140.0,      # distance [pc]
    mstar=0.8,           # stellar mass [M_sun]
    rdust_90=80.0,       # 90 percent dust radius [au]
    band="band6",        # band6: C18O J=2-1; band7: C18O J=3-2
)

print(result["Mgas_pred_Msun"], result["gtd_pred_dimless"], result["Rc_pred_au"])
```

## Table Inference

For tables, use `infer.from_dataframe(...)` with columns for continuum flux density, $\mathrm{C^{18}O}$ integrated flux, distance, stellar mass, and `Rdust_90_au`. The returned table includes:

- `Mgas_pred_Msun`
- `Mdust_pred_Msun`
- `gtd_pred_dimless`
- `Rc_pred_au`
- `is_outside_grid`
- `outside_features`
- `nn_dist`

Always check the grid-domain diagnostics. A flagged prediction means the observed inputs fall outside the model training domain and should be treated as an extrapolation warning, not as a precise physical constraint.

## Bundled Models

The package bundles two trained artifacts:

| Model | Line input | Validation note |
|---|---|---|
| `band6` | $\mathrm{C^{18}O}$ J=2-1 | test $R^2 \approx 0.9824$ |
| `band7` | $\mathrm{C^{18}O}$ J=3-2 | test $R^2 \approx 0.9835$ |

The models were trained from the DiskMINT-GARDEN grid wrapper artifacts prepared on June 23, 2026. Predictions outside the training domain are clipped to the grid limits and should be treated as extrapolation warnings.

For a minimal runnable example, see {doc}`../examples/example_garden`.
