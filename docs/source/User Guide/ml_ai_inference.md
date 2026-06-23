# ML / AI-Assisted Inference

DiskMINT is being extended with machine-learning and AI tools to accelerate disk modeling and assist users throughout the workflow. Two tools are currently available or in development.

---

## DiskMINT-Nursery — AI Assistant Skill

**Status: Available (experimental)**

[DiskMINT-Nursery](https://github.com/DingshanDeng/DiskMINT-Nursery) is a companion skill for [Claude Code](https://claude.ai/code) and [OpenAI Codex CLI](https://developers.openai.com/codex/cli) that guides you through the full DiskMINT workflow — from installation to scientific results. It provides three features:

- Guided installation and environment verification
- Runtime assistance for model setup, parameter selection, and output interpretation
- Error diagnosis and support escalation

**See {doc}`../AI Features/nursery_tutorial` for a complete tutorial.**

---

## DiskMINT-GARDEN

**Status: beta in DiskMINT v1.7.0-beta**

DiskMINT-GARDEN (Grid of Astrochemical Radiative Disk EmissioN) is a grid of DiskMINT models covering a range of stellar masses, disk masses, gas-to-dust ratios, and characteristic radii. An XGBoost regressor trained on this grid will provide:

- **Inverse inference:** estimates of physical disk properties — gas mass, gas-to-dust ratio, and characteristic radius $R_c$ — from observed fluxes

Install the optional ML dependencies before using this beta API:

```bash
pip install "diskmint[garden]"
```

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

For tables, use `infer.from_dataframe(...)` with columns for continuum flux density, $\mathrm{C^{18}O}$ integrated flux, distance, stellar mass, and `Rdust_90_au`. The returned table includes `Mgas_pred_Msun`, `Mdust_pred_Msun`, `gtd_pred_dimless`, `Rc_pred_au`, and grid-domain diagnostics: `is_outside_grid`, `outside_features`, and `nn_dist`.

The beta package bundles two trained artifacts:

| Model | Line input | Validation note |
|---|---|---|
| `band6` | $\mathrm{C^{18}O}$ J=2-1 | test $R^2 \approx 0.9824$ |
| `band7` | $\mathrm{C^{18}O}$ J=3-2 | test $R^2 \approx 0.9835$ |

The models were trained from the DiskMINT-GARDEN grid wrapper artifacts prepared on June 23, 2026. Predictions outside the training domain are clipped to the grid limits and should be treated as extrapolation warnings, not as precise physical constraints.
