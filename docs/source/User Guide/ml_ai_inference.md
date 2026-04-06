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

**Status: In development, will be available soon**

DiskMINT-GARDEN is a grid of DiskMINT models covering a range of stellar masses, disk masses, gas-to-dust ratios, and characteristic radii. An XGBoost regressor trained on this grid will provide:

- **Forward inference:** fast prediction of disk observables (millimeter continuum flux, $\mathrm{C^{18}O}$ line flux) from input disk parameters
- **Inverse inference:** estimates of physical disk properties — gas mass, gas-to-dust ratio, and characteristic radius $R_c$ — from observed fluxes

Once published alongside the DiskMINT grid paper, users will be able to run inference with a single call:

```python
import diskmint.infer as infer

result = infer.from_observations(
    flux_mm=...,       # millimeter continuum flux [Jy]
    flux_c18o=...,     # C18O line flux [Jy km/s]
    mstar=...,         # stellar mass [M_sun]
    distance=...,      # distance [pc]
)
```

Full documentation, API reference, and example notebooks reproducing the paper figures will be added to the {doc}`../AI Features/ai_ref_index` section at that time.
