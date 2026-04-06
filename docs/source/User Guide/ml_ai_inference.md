# ML / AI-Assisted Inference

*This section is a placeholder. Full content will be added once the machine-learning inference tools are released alongside DiskMINT.*

---

DiskMINT is being extended with machine-learning tools to accelerate disk modeling:

**XGBoost surrogate model (coming in a future release):** An XGBoost regressor trained on the DiskMINT-GARDEN model grid provides fast inference of disk observables from input parameters — and, conversely, allows estimating physical disk properties (mass, gas-to-dust ratio, characteristic radius) directly from observed fluxes. Once published, users will be able to run inference with a single call to `diskmint.infer()` and reproduce the paper figures from [Deng et al. (2025)](https://ui.adsabs.harvard.edu/abs/2025ApJ...995...98D).

**AI agent skill (available now):** The [DiskMINT-Nursery](https://github.com/DingshanDeng/DiskMINT-Nursery) skill for Claude Code can assist with installation, model setup, and result interpretation. See the {doc}`../AI Features/ai_ref_index` section for details.

This page will be expanded with tutorials, API documentation, and example notebooks after the ML tools are released.
