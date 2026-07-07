# Examples

This section provides worked examples of DiskMINT disk models for different stellar and disk configurations. Each example includes:

- A **parameter file** (`.csv`) describing the model setup
- A **main Python script** to run the full model pipeline
- **Jupyter notebooks** for post-processing and visualizing the outputs

The examples are located in the `examples/example_diskmint_models/` directory of the repository.

---

## Available Examples

::::{grid} 1
:::{grid-item-card} Disk Around a 0.5 M☉ Star — Most Comprehensive Example
:link: example_0p5ms_star
:link-type: doc

The most complete worked example available. Models a protoplanetary disk around a 0.5 M☉ T Tauri star (parameters inspired by HH30), covering the full pipeline: VHSE, dust settling, SED, and the chemical network. Includes three annotated Jupyter notebooks for post-run analysis (density/thermal structure, SED, and chemical output).
:::
::::

::::{grid} 1
:::{grid-item-card} Disk Around RU Lup
:link: example_RULup
:link-type: doc

The best-fit DiskMINT model for RU Lup as published in [Deng et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023ApJ...954..165D). Demonstrates VHSE, the chemical network, and line radiative transfer with LIME.
:::
::::

::::{grid} 1
:::{grid-item-card} Fitting Continuum and C18O Radial Profiles
:link: example_radial_profile_fit
:link-type: doc

A target-agnostic workflow example for iteratively adjusting the dust surface density and gas-to-dust mass ratio to match continuum and $\mathrm{C^{18}O}$ radial profiles.
:::
::::

::::{grid} 1
:::{grid-item-card} DiskMINT-GARDEN Inference Quickstart
:link: example_garden
:link-type: doc

A lightweight example for `diskmint.garden`: load the bundled surrogate model and infer gas mass, dust mass, gas-to-dust ratio, and characteristic radius from observed continuum and $\mathrm{C^{18}O}$ fluxes. Includes a runnable smoke script and an annotated notebook.
:::
::::

```{toctree}
:hidden:
:maxdepth: 1

example_0p5ms_star
example_RULup
example_radial_profile_fit
example_garden
```
