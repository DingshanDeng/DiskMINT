# Examples

This section provides worked examples of DiskMINT disk models for different stellar and disk configurations. Each example includes:

- A **parameter file** (`.csv`) describing the model setup
- A **main Python script** to run the full model pipeline
- **Jupyter notebooks** for post-processing and visualizing the outputs

The examples are located in the `examples/example_diskmint_models/` directory of the repository.

---

## Available Examples

::::{grid} 1
:::{grid-item-card} Disk Around RU Lup
:link: example_RULup
:link-type: doc

The best-fit model for RU Lup, a classical T Tauri star in Lupus, as presented in Deng et al. (2023). Demonstrates VHSE, the chemical network, and line radiative transfer with LIME.
:::
::::

::::{grid} 1
:::{grid-item-card} Disk Around a 0.5 M☉ Star
:link: example_0p5ms_star
:link-type: doc

A complete end-to-end example modeling a protoplanetary disk around a 0.5 solar-mass T Tauri star, with parameters inspired by HH30. Covers vertical hydrostatic equilibrium, dust settling, SED, and the chemical network.
:::
::::

```{toctree}
:hidden:
:maxdepth: 1

example_RULup
example_0p5ms_star
```