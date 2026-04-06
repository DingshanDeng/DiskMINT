# User Guide

This section covers the more advanced parts of the DiskMINT workflow, from building a model for a new target to running large grids and interpreting the physical assumptions behind the code.

It is intended for users who have already completed the {doc}`../Quick Start/quick_start_index` section or run one of the worked {doc}`../examples/example_index`.

---

## What You Will Find Here

- How to adapt DiskMINT to a new science target with your own parameter file, dust opacities, and stellar spectrum
- How to scale from a single model to many models on an HPC cluster
- How the key physical ingredients in DiskMINT connect to the model outputs
- How DiskMINT can be used for surface-density fitting and future ML-assisted inference workflows

---

## Guide Map

::::{grid} 1 1 2 2
:::{grid-item-card} Build Your Own Model
:link: build_your_own_model
:link-type: doc

Start from the example workflow, prepare a parameter CSV for a new target, and understand which model settings matter most in practice.
:::
:::{grid-item-card} Run on HPC
:link: hpc_guide
:link-type: doc

Use the example SLURM workflow to generate and submit large model grids efficiently on a cluster.
:::
:::{grid-item-card} Physics Background
:link: physics_background
:link-type: doc

Review the physical ingredients behind DiskMINT, including VHSE, dust settling, grain sizes, and surface density structure.
:::
:::{grid-item-card} Thermochemical Processes
:link: thermochemical_processes
:link-type: doc

See how vertical structure, isotope-selective photodissociation, and grain-surface chemistry work together to shape $\mathrm{C^{18}O}$ emission.
:::
:::{grid-item-card} Surface Density Fitting
:link: surface_density_fitting
:link-type: doc

Learn the profile-fitting workflow used to iteratively constrain radial dust and gas surface densities from resolved observations.
:::
:::{grid-item-card} ML / AI-Assisted Inference
:link: ml_ai_inference
:link-type: doc

Explore the developing AI and machine-learning tools around DiskMINT, including DiskMINT-Nursery and the DiskMINT-GARDEN inference direction.
:::
::::

---

## Suggested Reading

### If you are modeling a single target or are new to DiskMINT
start with {doc}`build_your_own_model`, then use {doc}`physics_background` and {doc}`thermochemical_processes` as reference pages while tuning the model.

### If you are preparing large parameter studies
read {doc}`hpc_guide` after the build-your-own-model page.

### If your goal is interpretation rather than setup
read {doc}`surface_density_fitting` and {doc}`ml_ai_inference` summarize the broader analysis directions that build on the core DiskMINT workflow.

```{toctree}
:hidden:
:maxdepth: 1

build_your_own_model
hpc_guide
physics_background
thermochemical_processes
surface_density_fitting
ml_ai_inference
