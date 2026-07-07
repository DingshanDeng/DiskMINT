# ML & AI Assistant

DiskMINT supports two categories of ML and AI tools: an **AI assistant skill** for guided installation, modeling, and support, and a **machine-learning inference API** for fast disk parameter estimation. This section covers both, and provides the structured reference files used by the AI assistant.

---

## DiskMINT-GARDEN -- ML Inference

DiskMINT-GARDEN (Grid of Astrochemical Radiative Disk EmissioN) is a grid of DiskMINT models spanning a range of stellar masses, disk masses, gas-to-dust ratios, and characteristic radii. The v1.7.0 package includes a predict-only XGBoost inference API:

- **Inverse inference:** estimates of physical disk properties (gas mass, gas-to-dust ratio, characteristic radius) from observed fluxes

Install with `pip install "diskmint[garden]"`, then call `diskmint.garden.infer.from_observations(...)` or `diskmint.garden.infer.from_dataframe(...)`. Required observed inputs are millimeter continuum flux density, $\mathrm{C^{18}O}$ integrated flux, distance, stellar mass, and the 90 percent dust radius. See {doc}`../user_guide/diskmint_garden` for API and domain-warning details.

---

## DiskMINT-Nursery — AI Assistant Skill

[DiskMINT-Nursery](https://github.com/DingshanDeng/DiskMINT-Nursery) is a companion skill for supported AI coding agents, for example [Claude Code](https://claude.ai/code) and [OpenAI Codex CLI](https://developers.openai.com/codex/cli). Once installed, it activates when you explicitly mention DiskMINT or `import diskmint` in your conversation — intentionally narrow, so it does not interfere with users of other thermochemical codes (DALI, ProDiMo, etc.). Some assistants may also support direct skill invocation. It provides three capabilities:

**Feature 1 — Installation & Onboarding**
Check and configure your full DiskMINT environment: Python package, Fortran chemistry network, RADMC-3D, gfortran version, and `DISKMINT_BIN_DIR`. Handles platform-specific issues (ARM Mac, legacy gfortran) and prints copy-paste commands for anything that needs `sudo`.

**Feature 2 — Runtime Assistant**
Answer questions about DiskMINT, help set up and run models, and interpret outputs. The skill reads the structured reference files below to ground its answers in your actual installation — parameters, file formats, pipeline steps, and known failure modes.

**Feature 3 — Support Escalation**
Diagnose unresolved errors against a known error reference, collect environment details, and draft a support email or GitHub issue report when a problem cannot be fixed automatically.

See {doc}`nursery_tutorial` for full installation and usage instructions.

---

## AI Agent Reference

The pages below are compact, structured reference files used by DiskMINT-Nursery to look up information during a session. They are also useful as quick-lookup pages for human readers — covering every parameter, the full pipeline, output file formats, and known errors.

::::{grid} 1 1 1 1

:::{grid-item-card} DiskMINT-Nursery Tutorial
:link: nursery_tutorial
:link-type: doc

Step-by-step guide to installing and using the DiskMINT-Nursery AI assistant skill with supported agents, including example workflows for Claude Code and Codex — covering installation & onboarding, runtime assistance, and support escalation.
:::

:::{grid-item-card} Installation Reference
:link: install_reference
:link-type: doc

All external tools, version requirements, install commands, environment variables, and platform-specific notes. Start here to set up or verify your environment.
:::

:::{grid-item-card} Parameters Reference
:link: parameters_reference
:link-type: doc

Every `Parameters` class attribute: physical meaning, units, default value, valid range, and interactions with other parameters.
:::

:::{grid-item-card} Workflow Reference
:link: workflow_reference
:link-type: doc

Step-by-step pipeline description — what each stage calls, what files are read and written, and how the stages connect.
:::

:::{grid-item-card} Output Format Reference
:link: output_format_reference
:link-type: doc

Directory structure of model output, file formats, column definitions for `.chem` files, and how to load each file type in Python.
:::

:::{grid-item-card} Error Reference
:link: error_reference
:link-type: doc

Known errors, root causes, and fixes — covering RADMC-3D crashes, chemistry convergence failures, gfortran issues, and common configuration mistakes.
:::

:::{grid-item-card} Prepared Prompts
:link: prompts_reference
:link-type: doc

Ready-to-use prompts for each DiskMINT-Nursery feature — copy into your AI assistant to get guided help for installation, model runs, and error diagnosis.
:::

::::

```{toctree}
:hidden:
:maxdepth: 1

nursery_tutorial
install_reference
parameters_reference
workflow_reference
output_format_reference
error_reference
prompts_reference
```
