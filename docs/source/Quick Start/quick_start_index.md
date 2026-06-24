# Quick Start

Install DiskMINT and run your first disk model in minutes.

---

## Installation

**KEY REQUIREMENT**: `DiskMINT` relies on `RADMC-3D` (currently built on [`RADMC-3D` v2.0](https://github.com/dullemond/radmc3d-2.0)) to do radiative transfer and generate the thermal distributions.
So please make sure the `RADMC-3D` is properly installed according to the instructions on their website [their website](https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/index.html).

As of v1.6.0+, DiskMINT can be installed in three steps:

1. Download the code use the `git clone` or download this repo as zip from the latest [releases](https://github.com/DingshanDeng/DiskMINT/releases) into your local directory `Yourpath/DiskMINT/`.
2. Open terminal, go inside your path `cd Yourpath/DiskMINT/`, type `make install`. This should install both the `Python` and `Fortran` modules into your machine. 
3. Start Using! 

::::{grid} 1 1 2 2
:::{grid-item-card} Full Installation Guide
:link: Installation
:link-type: doc

Detailed instructions for both quick and manual installation, including Fortran compilation and platform-specific notes.
:::
:::{grid-item-card} External Requirements
:link: requirements
:link-type: doc

What external tools are needed (RADMC-3D, dust opacities, line RT code) and how to install each one.
:::
::::

---

## AI-Assisted Installation (Experimental)

```{warning}
DiskMINT-Nursery is experimental.
Install it from its own repository:
https://github.com/DingshanDeng/DiskMINT-Nursery
More information is described in {doc}`ML & AI Assistant <../AI Features/ai_ref_index>`
```

Try out [DiskMINT-Nursery](https://github.com/DingshanDeng/DiskMINT-Nursery), the new AI coding asistant that can walk you through the installation interactively.
This is an AI agent skill for [Claude Code](https://claude.ai/code),
[OpenAI Codex CLI](https://github.com/openai/codex), and compatible assistants.

Once the skill is installed, open your AI assistant in any directory and say:

> "Help me install DiskMINT from scratch."

The assistant will check each required tool in order (conda, DiskMINT Python package,
Fortran chemistry network, RADMC-3D, gfortran), install what it can automatically, and
print exact copy-paste commands for anything that requires elevated permissions.

After the installation, DiskMINT-Nursery can also walk you through setting up your first DiskMINT model.

---

## Run Your First Example

Follow a step-by-step walkthrough to run DiskMINT on the RU Lup example and understand the model outputs.

::::{grid} 1
:::{grid-item-card} Running Your First Model
:link: quickstart
:link-type: doc

Run the RU Lup example, understand the key model options, and learn how to set up dust opacities and the stellar spectrum for your own target.
:::
::::

```{toctree}
:hidden:
:maxdepth: 1

Installation
requirements
quickstart
```
