# AI Features

*This section is under construction.*

This section provides compact, structured reference material for `DiskMINT`.
It is written to be machine-scannable and is used by [DiskMINT-Nursery](https://github.com/DingshanDeng/DiskMINT-Nursery) — the AI assistant skill for `DiskMINT` — to quickly locate information and help users.

All users will also find these pages useful as quick-lookup references. For narrative explanations and step-by-step tutorials, see the {doc}`../User Guide/user_guide_index`.

---

::::{grid} 1 1 2 2

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

install_reference
parameters_reference
workflow_reference
output_format_reference
error_reference
prompts_reference
```
