# DiskMINT-Nursery Tutorial

[DiskMINT-Nursery](https://github.com/DingshanDeng/DiskMINT-Nursery) is an AI assistant skill for supported coding agents, for example [Claude Code](https://claude.ai/code) and [OpenAI Codex CLI](https://developers.openai.com/codex/cli). It helps you work with DiskMINT from first installation through model setup, output interpretation, and error diagnosis. It works by reading the structured AI reference files built into the DiskMINT documentation and using them to ground its answers in your actual installation.

---

## Prerequisites

Before installing the skill, make sure you have:

- DiskMINT installed and importable:
  ```bash
  python -c "import diskmint; print('OK')"
  ```
- A supported AI coding agent installed, for example [Claude Code](https://claude.ai/code) or [OpenAI Codex CLI](https://developers.openai.com/codex/cli)
- A terminal with `git` and `make`

---

## Installing the Skill

**Step 1.** Clone DiskMINT-Nursery:

```bash
git clone https://github.com/DingshanDeng/DiskMINT-Nursery.git
cd DiskMINT-Nursery
```

**Step 2.** Install the skill:

```bash
make install
```

This installs the skill for the supported local assistants. By default that includes Claude Code and Codex, copying the skill to `~/.claude/skills/diskmint-nursery/` and `~/.codex/skills/diskmint-nursery/`. To install for one platform only:

```bash
make install-claude   # Claude Code only
make install-codex    # Codex CLI only
```

**Step 3.** Restart your AI assistant to load the skill.

To verify the install, here are example checks for the currently supported assistants:

```bash
ls ~/.claude/skills/diskmint-nursery/SKILL.md    # Claude Code
ls ~/.codex/skills/diskmint-nursery/SKILL.md     # Codex CLI
```

---

## How the Skill Activates

The skill requires an **explicit DiskMINT reference** to activate. This is intentional: many disk modelers use other thermochemical codes (DALI, ProDiMo, ANDES, RAC2D, ...) and the skill should not interfere with their workflows.

The skill activates when your message explicitly:
- mentions **DiskMINT** by name, or
- references the Python package (`import diskmint`, `diskmint.model`, etc.)

The skill does **not** activate for generic disk modeling questions about RADMC-3D, VHSE, CO chemistry, or `.chem` files unless you also say DiskMINT. If you are a DiskMINT user, simply include "DiskMINT" in your message and the skill will engage.

Some assistants may also support direct skill invocation. For example, in Claude Code you can invoke the skill directly with:

```
/diskmint-nursery
```

In Codex, the reliable path is to mention DiskMINT explicitly in your message.

---

## Feature 1 — Installation & Onboarding

Use this feature to set up or verify your full DiskMINT environment.

**When to use:**
- First-time install on a new machine
- Setting up DiskMINT on an HPC cluster
- After a system update that may have broken gfortran or RADMC-3D
- If `disk_main` is not found or `DISKMINT_BIN_DIR` is not set

**How to invoke:**

Describe what you need in natural language:

> "Help me install and verify my full DiskMINT environment."

> "I'm setting up DiskMINT on a new Linux cluster. Walk me through the installation."

> "My DiskMINT chemistry network isn't compiling. Can you check my gfortran setup?"

**What the skill checks, in order:**

| Component | Verification |
|---|---|
| conda | `conda --version` |
| conda environment active | `conda info \| grep "active environment"` |
| DiskMINT Python package | `python -c "import diskmint; print('OK')"` |
| Fortran chemistry network | `$DISKMINT_BIN_DIR/disk_main` exists and runs |
| `DISKMINT_BIN_DIR` | variable is set and points to `chemistry/bin/` |
| RADMC-3D v2.0 | `radmc3d info` |
| gfortran 10+ | `gfortran --version` |
| optool | `optool --version` (optional) |

For missing components, the skill installs what it can automatically (no `sudo` required) and prints copy-paste commands for anything that needs elevated permissions.

---

## Feature 2 — Runtime Assistant

Use this feature when setting up a new model, adjusting parameters, running the pipeline, or interpreting outputs.

**When to use:**
- Setting up a parameter CSV for a new target
- Choosing dust opacity files or running `wrapper_optool_opac()`
- Understanding what a parameter controls and what value to use
- Loading and plotting a `.chem` output file
- Debugging a slow or non-converging VHSE run

**How to invoke:**

Ask naturally — the skill reads the appropriate reference file before answering:

> "What does `n_vhse_loop` control in DiskMINT? What value should I use?"

> "How do I configure dust opacities with optool for a new target? I don't have custom grain files."

> "How do I load and plot the $\mathrm{C^{18}O}$ abundance from a `.chem` file?"

> "My VHSE iterations aren't converging after 20 loops. What parameters should I check?"

> "Walk me through setting up a new DiskMINT model for a 1 M☉ T Tauri star."

**How the skill looks up information:**

For each question type, the skill reads the relevant reference file before answering, ensuring its response is grounded in your actual DiskMINT version:

| Question type | Reference file |
|---|---|
| Install / environment setup | {doc}`install_reference` |
| Parameters / CSV setup | {doc}`parameters_reference` |
| Pipeline / workflow | {doc}`workflow_reference` |
| Output files / `.chem` columns | {doc}`output_format_reference` |
| Errors / crashes | {doc}`error_reference` |

---

## Feature 3 — Support Escalation

Use this feature when an error cannot be resolved with standard fixes, or when you need to contact the DiskMINT author.

**When to use:**
- An error persists after trying the suggested fixes
- The error is not in the error reference
- You want to open a GitHub issue with complete context
- You need to send a support email to the DiskMINT author

**How to invoke:**

> "I'm getting this error running my DiskMINT model and can't resolve it: [paste error message]. Help me diagnose it."

> "My DiskMINT model won't run. Can you help me collect my environment info and draft a support email?"

**What the skill does:**
1. Checks {doc}`error_reference` for known fixes
2. Tries each suggested fix in sequence
3. If unresolved, collects your environment details (DiskMINT version, OS, gfortran version, full error log)
4. Drafts a complete support email addressed to dingshandeng@gmail.com, or formats the information for a GitHub issue at [DiskMINT/issues](https://github.com/DingshanDeng/DiskMINT/issues)

---

## Uninstalling

```bash
cd DiskMINT-Nursery
make uninstall
```

Or for a single platform:

```bash
make uninstall-claude
make uninstall-codex
```

---

## Notes

- The skill reads reference files from `docs/source/ai_features/` inside your local DiskMINT installation. If those files are not present (e.g., an older DiskMINT version installed without docs), the skill falls back to [diskmint.readthedocs.io](https://diskmint.readthedocs.io).
- The skill never modifies your parameter CSV without showing a diff and asking for confirmation first.
- The skill does not run `sudo` commands automatically — commands requiring elevated permissions are always printed for you to copy-paste.
- For bugs or feedback on the skill itself, open an issue at [DiskMINT-Nursery/issues](https://github.com/DingshanDeng/DiskMINT-Nursery/issues).
