# Installation Reference

Compact reference for setting up the full `DiskMINT` environment.
For a narrative walkthrough, see {doc}`../Quick Start/Installation` and {doc}`../Quick Start/requirements`.

---

## Checklist

Run these to check what is already installed before doing anything else:

```bash
conda info --envs                                          # list conda environments
echo $CONDA_DEFAULT_ENV                                    # check active env (empty or 'base' = no dedicated env)
python -c "import diskmint.model; print('DiskMINT OK')"   # DiskMINT Python package
echo $DISKMINT_BIN_DIR                                     # chemistry binary pointer (must be non-empty)
$DISKMINT_BIN_DIR/disk_main --version 2>/dev/null || echo "chemistry binary not found"
radmc3d info                                               # RADMC-3D binary
which optool && optool --version                          # optool (optional but recommended)
gfortran --version | head -1                              # gfortran (need 10+)
```

---

## 0. Python Environment (conda)

**What it is:** A dedicated conda environment keeps DiskMINT dependencies isolated and reproducible.

**Install conda** (if not already available):
- Miniconda (recommended): https://docs.anaconda.com/miniconda/
- Anaconda (full distribution): https://docs.anaconda.com/anaconda/install/

**Create and activate the environment** (choose any name you like; `diskmint_stable` is used in the DiskMINT documentation as an example):

```bash
conda create -n <your_env_name> python=3.11
conda activate <your_env_name>
```

**Verify:**
```bash
echo $CONDA_DEFAULT_ENV   # should print your environment name (not empty, not 'base')
```

All subsequent installation steps assume a dedicated conda environment is active. If `import diskmint` fails, the most common cause is that the base environment or the wrong environment is active — run `conda activate <your_env_name>` and try again.

---

## 1. DiskMINT

**What it is:** Python3-Fortran thermochemical disk modeling package.

**Install (v1.6.0+):**
```bash
git clone https://github.com/DingshanDeng/DiskMINT.git Yourpath/DiskMINT
cd Yourpath/DiskMINT
make install          # compiles Fortran chemistry + pip install -e .
```

**Targets:**
| Command | Effect |
|---|---|
| `make install` | Full install: Fortran compile + Python editable install |
| `make install_python` | Python only (no Fortran recompile) |
| `make chemistry` | Fortran chemistry only |
| `make link_bin` | Symlink `disk_main`, `disk_extract` to `~/.local/bin` |
| `make clean` | Remove compiled artifacts |
| `make uninstall` | Full removal |

**Verify:**
```bash
python -c "import diskmint.model; print('DiskMINT OK')"
```

---

## 2. Fortran Chemistry Network

**What it is:** Reduced chemical network (Ruaud, Gorti & Hollenbach 2022) for $\mathrm{C^{18}O}$ modeling. Compiled separately from the Python package.

**Requires:** `gfortran` version **10 or higher**.

**Check gfortran version:**
```bash
gfortran --version | head -1
# Output should show version 10, 11, 12, 13, or higher
```

**If version < 10:**
- macOS: `brew install gcc` (installs gfortran-13 or similar alongside)
- Linux: `sudo apt install gfortran-10` or `sudo dnf install gcc-gfortran`
- ARM Mac: **do not use the native Darwin gfortran stub** — install via `brew install gcc` and use the versioned binary (e.g. `gfortran-13`)

**Patch the Makefile if needed** (when default `gfortran` is not version 10+):
```bash
# In Yourpath/DiskMINT/chemistry/src/Makefile, change:
FC = gfortran
# to the versioned binary, e.g.:
FC = gfortran-13
```

**Compile:**
```bash
cd Yourpath/DiskMINT/chemistry/src
# If upgrading: remove old build artifacts first
rm -f *.o *.mod
make
```

**Set `DISKMINT_BIN_DIR`** (tells DiskMINT where the compiled binary lives):
```bash
export DISKMINT_BIN_DIR="Yourpath/DiskMINT/chemistry/bin"
# Add to ~/.zshrc or ~/.bashrc to make permanent:
echo 'export DISKMINT_BIN_DIR="Yourpath/DiskMINT/chemistry/bin"' >> ~/.zshrc
source ~/.zshrc
```

**Verify:**
```bash
echo $DISKMINT_BIN_DIR          # should be non-empty
ls $DISKMINT_BIN_DIR/disk_main  # binary should exist
```

---

## 3. RADMC-3D

**What it is:** Dust and gas radiative transfer code. Required for all DiskMINT model runs.
**Tested version:** v2.0

**Install:**
```bash
git clone https://github.com/dullemond/radmc3d-2.0.git
cd radmc3d-2.0/src
make
```

**Add to PATH:**
```bash
export PATH="$PATH:/your/path/to/radmc3d-2.0/bin"
# Add to ~/.zshrc or ~/.bashrc to make permanent
echo 'export PATH="$PATH:/your/path/to/radmc3d-2.0/bin"' >> ~/.zshrc
source ~/.zshrc
```

**Platform notes:**
- Linux (Intel): standard `make` works
- Intel Mac: standard `make` works
- ARM Mac: may need `FFLAGS` adjustments; install `gcc` via Homebrew first

**Verify:**
```bash
radmc3d info
```

**Full documentation:** https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/index.html

---

## 4. optool

**What it is:** Computes dust opacity files from grain composition and size distribution. Required to generate custom `RADMC-3D`-format opacity files.
**Not needed** if you are using pre-computed opacity files from the examples.

**Install:**
```bash
git clone https://github.com/cdominik/optool.git
cd optool
make
export PATH="$PATH:/your/path/to/optool"
echo 'export PATH="$PATH:/your/path/to/optool"' >> ~/.zshrc
source ~/.zshrc
```

**Verify:**
```bash
optool --version
```

**DiskMINT wrapper:** `examples/example_utils/diskmint_utils.py` provides `wrapper_optool_opac()` to call optool with DIANA-standard grain parameters and output directly to RADMC-3D format. See {doc}`../User Guide/user_guide_index` for usage.

---

## 5. Line Radiative Transfer (choose one)

**What it is:** Converts DiskMINT density/temperature output into synthetic spectral line images.

### RADMC-3D (recommended, LTE only)
Already installed in step 3. No additional setup needed.
Handles LTE line radiative transfer natively.

### LIME (non-LTE)
Use when non-LTE excitation is needed (e.g., low-density regions, envelope).
Tested version: v1.9.5

```bash
git clone https://github.com/lime-rt/lime.git
# Follow LIME build instructions at https://github.com/lime-rt/lime
```

---

## Platform Quick Reference

| Platform | gfortran source | RADMC-3D notes |
|---|---|---|
| Linux (Intel) | `apt install gfortran` or system default | Standard `make` |
| Intel Mac | `brew install gcc` | Standard `make` |
| ARM Mac | `brew install gcc` → use `gfortran-13` | May need `FFLAGS`; avoid Darwin stub |

---

## Full Verification Block

```bash
python -c "import diskmint.model; print('DiskMINT OK')"
echo "DISKMINT_BIN_DIR = $DISKMINT_BIN_DIR"
ls $DISKMINT_BIN_DIR/disk_main
radmc3d info
which optool && optool --version
gfortran --version | head -1
```

All lines should return without errors before running a model.
