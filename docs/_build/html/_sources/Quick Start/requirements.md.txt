# Requirements

`DiskMINT` is a Python-Fortran hybrid package that coordinates several external tools to build a complete disk model.
This page explains what each tool does and how to install it.

---

## RADMC-3D (Required)

`DiskMINT` relies on [`RADMC-3D` v2.0](https://github.com/dullemond/radmc3d-2.0) to compute radiative transfer and generate the dust and gas thermal structure of the disk.
Without `RADMC-3D`, the model cannot run.

**Install steps:**

1. Clone the repository:
   ```bash
   git clone https://github.com/dullemond/radmc3d-2.0.git
   ```
2. Go into the source directory and compile:
   ```bash
   cd radmc3d-2.0/src
   make
   ```
3. Add the `radmc3d` binary to your `PATH`. For example, add the following line to your `~/.bashrc` or `~/.zshrc`:
   ```bash
   export PATH="$PATH:/your/path/to/radmc3d-2.0/bin"
   ```
4. Verify the installation:
   ```bash
   radmc3d info
   ```

For full documentation and troubleshooting, see the [RADMC-3D manual](https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/index.html).

---

## Dust Opacity Files (Required)

`DiskMINT` models need dust opacity files in `RADMC-3D` format — one for each grain size population — in order to compute the Spectral Energy Distribution (SED) and the dust thermal structure.
You need to bring your own opacity files for your target.

We recommend generating them with one of these two tools:

### Option 1: `optool`

[`optool`](https://github.com/cdominik/optool) computes dust opacities from grain composition and size distribution, supporting a wide range of materials and mixing rules.

**Install steps:**

1. Clone and compile:
   ```bash
   git clone https://github.com/cdominik/optool.git
   cd optool
   make
   ```
2. Add to your `PATH`:
   ```bash
   export PATH="$PATH:/your/path/to/optool"
   ```
3. Verify:
   ```bash
   optool --version
   ```

### Option 2: `dsharp_opac` (Python package)

[`dsharp_opac`](https://github.com/birnstiel/dsharp_opac) provides pre-computed DSHARP dust opacities as a Python package, which is a convenient choice if you want to use the same opacity model as the DSHARP survey.

**Install:**
```bash
pip install dsharp-opac
```

Example opacity files generated with both tools are provided in the `/examples/` directory of `DiskMINT` as a starting point.

---

## Line Radiative Transfer Code (Required for Synthetic Line Emission)

`DiskMINT` outputs the density and temperature structure of the disk. To produce synthetic line emission images (e.g., CO, $\mathrm{C^{18}O}$), you need a separate line radiative transfer code.
Two options are tested and supported:

### Option 1: `RADMC-3D` (LTE only)

`RADMC-3D` includes a built-in line radiative transfer module that supports Local Thermal Equilibrium (LTE). Since you already need `RADMC-3D` for the dust continuum step, this requires no additional installation.
This is the recommended approach for most users.

### Option 2: `LIME` (LTE and non-LTE)

[`LIME` v1.9.5](https://github.com/lime-rt/lime) supports both LTE and non-LTE line radiative transfer.
It was used for the results published in [Deng et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023arXiv230702657D/abstract).
Use `LIME` if you need non-LTE excitation (e.g., for low-density envelope regions).

Installation instructions are on the [LIME GitHub page](https://github.com/lime-rt/lime).

---

## DiskMINT Utilities for External Tools

`DiskMINT` ships with utility functions that wrap these external tools so you do not need to call them manually.
These are provided in `examples/example_utils/diskmint_utils.py`:

- **`wrapper_optool_opac()`** — generates dust opacity files for multiple grain size bins by calling `optool` with the correct arguments and moving the output into the `RADMC-3D`-compatible format expected by `DiskMINT`. The default grain composition follows the DIANA standard opacity (pyroxene + carbon, with 25% porosity).

More information about the model input/output and how to set up the models are described in the {doc}`../User Guide/user_guide_index`.

---

## Verify Everything Is Working

After installing all tools, you can run the following to confirm your environment is set up correctly:

```bash
python -c "import diskmint.model; print('DiskMINT OK')"
radmc3d info
optool --version        # if using optool
echo $DISKMINT_BIN_DIR  # should point to the compiled DiskMINT chemistry binary
```

---

For more details on how each tool fits into a full model run, see the {doc}`../User Guide/user_guide_index`.
