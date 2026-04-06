# Error Reference

Known errors, root causes, and fixes — grouped by pipeline stage. Use this alongside the {doc}`install_reference` for installation-time errors and the escalation template for unresolved issues.

---

## Installation Errors

### `Cannot import diskmint`

**Symptom:**
```
ModuleNotFoundError: No module named 'diskmint'
```

**Cause:** DiskMINT is not installed or the conda environment is not active.

**Fix:**
```bash
conda activate diskmint_stable      # or your environment name
cd Yourpath/DiskMINT
make install
python -c "import diskmint.model; print('OK')"
```

---

### gfortran compile error: `Error: Rank mismatch`

**Symptom:** `make` fails in `chemistry/src/` with rank mismatch or type errors.

**Cause:** gfortran version < 10. Older gfortran does not support the Fortran 2008 features used in the chemistry network.

**Fix:**
1. Check version: `gfortran --version`
2. Install gfortran 10+:
   - Linux: `sudo apt install gfortran-12` (or via your distro's package manager)
   - Intel Mac: `brew install gcc`
   - ARM Mac: `brew install gcc` (do **not** use the native Darwin stub)
3. Edit `chemistry/src/Makefile`: change `FC = gfortran` to `FC = gfortran-12` (or your version)
4. Clean and recompile:
   ```bash
   cd chemistry/src
   rm -f *.o *.mod
   make
   ```

---

### `DISKMINT_BIN_DIR` not set / binary not found

**Symptom:**
```
[DiskMINT Warning] Binary directory not found at: ...
Please set: export DISKMINT_BIN_DIR='/path/to/diskmint/chemistry/bin'
```

**Cause:** The environment variable is not set, or points to a path that does not contain `disk_main`.

**Fix:**
```bash
export DISKMINT_BIN_DIR="Yourpath/DiskMINT/chemistry/bin"
ls $DISKMINT_BIN_DIR/disk_main     # verify the binary exists
```

Add to `~/.zshrc` or `~/.bashrc` to make permanent. See {doc}`../Quick Start/Installation` for details.

---

## RADMC-3D Errors

### `dust_temperature.dat` not produced

**Symptom:** RADMC-3D runs but `dust_temperature.dat` is never written; DiskMINT retries 3 times and prints a warning.

**Causes and fixes:**

| Cause | Fix |
|-------|-----|
| `radmc3d` not on PATH | `which radmc3d` — if missing, add `~/.local/bin` to PATH or reinstall RADMC-3D |
| `dust_density.inp` is all zeros | Check `mdiskd` and `ratio_g2d_global` in CSV — masses must be > 0 |
| Grid too coarse for `rin` | Decrease `rin` or increase `nr` so the inner sublimation zone is resolved |
| Wrong `dustopac.inp` | Check that `dustopacname_1` in CSV matches the actual `dustkappa_{name}_{i}.inp` files in `data_dir` |

---

### RADMC-3D crashes with `Segmentation fault`

**Causes and fixes:**

| Cause | Fix |
|-------|-----|
| Too many photon packages for available memory | Reduce `nphot` to `1e6`; increase RAM on HPC |
| `nthreads` exceeds physical cores | Set `nthreads` to the actual physical core count |
| Corrupt `amr_grid.inp` | Delete `data_dir/*.inp` and re-run `problem_setup` |

---

### `radmc3d mctherm` hangs indefinitely

**Cause:** Likely a deadlock in OpenMP with `nthreads > 1` on some macOS configurations.

**Fix:** Set `nthreads = 1` for local testing; use a Linux HPC for production runs with many threads.

---

## VHSE Convergence Errors

### `delta_rho/rho starts to increase` — divergence detected

**Symptom:** Print output shows mean and max differences growing between iterations.

**Cause:** The grid is too coarse or the model is in an unphysical region of parameter space.

**Fixes:**
- Increase `nr` (e.g., 30 → 50) and `ntheta` (e.g., 100 → 150)
- Increase `nphot` (e.g., `1e6` → `1e7`)
- Check that `rin` is at or inside the dust sublimation radius (typically `rin ≤ Rsub ≈ Rstar × (Tstar/1500 K)^2`)
- Reduce `n_vhse_loop` to 5–8 to stop early with a usable (if unconverged) result

---

### `new density identical to previous iteration` (redundant iteration)

**Symptom:** Many lines of `something wrong in this iteration, repeat this iteration`.

**Cause:** The solver is stuck — usually a numerical issue in an extreme parameter combination (very low or very high disk mass).

**Fix:** Check `mdiskd` and `ratio_g2d_global` are in physically reasonable ranges. A very optically thick disk (`mdiskd > 0.1 M☉`) or very thin disk (`mdiskd < 1e-6 M☉`) may need `nphot` tuning.

---

## Chemistry Network Errors

### `disk_main: command not found`

**Cause:** The chemistry binary was not compiled or `DISKMINT_BIN_DIR` is wrong.

**Fix:** See installation error above. After compiling, verify:
```bash
$DISKMINT_BIN_DIR/disk_main --version 2>/dev/null || echo "not found"
```

---

### Chemistry produces all-zero or NaN abundances

**Symptom:** `COendgrid-*.chem` contains only zeros or `nan`.

**Causes and fixes:**

| Cause | Fix |
|-------|-----|
| `G0Hab_set = 0` | UV field is required; set `G0Hab_set` to a physically motivated value (see {doc}`parameters_reference`) |
| `COinitgrid-*.dat` empty | Check that chemistry_setup ran before runchemistry; look for `COinitgrid` files in `data_dir` |
| Grid dimensions mismatch | Verify `nr_cyl_LIME` and `nz_cyl_LIME` are consistent with your grid |

---

### Chemistry network exits immediately

**Symptom:** `disk_main` starts and exits in < 1 second with no output.

**Cause:** Usually a missing input file (initial abundance grid) or a permission error on the binary.

**Fix:**
```bash
ls -la $DISKMINT_BIN_DIR/disk_main    # check executable bit
ls data/COinitgrid-*.dat              # check input file exists
```

---

## File and Path Errors

### `FileNotFoundError: model_parameters_*.csv`

**Cause:** The parameter CSV path is wrong.

**Fix:** Pass the full path and filename separately:
```python
para.read_parameters_from_csv(
    filename="model_parameters_1p0Msolar.csv",
    directory="/path/to/project/",
    extension=""
)
```

---

### `No files found for pattern: dustkappa_*.inp`

**Cause:** The `dustopacname_1` in the CSV does not match the actual filenames in `data_dir`.

**Fix:** List the opacity files and reconcile:
```bash
ls data/dustkappa_*.inp
# The name between 'dustkappa_' and '_N.inp' must match dustopacname_1
```

---

## Escalation

If an error is not listed here or cannot be resolved with the above steps, use the support escalation procedure:

1. Collect environment info and logs (see {doc}`../AI Features/install_reference`)
2. Check open issues: https://github.com/DingshanDeng/DiskMINT/issues
3. Draft a support email using the template in the DiskMINT-Nursery escalation reference, or open a new issue with the collected information
