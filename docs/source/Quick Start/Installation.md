# Installation

**Quicker Installation (as of v1.6.0+)** 

1. Download the code use the `git clone` or download this repo as zip from the latest [releases](https://github.com/DingshanDeng/DiskMINT/releases) into your local directory `Yourpath/DiskMINT/`.
2. Open terminal, go inside your path `cd Yourpath/DiskMINT/`, type `make install`. This should install both the `Python` and `Fortran` modules into your machine. 
3. Start Using! 

**Manual Installation (all versions)**

1. Download the code to your local directory `Yourpath/DiskMINT/`: clone this repo using `git clone` or other methods to download the part of the code you would like to use. We kindly note that you need to put the package in a directory that you have permission to read, write and execute.
2. To use the `Python3` module (`Yourpath/DiskMINT/diskmint/src/`): add the path to the package to your `Python3` path: you can do this by adding the `Yourpath/DiskMINT/src/` to your `Python3` environment `$PATH`, or write `sys.path.append('Yourpath/DiskMINT/src/')` in the `Python3` script when you want to call the package. 
3. To use the chemical network (`Yourpath/DiskMINT/chemistry/`) we provided: go to `Yourpath/DiskMINT/chemistry/src/`. Edit the `Makefile` in the `src` folder, and change the ` FC = gfortran` to your version of `gfortran`. **Before `make`**, (a) (especially when you need to upgrade the chemical network) If `.o` and `.mod` files already exist in the `/src`, remove them by `rm -i *.o` and `rm -i *.mod`; (b) Make sure that your version of `gfortran` is 10 or higher (you can use `gfortran -v` to check). Note that in some legacy Linux system, the default `gfortran` installed is not the latest version, please check [`gfortran`](https://gcc.gnu.org/fortran) website for newer versions. Then just use command `make` to compile the `Fortran` chemical network.
4. Start Using!

---

## Setting `DISKMINT_BIN_DIR`

After compiling the Fortran chemistry network (step 3 above), you need to tell DiskMINT where to find the compiled binaries by setting the `DISKMINT_BIN_DIR` environment variable:

```bash
export DISKMINT_BIN_DIR="Yourpath/DiskMINT/chemistry/bin"
```

To make this permanent, add it to your shell profile:

```bash
# For zsh (default on macOS):
echo 'export DISKMINT_BIN_DIR="Yourpath/DiskMINT/chemistry/bin"' >> ~/.zshrc
source ~/.zshrc

# For bash:
echo 'export DISKMINT_BIN_DIR="Yourpath/DiskMINT/chemistry/bin"' >> ~/.bashrc
source ~/.bashrc
```

**Verify:**

```bash
echo $DISKMINT_BIN_DIR          # should print the path you set
ls $DISKMINT_BIN_DIR/disk_main  # compiled binary should exist
```

If `DISKMINT_BIN_DIR` is not set, DiskMINT will attempt to locate the binaries relative to the package installation directory. This fallback works for most cases after `make install`, but setting the variable explicitly is recommended to avoid path issues.

> **Note:** If you used `make install` (Quicker Installation), the `Makefile` also runs `make link_bin`, which symlinks `disk_main` and `disk_extract` into `~/.local/bin`. You still need to set `DISKMINT_BIN_DIR` pointing to `Yourpath/DiskMINT/chemistry/bin` for the Python package to locate the binaries correctly.

---

## Requirements

*The code is built on a `Linux` machine with `Intel` CPU and tested on both `Intel-Mac` and `ARM-Mac`, and we note that to run `gfortran` on your ARM Mac, be sure to install `gcc` and not use the native `Darwin Mac` Version.*

*The chemical network is built and tested on `gfortran-10+`.*

**Python dependencies.** The `Python3` module of `DiskMINT` relies on several auxiliary but widely used `Python3` packages:

```
    numpy, scipy: used for array operations and calculations

    # Python 3 standard library
    os, sys, copy, shutil: used for handling files and paths
    subprocess (for chemical network): used for executing Fortran code 
    if you want to use the chemical network we provided
```

*The code was originally built with `Python 3.8`, and later tested on `Python 3.9-3.13`.*

[**!IMPORTANT**] Compatibility Note: If you are using a version of DiskMINT prior to v1.5.0, you must manually install radmc3dPy. Please note that `radmc3dPy` has its own dependencies. For detailed installation instructions, please refer to the official RADMC-3D documentation, and here is a link to [their website](https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/index.html).

For external tool requirements (RADMC-3D, dust opacities, line radiative transfer code) and their installation instructions, see {doc}`requirements`.


