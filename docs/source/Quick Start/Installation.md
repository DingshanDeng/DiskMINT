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

## Requirements:

*The code is built on a `Linux` machine with `Intel` CPU and tested on both `Intel-Mac` and `ARM-Mac`, and we note that to run `gfortran` on your ARM Mac, be sure to install `gcc` and not use the native `Darwin Mac` Version.*

*The chemical network is built and tested on `gfortran-10+`.*

**RADMC-3D.** `DiskMINT` relies on `RADMC-3D` (currently built on [`RADMC-3D` v2.0](https://github.com/dullemond/radmc3d-2.0)) to do radiative transfer and generate the thermal distributions.
So please make sure the `RADMC-3D` is properly installed according to the instructions on their website [their website](https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/index.html).

**Python Module.** The `Python3` module of `DiskMINT` relies on several auxiliary but widely used `Python3` packages:

```
    numpy, scipy: used for array operations and calculations

    # Python 3 standard library
    os, sys, copy, shutil: used for handling files and paths
    subprocess (for chemical network): used for executing Fortran code 
    if you want to use the chemical network we provided
```

*The code was originally built with `Python 3.8`, and later tested on `Python 3.9-3.13`.*

[**!IMPORTANT**] Compatibility Note: If you are using a version of DiskMINT prior to v1.5.0, you must manually install radmc3dPy. Please note that `radmc3dPy` has its own dependencies. For detailed installation instructions, please refer to the official RADMC-3D documentation, and here is a link to [their website](https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/index.html).

### Other dependences:

**Line Radiative Transfer Code.** The `DiskMINT` module will generate the thermal and density distributions of the disk, which can then be used to do line radiative transfer and make synthetic images.
We do not provide the code to do the line radiative transfer and synthetic imaging (*But do not worry, examples are provided*). 
To get the synthetic image, please use `LIME` -- supports LTE (Local Thermal Equilibrium) and non-LTE -- or the `RADMC-3D` -- only supports LTE -- line radiative transfer module. 
`LIME (v1.9.5)` and `RADMC-3D (v2.0)` are tested on `DiskMINT` model.
For the example results as presented in [Deng et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023arXiv230702657D/abstract), we used [`LIME` v1.9.5](https://github.com/lime-rt/lime) to do non-LTE line radiative transfer.

We also note that to compare with the observation, a proper **image convolution** based on the beam of observation is necessary, and including the continuum emission from the dust in the synthetic image would be helpful in many cases (the continuum can then be subtracted with the same method that is adopted for the observation dataset, or directly compare the synthetic image with continuum to the observation dataset that also contains the continuum).

**Dust Opacities.** Because one major step in the model is to match the Spectral Energy Distribution (SED) of the disk, which requires different dust opacities for each disk.
Therefore, to properly use the model, you need to bring your own dust opacity files (in `RADMC-3D` format), and we also provide a few examples here in the `/examples/`
We recommend using [`optool`](https://github.com/cdominik/optool) or [`dsharp_opac`](https://github.com/birnstiel/dsharp_opac) to generate the opacities for multiple size distributions that are required to correctly calculate the gas thermal structure from gas-dust thermal equilibrium.


