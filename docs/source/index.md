# DiskMINT

**Disk Model for INdividual Targets**

**Contributors:**

Dingshan Deng (dingshandeng@arizona.edu),
The University of Arizona

Maxime Ruaud,
SETI Institute

Uma Gorti,
SETI Institute

Ilaria Pascucci,
The University of Arizona

## About DiskMINT

`DiskMINT` (Disk Model for INdividual Targets) is a tool for modeling individual disks and deriving more robust disk mass estimates.
`DiskMINT` is a `Python3`-`Fortran` code built on [`RADMC-3D` v2.0](https://github.com/dullemond/radmc3d-2.0) for the continuum (and gas line) radiative transfer, and it includes a reduced chemical network suggested in [Ruaud, Gorti & Hollenbach (2022)](https://ui.adsabs.harvard.edu/abs/2022ApJ...925...49R/abstract) to determine the ${\rm C^{18}O}$ emission.
For more information, see [Deng et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023arXiv230702657D/abstract).

### The Code

The code contains two parts:

(A) a `Python3` module (located in `DiskMINT/diskmint/src/`) that is capable of generating a self-consistent 2D disk structure that satisfies VHSE (Vertical Hydrostatic Equilibrium).

(B) a `Fortran` code (located in `DiskMINT/chemistry/`) of the reduced chemical network focusing on modeling C18O in the disk, and it contains the main chemical processes that are necessary for the C18O modeling: a. the isotopologue-selective photodissociation; b. the grain-surface chemistry where the CO converting to CO2 ice is the main reaction.

`DiskMINT` is still a developing tool, and we aim to build it as an easy-to-update and easy-to-use tool. 
Feel free to put any questions, suggestions, and/or comments on the `GitHub` page or directly email the author (dingshandeng@arizona.edu).

`DiskMINT` is open source and released under the [MIT License](https://opensource.org/licenses/MIT). You are free to use, modify, and distribute this codebase in accordance with the terms and conditions of the license.

---

## REQUIREMENTS:

**RADMC-3D.** `DiskMINT` relies on `RADMC-3D` (currently built on [`RADMC-3D` v2.0](https://github.com/dullemond/radmc3d-2.0)) to do radiative transfer and generate the thermal distributions.
So please make sure the `RADMC-3D` is properly installed.

**Python Module.** The `Python3` module of `DiskMINT` relies on several auxiliary but widely used `Python3` packages:

```
    numpy, scipy: used for array operations and calculations

    # Python 3 standard library
    os, sys, copy, shutil: used for handling files and paths
    subprocess (for chemical network): used for executing Fortran code 
    if you want to use the chemical network we provided
```

> Acknowlegement: In this recent update (v1.5.0), the package [`radmc3dPy`](https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/index.html) is no longer required. New stand-alone functions are added to `DiskMINT` to handle the input/output files for `RADMC-3D`. These functions are based on the original `radmc3dPy` scripts that are developed by Dr. Attila Juhasz and later maintained by Prof. Dr. Cornelis P. Dullemond. If you are still using the old version of `DiskMINT` (any of the pre-release versions before v1.5.0), please make sure to install `radmc3dPy`. We note that there are also a few other packages that are required by `radmc3dPy`, such as `Matplotlib`, so please check [their website](https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/index.html) to properly install the package.

**Line Radiative Transfer Code.** The `DiskMINT` module will generate the thermal and density distributions of the disk, which can then be used to do line radiative transfer and make synthetic images.
We do not provide the code to do the line radiative transfer and synthetic imaging (*But do not worry, examples are provided*). 
To get the synthetic image, please use `LIME` -- supports LTE (Local Thermal Equilibrium) and non-LTE -- or the `RADMC-3D` -- only supports LTE -- line radiative transfer module. 
`LIME (v1.9.5)` and `RADMC-3D (v2.0)` are tested on `DiskMINT` model.
For the example results as presented in [Deng et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023arXiv230702657D/abstract), we used [`LIME` v1.9.5](https://github.com/lime-rt/lime) to do non-LTE line radiative transfer.

We also note that to compare with the observation, a proper **image convolution** based on the beam of observation is necessary, and including the continuum emission from the dust in the synthetic image would be helpful in many cases (the continuum can then be subtracted with the same method that is adopted for the observation dataset, or directly compare the synthetic image with continuum to the observation dataset that also contains the continuum).

**Dust Opacities.** Because one major step in the model is to match the Spectral Energy Distribution (SED) of the disk, which requires different dust opacities for each disk.
Therefore, to properly use the model, you need to bring your own dust opacity files (in `RADMC-3D` format), and we also provide a few examples here in the `/examples/`
We recommend using [`optool`](https://github.com/cdominik/optool) or [`dsharp_opac`](https://github.com/birnstiel/dsharp_opac) to generate the opacities for multiple size distributions that are required to correctly calculate the gas thermal structure from gas-dust thermal equilibrium.

**In development.** We will also release a utility file/functions to handle the input/output of the line radiative transfer and dust opacities in late 2025 together with a install-friendly version and a detailed documentation.

---

## Citations
If you use `DiskMINT` as part of your research, please cite 

```
@ARTICLE{2023ApJ...954..165D,
       author = {{Deng}, Dingshan and {Ruaud}, Maxime and {Gorti}, Uma and {Pascucci}, Ilaria},
        title = "{DiskMINT: A Tool to Estimate Disk Masses with CO Isotopologues}",
      journal = {\apj},
     keywords = {Protoplanetary disks, Astrochemistry, Chemical abundances, CO line emission, Planet formation, 1300, 75, 224, 262, 1241, Astrophysics - Earth and Planetary Astrophysics, Astrophysics - Instrumentation and Methods for Astrophysics, Astrophysics - Solar and Stellar Astrophysics},
         year = 2023,
        month = sep,
       volume = {954},
       number = {2},
          eid = {165},
        pages = {165},
          doi = {10.3847/1538-4357/acdfcc},
archivePrefix = {arXiv},
       eprint = {2307.02657},
 primaryClass = {astro-ph.EP},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2023ApJ...954..165D},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

## Community Guidelines
We welcome contributions, issue reports, and questions about `DiskMINT`! If you encounter a bug or issue, check out the [Issues page](https://github.com/DingshanDeng/DiskMINT/issues) and provide a report with details about the problem and steps to reproduce it. For general support, usage questions and suggestions, you can start a discussion in [Discussions page](https://github.com/DingshanDeng/DiskMINT/discussions), and of course feel free to send emails directly to us. If you want to contribute, feel free to fork the repository and create pull requests here. `DiskMINT` is licensed under MIT license, so feel free to make use of the source code in any part of your own work/software.

```{toctree}
:hidden:
:titlesonly:
:maxdepth: 2

About DiskMINT <self>
Quick Start <quickstart>
Advanced Usage <advanced>
API Reference <api>
