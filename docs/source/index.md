# DiskMINT

**Disk Model for INdividual Targets**

<p align="center">
  <img src="_static/assets/images/card-software-transparent.png" alt="" style="width:38.2%; height:auto;">
</p>

**Contributors:**

Dingshan Deng (dingshandeng@arizona.edu),
The University of Arizona

Maxime Ruaud,
SETI Institute

Uma Gorti,
SETI Institute

Ilaria Pascucci,
The University of Arizona

---

Acknowledgement

`DiskMINT` now includes standalone functions for handling RADMC-3D input and output files. These functions are based on the original `radmc3dPy` scripts developed by Dr. Attila Juhasz and maintained by Prof. Dr. Cornelis P. Dullemond.

We also thank the anonymous referee, science editors, and data editors at AAS Journals for their helpful suggestions and comments, which significantly improved this work. 
Finally, we express our gratitude to the many people in our community whose feedback and support have been instrumental to the growth of this project.

## About DiskMINT

`DiskMINT` (Disk Model for INdividual Targets) is a tool for modeling individual disks and deriving more robust disk mass estimates.
`DiskMINT` is a `Python3`-`Fortran` code built on [`RADMC-3D` v2.0](https://github.com/dullemond/radmc3d-2.0) for the continuum (and gas line) radiative transfer, and it includes a reduced chemical network suggested in [Ruaud, Gorti & Hollenbach (2022)](https://ui.adsabs.harvard.edu/abs/2022ApJ...925...49R/abstract) to determine the $\mathrm{C^{18}O}$ emission.
For more information, see [Deng et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023arXiv230702657D/abstract) and [Deng et al. (2025)](https://ui.adsabs.harvard.edu/abs/2025ApJ...995...98D).

For more details about `DiskMINT`, see [Documentation](https://diskmint.readthedocs.io/en/latest/).

### The Code

The code contains two parts:

(A) a `Python3` module (located in `DiskMINT/diskmint/src/`) that is capable of generating a self-consistent 2D disk structure that satisfies VHSE (Vertical Hydrostatic Equilibrium).

(B) a `Fortran` code (located in `DiskMINT/chemistry/`) of the reduced chemical network focusing on modeling $\mathrm{C^{18}O}$ in the disk, and it contains the main chemical processes that are necessary for the $\mathrm{C^{18}O}$ modeling: a. the isotope-selective photodissociation; b. the grain-surface chemistry where the CO converting to $\mathrm{CO_2}$ ice is the main reaction.

`DiskMINT` is still a developing tool, and we aim to build it as an easy-to-update and easy-to-use tool. 
Feel free to put any questions, suggestions, and/or comments on the `GitHub` page or directly email the author (dingshandeng@arizona.edu).

`DiskMINT` is open source and released under the [MIT License](https://opensource.org/licenses/MIT). You are free to use, modify, and distribute this codebase in accordance with the terms and conditions of the license.

---

## Quick START

1. Download the code use the `git clone` or download this repo as zip into your local directory `Yourpath/DiskMINT/`.
2. Open terminal, go inside your path `cd Yourpath/DiskMINT/`, type `make install`. This should install both the `Python` and `Fortran` modules into your machine. 
3. Start Using! 

For more information, please check our {doc}`Quick Start Guide </quick_start/quick_start_index>`.

---

## Citations
If you use `DiskMINT` as part of your research, please cite 

```
@ARTICLE{2025ApJ...995...98D,
       author = {{Deng}, Dingshan and {Gorti}, Uma and {Pascucci}, Ilaria and {Ruaud}, Maxime},
        title = "{DiskMINT: Self-consistent Thermochemical Disk Models with Radially Varying Gas and Dust{\textemdash}Application to the Massive, CO-Rich Disk of IM Lup}",
      journal = {\apj},
     keywords = {Protoplanetary disks, Astrochemistry, Chemical abundances, CO line emission, Planet formation, 1300, 75, 224, 262, 1241, Earth and Planetary Astrophysics},
         year = 2025,
        month = dec,
       volume = {995},
       number = {1},
          eid = {98},
        pages = {98},
          doi = {10.3847/1538-4357/ae0e66},
archivePrefix = {arXiv},
       eprint = {2509.15487},
 primaryClass = {astro-ph.EP},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2025ApJ...995...98D},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

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
Quick Start <quick_start/quick_start_index>
User Guide <user_guide/user_guide_index>
Examples <examples/example_index>
ML & AI Assistant <ai_features/ai_ref_index>
API Reference <api>
```
