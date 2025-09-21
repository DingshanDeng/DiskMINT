# Quick Start

## INSTALLATION

1. Download the code to your local directory `Yourpath/DiskMINT/`: clone this repo using `git clone` or other methods to download the part of the code you would like to use. We kindly note that you need to put the package in a directory that you have permission to read, write and execute.
2. To use the `Python3` module (`Yourpath/DiskMINT/diskmint/src/`): add the path to the package to your `Python3` path: you can do this by adding the `Yourpath/DiskMINT/src/` to your `Python3` environment `$PATH`, or write `sys.path.append('Yourpath/DiskMINT/src/')` in the `Python3` script when you want to call the package. 
3. To use the chemical network (`Yourpath/DiskMINT/chemistry/`) we provided: go to `Yourpath/DiskMINT/chemistry/src/`. Edit the `Makefile` in the `src` folder, and change the ` FC = gfortran` to your version of `gfortran`. **Before `make`**, (a) (especially when you need to upgrade the chemical network) If `.o` and `.mod` files already exist in the `/src`, remove them by `rm -i *.o` and `rm -i *.mod`; (b) Make sure that your version of `gfortran` is 10 or higher (you can use `gfortran -v` to check). Note that in some legacy Linux system, the default `gfortran` installed is not the latest version, please check [`gfortran`](https://gcc.gnu.org/fortran) website for newer versions. Then just use command `make` to compile the `Fortran` chemical network.

*The code is built on a `Linux` machine with `Intel` CPU and tested on both `Intel-Mac` and `ARM-Mac`, and we note that to run `gfortran` on your ARM Mac, be sure to install `gcc` and not use the native `Darwin Mac` Version.*

## Try the examples

There are a few examples we provided to carry out disk modeling with DiskMINT in `Yourpath/DiskMINT/examples/` directory. For example, you can start running the example of the disk around RU Lup by following the steps below

1. Go to the `Yourpath/DiskMINT/examples/example_RULup/` directory.
2. Open the `Python3` script `example_model_RULup.py`. This script calls the functions from `diskmint` module to build your disk model, calls `RADMC-3D` to do radiative transfer and calculate the thermal structure, and also calls the chemical network we provide to calculate the molecular abundances.
3. Edit the `example_model_RULup.py` following the comments in the script. (a) Edit `package_position = "Yourpath/DiskMINT/src/"` to be the position where you put the module, and you can ignore this line (safe to comment out or delete) if you already add the `Python3` module to your `$PATH`. This is for enabling `import diskmint` in `Python3`; (b) Edit a list of directories and names, including where you want to run the code (`working_dir`), where you want to save the files (`save_dir`), what is the name to be saved for this model (`name_of_this_model`), name of the parameter file (`file_parameters`), and where is the chemical network (`chem_code_dir`); (c) Set up a list of the options of the model, including whether you want to read in your premade dust opacity info (`bool_MakeDustKappa`), to compute the SED (`bool_SED`), to solve the vertical hydrostatic equilibrium (`bool_VHSE`), to run the chemical network (`bool_chemistry`), and to save a copy of the model to your `save_dir` (`bool_savemodel`).
4. Then, you should be able to run the example model of RU Lup using `python3 example_model_RULup.py`.
5. The final output files are the density and thermal distribution in `RADMC-3D` format that can be read by `radmc3dPy` in the `data` folder. Also, if you are using the chemical network, there will be the `COinitgrid-*.dat` (initial grid for CO abundances), `COinitgrid-GSinit_*.dat` (initial grid for CO abundances with pre-set Gaussian vertical distribution before solving VHSE) and `COendgrid-*.chem` (the CO abundances after chemical network). The columns in these `COendgrids-*.chem` are {r[cm], z[cm], log10(abundanceH2)[relative to H nuclei], log10(abundanceC18O)[relative to H nuclei]}.
6. Then, you can use the line radiative transfer code to do the following line radiative transfer. We provide a simple wrapper for running `LIME` as an example in the `example_RULup/LIME_example_RULup/` folder.

## Build Your First DiskMINT Model

For using the model on another target, or if you want to play with it a bit, there are a few things that can be changed.

1. **Edit Parameter File.** Open (and edit) the file `example_RULup_parameters.dat` that sets all the parameters for this model. The data file structure follows `.csv` format (comma `,` as separator) and ignores the line starting with `#`. The columns are with the sequence `name(str), value(str/int/float), units(str), valuetype, description`. The parameters set in this `file_parameters` would be read in the `example_model_RULup.py` with the function `para.read_parameters_from_csv(filename=file_parameters)`, and the parameters can then be changed with `para.edit_parameter(para_name, new_value=new_value)`. So, some parameters already set in the `file_parameters`, such as the model name, are edited in the `Python3` script with `para.edit_parameter("chemical_save_name", new_value = name_of_this_model)`.

2. **Use New Dust Opacity - Part A.** To change the dust opacity (you need to make your own `RADMC-3D` format opacity files using other software such as `optool` or `dsharp_opac`), you need to edit the name of the opacity file (`dustopacname_1`) in the `file_parameters`. Currently, the code supports two different dust compositions, `dustopacname_1` and `dustopacname_2`, and the option of using only one `dustopacname_1` is well-tested. For each dust composition, it needs to be separated into different dust size bins, with the name sequence of `'dustkappa_'+dustopacname_1+'_'+str(i_dust)+'.inp`, where `i_dust` is an integer from `i_dust = 0` to the number of the total size bins `i_dust = (nr_dust_1 - 1)`, and `0` represents the smallest one while the `nr_dust_1 - 1` for the largest size. The `dustkappa_xxx.inp` files need to be put in the `data_dir` that was set in the `model.py` file.

3. **Use New Dust Opacity - Part B.** Besides the `dustkappa_xxx.inp` files, there are also a few other files related to the dust opacity information that need to be provided, including:

   1. `aave.inp` (1-d array): average dust radius of dust species in the unit of `cm`.

   2. `ndsd.inp` (1-d array): `nd(a) * np.pi * a**2`. Where `nd` is the number density of the dust with size `a`.

   3. `frac.inp` (1-d array): mass fractions of dust species -- the mass of the dust with this composition and size divided by the total mass -- the sum needs to be `1`.

   4. `frac_nb.inp` (1-d array): the number fractions of dust species -- number density of the dust with this composition and size divided by the total number density -- the sum needs to be `1`.

      You can follow the example files in the `data` folder to make your own files for the dust composition you used, and for example, if you already set up the `aave_array = np.array([1e-5, 1e-4, 1e-3])`, then you can save it with `np.savetxt('aave.inp', aave_array)`.
      We also provide the function in the code for *making these four `.inp` files* that for the dust opacity (currently only support single composition with `dustopacname_1`) with a simple power law distribution ($nd(a) \sim a^{-\mathrm{pla}}$ with $a_{\mathrm{min}} < a < a_{\mathrm{max}}$) and a uniformly distributed size in log-space (`ai = np.logspace(np.log10(para.amin_1), np.log10(para.amax_1), para.nr_dust_1 + 1)`).

      To use this, you need to turn on the `bool_MakeDustKappa = True` and correctly change the parameters in the `file_parameters` (Section `4. Dust Kappa Setup`), including `dustopacname_1`, `nr_dust_1`, `amin_1`, `amax_1`, `amin_all`, `amax_all`, `pla_dustsize`, and `rhobulk`.

4. **Star (+UV) Spectrum.** In addition to the dust opacity information, the stellar spectrum (including the UV spectrum) is also required. It needs to be in the `RADMC-3D` format, following the example of `RULupSpecModel_wuv.inp` in the `data` folder. The UV spectrum is required to estimate the UV field close to the star, and it is important for CO chemistry. If you have any UV observations on the target, you can estimate the UV spectrum by fitting a blackbody profile. Or if you do not want to do the fitting or not having the UV observations, scaling the UV spectrum on TW Hya to your target should also give a good estimation. *If you are not using the chemical network part of the code, then only providing the stellar spectrum would still be fine.*

With all of these files and parameters being correctly set up, you should be able to start running your first model.
Also, remember to change the script you have (not provided here) for doing line radiative transfer and making the synthetic image for the new target or new observation to compare with.
