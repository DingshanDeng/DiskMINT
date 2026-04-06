# Running Your First Model

## Try the examples of RU Lup

> For the full example description including all input files and the LIME workflow, see the {doc}`RU Lup example page <../examples/example_RULup>`.

There are a few examples we provided to carry out disk modeling with DiskMINT in `Yourpath/DiskMINT/examples/` directory. For example, you can start running the example of the disk around RU Lup by following the steps below

1. Go to the `Yourpath/DiskMINT/examples/example_RULup/` directory.
2. Open the `Python3` script `example_model_RULup.py`. This script calls the functions from `diskmint` module to build your disk model, calls `RADMC-3D` to do radiative transfer and calculate the thermal structure, and also calls the chemical network we provide to calculate the molecular abundances.
3. Edit the `example_model_RULup.py` following the comments in the script. 
   1. Edit `package_position = "Yourpath/DiskMINT/src/"` to be the position where you put the module, and you can ignore this line (safe to comment out or delete) if you already add the `Python3` module to your `$PATH`. This is for enabling `import diskmint` in `Python3`; 
   2. Edit a list of directories and names, including where you want to run the code (`working_dir`), where you want to save the files (`save_dir`), what is the name to be saved for this model (`name_of_this_model`), name of the parameter file (`file_parameters`), and where is the chemical network (`chem_code_dir`); 
   3. Set up a list of the options of the model, including whether you want to read in your premade dust opacity info (`bool_MakeDustKappa`), to compute the SED (`bool_SED`), to solve the vertical hydrostatic equilibrium (`bool_VHSE`), to run the chemical network (`bool_chemistry`), and to save a copy of the model to your `save_dir` (`bool_savemodel`).
4. Then, you should be able to run the example model of RU Lup using `python3 example_model_RULup.py`. 
   1. We also highly recommend running the code with a log file. For example, using `python3 example_model_RULup.py 2>&1 | tee -a diskmint_example_rulup.log`
5. The final output files are the density and thermal distribution in `RADMC-3D` format that can be read by `diskmint` (or `radmc3dPy` ) in the `data` folder. Also, if you are using the chemical network, there will be the `COinitgrid-*.dat` (initial grid for CO abundances), `COinitgrid-GSinit_*.dat` (initial grid for CO abundances with pre-set Gaussian vertical distribution before solving VHSE) and `COendgrid-*.chem` (the CO abundances after chemical network). The columns in these `COendgrids-*.chem` are {r[cm], z[cm], log10(abundanceH2)[relative to H nuclei], log10(abundanceC18O)[relative to H nuclei]}.
6. Then, you can use the line radiative transfer code to do the following line radiative transfer. We provide a simple wrapper for running `LIME` as an example in the `example_RULup/LIME_example_RULup/` folder.

## Adapting the Example to a New Target

The RU Lup example above is meant to get you to a first successful run quickly. If you want to modify the parameter file, change dust opacities, prepare the stellar and UV spectrum, or understand the chemistry outputs before applying DiskMINT to a new target, use the {doc}`User Guide <../User\ Guide/build_your_own_model>`. More details on what it means for different parameters are also described there.

## Build Your Own Model and Check Other Examples

::::{grid} 2
:gap: 3

:::{grid-item-card} Build Your Own Model
:link: ../User\ Guide/user_guide_index
:link-type: doc

Learn how to build your own DiskMINT model step by step in the User Guide.
:::

:::{grid-item-card} Browse Other Examples
:link: ../examples/example_index
:link-type: doc

Explore additional model examples and ready-to-run pipelines in the Examples section.
:::

::::
