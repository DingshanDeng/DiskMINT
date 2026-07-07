#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author: Dingshan Deng @ University of Arizona
# contact: dingshandeng@arizona.edu
# created: 07/07/2026

"""
Simple example: install an externally derived dust surface density radial
profile into a DiskMINT model (single flat run, no CLI arguments).

Run from this directory:
python -u 0-model_radial_profile_fit_simple.py 2>&1 | tee -a output_radial_profile_fit_simple.log
"""

import os, shutil, glob
import numpy as np
import time

import diskmint.model as model
import diskmint.execute as exe

time0 = time.time()

#
# I. Set up the directories and the names
#
# (0) Working directory: this is the place that keeps all the example files together
working_dir = os.path.dirname(os.path.abspath(__file__))
#
# (1) Data directory: where the RADMC-3D and chemistry model input/output files are kept
data_dir = os.path.join(working_dir, 'data')
#
# (2) Input files directory: dustkappa_* files, "aave.inp", "fracs.inp",
# "fracs_numb.inp", "ndsd.inp", and the stellar spectrum para.fmodel_filename
read_dir = os.path.join(working_dir, 'input')
#
# (3) Saving files directory: where the chemistry and RADMC-3D outputs are saved
save_dir = os.path.join(working_dir, 'output')
#
# (4) The parameter file: shared with 1-model_radial_profile_fit_advanced.py,
# not duplicated here
file_parameters = 'radial_profile_fit_parameters'
#
# (5) Chemistry code directory: name of the bundled chemical network template
chem_code_dir = 'reducedRGH22'

#
# II. The Options of the Model
#
# whether to create new dust kappa files (fracs, aave, etc.) instead of only reading them
bool_MakeDustKappa = True
#
# whether to compute the SED
bool_SED = True
#
# whether to use the VHSE distribution
bool_VHSE = True
# the maximum number of loops to solve the VHSE structure
# (the iteration stops early once it converges, normally <10 iterations)
n_vhse_loop = int(10)
#
# whether to enable dust settling
bool_dust_settling = True
#
# whether to run the chemical network
bool_chemistry = True
#
# whether to save the radmc3d model results in another folder
bool_savemodel = True

#
# III. Read in the parameters and set up the reference radial profile
#
para = model.Parameters()
para.read_parameters_from_csv(filename=file_parameters, directory=working_dir, extension='.csv')

#
# chemical_save_dir must be overridden to an absolute path: the CSV stores it
# as the relative string "output", but this script os.chdir()s into data_dir
# before running, so without this override the outputs would land in
# data/output/ instead of the shared top-level output/
para.edit_parameter("chemical_save_dir", new_value=save_dir)
#
# give this run its own name so it does not collide with the advanced
# script's model_A/model_B outputs under output/
para.edit_parameter("chemical_save_name", new_value="example_radial_profile_fit_simple")

#
# mdiskd, pl_sufdens, pl_tapoff, Rtap, and ratio_g2d_global are already set to
# the values used by this example inside radial_profile_fit_parameters.csv,
# so they are intentionally left untouched here (single source of truth).
# NOTE: once sigmad_ref is installed below, pl_sufdens/pl_tapoff/Rtap only
# describe the power-law fallback DiskMINT would use if no reference table
# were given -- they no longer set the actual dust surface density profile.
#

#
# Load an externally derived dust surface density radial profile and install
# it as the reference used to build the model's density structure. The file
# is a plain two-column ASCII table: [radius_cm, sigma_dust_g_cm-2].
#
sigmad_ref_read = np.loadtxt(os.path.join(working_dir, 'model_inputs', 'sigma_reference_template.dat'))
para.edit_parameter("sigmad_ref", new_value=sigmad_ref_read)
#
# A radially-varying gas-to-dust ratio can be supplied the exact same way, by
# loading model_inputs/ratio_g2d_reference_template.dat with np.loadtxt and
# calling para.edit_parameter("ratio_g2d", new_value=...). This example keeps
# ratio_g2d_global constant (from the CSV) for simplicity; see
# 1-model_radial_profile_fit_advanced.py --mode model_B_radial_profile_fit
# for the radially-varying gas-to-dust ratio worked example.
#

#
# Set up the working directory and change into it
#
if not os.path.exists(data_dir):
    os.makedirs(data_dir)
    print(f"Directory '{data_dir}' created.")
print('current working directory:\n' + os.path.dirname(os.path.abspath(__file__)))
os.chdir(data_dir)
print('current working directory (after changed):\n' + os.path.dirname(os.path.abspath(__file__)))

#
# IV. Copy the input files into the data directory
#
files_input_list = [
    para.fmodel_filename,
    f"dustkappa_{para.dustopacname_1}_*",
    f"dustkappa_{para.dustopacname_2}_*",
    "aave.inp", "fracs.inp", "fracs_numb.inp", "ndsd.inp",
]

for file_pattern in files_input_list:
    source_paths = glob.glob(os.path.join(read_dir, file_pattern))
    for source_path in source_paths:
        if os.path.isfile(source_path):
            destination_path = os.path.join(data_dir, os.path.basename(source_path))
            shutil.copy(source_path, destination_path)
            print(f"Copied {source_path} to {destination_path}")
        else:
            print(f"File not found: {source_path}")

#
# V. Build and run the model
#
print('--------------------')
print('%s' % (para.chemical_save_name))
print('--------------------')

mint = model.Mint(para, file_dir=data_dir)
mint.bool_MakeDustKappa = bool_MakeDustKappa
mint.bool_SED = bool_SED
mint.bool_VHSE = bool_VHSE
mint.n_vhse_loop = n_vhse_loop
mint.bool_dust_settling = bool_dust_settling
mint.bool_chemistry = bool_chemistry
mint.bool_savemodel = bool_savemodel
mint.chem_code_dir = chem_code_dir

exe.runmodel(mint, test_alliteration=False)

print('--------------------')
print('%s' % (para.chemical_save_name))
print('END at %.2f [s]' % (time.time() - time0))
print('cost %.2f HOURS' % ((time.time() - time0) / 3600.0))
print('--------------------')
