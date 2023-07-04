#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author: Dingshan Deng @ University of Arizona
# contact: dingshandeng@arizona.edu
# created: 05/22/2023

"""
The model example of RU Lup best-fit model
for DiskMINT v1.0beta
"""

"""
**NOTE**
Please edit the following `Yourpath` and `Yourpath_saving_files`
according to your system directories.

Please follow the instructions in the README file to use the DiskMINT code
and this example.

Go through this Python3 script before running the code is highly recommended
"""

import os, sys, copy
import numpy as np

#
# Temporarily add the package position to the Python Path
# If it is already added to your Python path, then you can ignore this line
#
## Please Make sure radmc3dPy is installed, and the PATH is added to Python Path
#  radmc3dPy_position = "Yourpath/radmc3d-2.0/python/radmc3dPy/"
#  sys.path.append(radmc3dPy_position)
#
## Import the diskmint Python3 module
package_position = "Yourpath/DiskMINT/"
sys.path.append(package_position)
import diskmint.constants as const
import diskmint.model as model
import diskmint.disk_density as dd

#
# I. Set up the directories and the names
#
# (0) Working Direcotry: this is the place that you want to keep every files together
working_dir = os.path.join(package_position, 'examples', 'example_RULup')
#
# (1) Data directory:
# where you want to save all of the model input and output files?
# it would be a lot of data files so make sure there would be enough space
# it can be a different place other than where you are running this script
data_dir = os.path.join(working_dir, 'data')
#
# (2) Saving files directory:
# a place where a copy of the chemical and RADMC-3D input/output files will be saved
# NOTE: This dir can also be changed in the parameters.dat file which is called 'chemical_save_dir'
save_dir = "Yourpath_saving_files"
#
# (3) Name of this model:
# Properly name each model you run to prevent any confusions!
# Change name for each run for the model is recommended!
# NOTE: This name can also be changed in the parameters.dat file which is called 'chemical_save_name'
name_of_this_model = 'example_model_RULup_VHSE_bestfit'
#
# (4) The parameter file:
# The code takes a parameter file to set up the model parameters
# please check example_RULup_parameters.dat for detail
file_parameters = os.path.join(working_dir, 'example_RULup_parameters.dat')
#
# (5) (Optional) Chemistry code directory:
# where is the checmial network, only needed if you want to use the chemical network we provide
chem_code_dir = os.path.join(package_position, 'chemistry', 'reducedRGH22')

#
# II. The Options of the Model
#
# whether created new Dust Kappa File
# so that want to calculate dust parameters
# such as fracs and aave instead of reading in
bool_MakeDustKappa = True
#
# whether to compute the SED
bool_SED           = True
#
# whether to use the VHSE distribution
bool_VHSE          = True
# the number of loops it tries to go through to solve the VHSE
# the iteration will automatically stop if it converges (normally <10 iterations)
# but you still would like to set a maximum limit here
n_vhse_loop        = int(20)
#
# whether to run the chemical network
bool_chemistry     = True
#
# whether to save the radmc3d model results in another folder
bool_savemodel     = True

#
# Setting up the path for save data and run RADMC-3D
#
# Get the absolute path to the directory
current_dir = os.path.dirname(os.path.abspath(__file__))
print('current working directory:\n'+current_dir)
# Change the current working directory to the path containing the radmc3d model files
target_dir = data_dir
# create the dir if it does not exist
if not os.path.exists(target_dir):
    os.makedirs(target_dir)
    print(f"Directory '{target_dir}' created.")
os.chdir(target_dir)
current_dir = os.path.dirname(os.path.abspath(__file__))
print('current working directory (after changed):\n'+current_dir)

#
# In the current version, the default parameters are the best fit VHSE model for RU Lup,
# and the default parameters are set up when the Parameters() function is called
#
para = model.Parameters()

#
# the model parameters can also be read from a pre-set up data files (following the example below)
# the data file is in csv format but can also accept comments that starts with '#' in that line
# as a example, see test_parameters.dat
# or the parameter data file that would be write out later in this script named as para.chemical_save_name+'_parameters.dat'
#
para.read_parameters_from_csv(filename=file_parameters)

#
# (Optional) Model Parameters can also be changed here in the script
# To do that, use the script below
# what to use to run the model
"""
para.edit_parameter("nphot", new_value=1e7) # photon number in RADMC3D
para.edit_parameter("nthreads", new_value=24) # threads number in RADMC3D
para.edit_parameter("ratio_g2d_global", new_value=30.0) # the global gas-to-dust mass ratio for the model
"""

# Have a look on all used parameters
print(para.parameters)

# the directory to save a copy of the chemical input/output files
para.edit_parameter("chemical_save_dir", new_value = save_dir)

# the name of this current model
para.edit_parameter("chemical_save_name", new_value = name_of_this_model)

#
# make a copy of the model parameters that is used
#
para.write_parameters_to_csv(os.path.join(working_dir, para.chemical_save_name+'_parameters.dat'))

#
# set up the DiskMINT mint model object using the parameters we just set
# also tell the object where is the data directory for running all the models
#
mint = model.Mint(para, file_dir=data_dir)

#
# read in the bools we set
#
mint.bool_MakeDustKappa = bool_MakeDustKappa
mint.bool_SED           = bool_SED
mint.bool_VHSE          = bool_VHSE
mint.n_vhse_loop        = n_vhse_loop
mint.bool_chemistry     = bool_chemistry
mint.bool_savemodel     = bool_savemodel
mint.chem_code_dir      = chem_code_dir
print(mint.chem_code_dir)

#
# have a look on the set up g2d in this model on the screen to double check
# (or print out other parameters as needed)
#
# print('the setup g2d global value is %.2f'%(mint.parameters.ratio_g2d_global))

#
# Run model
#
dd.runmodel(mint)
