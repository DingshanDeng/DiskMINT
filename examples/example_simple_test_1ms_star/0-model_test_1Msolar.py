#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author: Dingshan Deng @ University of Arizona
# contact: dingshandeng@arizona.edu
# created: 01/24/2025

"""
model the disk similar to RGH22
"""

import os, sys, copy, shutil, glob
import numpy as np
import time

time0 = time.time()

#
# Temporarily add the package position to the Python Path
# If it is already added to your Python path, then you can ignore this line
#
# package_position = "/home/dingshandeng/github/DiskModeling/0-DiskMINT/"
package_position = "/home/dingshandeng/github/DiskMINT/"
sys.path.append(os.path.join(package_position, 'src'))
import diskmint.constants as const
import diskmint.model as model
import diskmint.disk_density as dd
import diskmint.execute as exe

#
# I. Set up the directories and the names
#
# (0) Working Direcotry: this is the place that you want to keep every files together
working_dir = current_dir = os.getcwd()
#
# (1) Data directory:
# where you want to save all of the model input and output files?
# it would be a lot of data files so make sure there would be enough space
# it can be a different place other than where you are running this script
data_dir = os.path.join(working_dir, 'data')
#
# (2) Input files directory:
# the place to put all the dustkappa_* files, "aave.inp", "fracs.inp", "fracs_numb.inp", "ndsd.inp",
# and the model spectra para.fmodel_filename
read_dir = os.path.join(working_dir, 'input')
#
# (2) Saving files directory:
# a place where a copy of the chemical and RADMC-3D input/output files will be saved
# NOTE: This dir can also be changed in the parameters.dat file which is called 'chemical_save_dir'
save_dir = os.path.join(working_dir, 'output')
#
# (3) Name of this model:
# Properly name each model you run to prevent any confusions!
# Change name for each run for the model is recommended!
# NOTE: This name can also be changed in the parameters.dat file which is called 'chemical_save_name'
# name_of_this_model = 'testmodel_IMLup_VHSE_fitSED_Rtap%.0f_modgtd%.0f_Jul4_t3'%(para.get_parameter_value('Rtap'), para.get_parameter_value('ratio_g2d_global'))
"""
Name will be defined later in this model
"""
#
# (4) The parameter file:
# The code takes a parameter file to set up the model parameters
file_parameters = 'model_parameters_1msolar'
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
bool_MakeDustKappa           = True
#
# whether to compute the SED
bool_SED                     = True
#
# whether to use the VHSE distribution
bool_VHSE                    = True
# the number of loops it tries to go through to solve the VHSE
# the iteration will automatically stop if it converges (normally <10 iterations)
# but you still would like to set a maximum limit here
n_vhse_loop                  = int(20) # int(10)
# 
# whether to enalble dust settling
# (Developing)
bool_dust_settling           = True
#
# whether to run the chemical network
bool_chemistry               = True
#
# whether to save the radmc3d model results in another folder
bool_savemodel               = True

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
# in the current version, the default parameters are the best fit VHSE model for RU Lup
#
para = model.Parameters()

#
# the model parameters can also be read from a pre-set up data files (following the example below)
# the data file is in csv format but can also accept comments that starts with '#' in that line
# as a example, see test_parameters.dat
# or the parameter data file that would be write out later in this script named as para.chemical_save_name+'_parameters.dat'
#
para.read_parameters_from_csv(filename=file_parameters, directory=working_dir, extension='.csv')

#
# (Optional) Model Parameters can also be changed here in the script
# To do that, use the script below
# what to use to run the model
"""
# NOTE: these are defined later in this script
para.edit_parameter("nphot", new_value=1e7) # photon number in RADMC3D
para.edit_parameter("nthreads", new_value=24) # threads number in RADMC3D
para.edit_parameter("ratio_g2d_global", new_value=30.0) # the global gas-to-dust mass ratio for the model
"""

#
# Copy the input files into the data directory
#

# List of input files, including wildcard patterns
files_input_list = [
    para.fmodel_filename, 
    f"dustkappa_{para.dustopacname_1}_*",
    f"dustkappa_{para.dustopacname_2}_*",
    "aave.inp", "fracs.inp", "fracs_numb.inp", "ndsd.inp",
]

# Copy the files
for file_pattern in files_input_list:
    # Resolve the wildcard pattern into actual file paths
    source_paths = glob.glob(os.path.join(read_dir, file_pattern))
    
    # Copy each resolved file
    for source_path in source_paths:
        # Ensure the file exists
        if os.path.isfile(source_path):
            destination_path = os.path.join(data_dir, os.path.basename(source_path))
            shutil.copy(source_path, destination_path)
            print(f"Copied {source_path} to {destination_path}")
        else:
            print(f"File not found: {source_path}")

i_model = 0

model_date_input = '2025901'
model_mark_input = 'grid1msolar'

i_start = 0
i_model = i_start-1

# This is an example of building up a few models in one scripts
for ratio_mdisk_to_star, \
    ratio_g2d_t, \
    in zip([0.01, 0.01, 0.01],
           [10.0, 100.0, 1000.0], 
           ):
 
    # record the model number
    i_model = int(i_model + 1)
    
    # change some variables here    
    rtap_set = 100
    para.edit_parameter("Rtap", new_value = rtap_set)
    
    model_date = model_date_input.copy()
    model_mark = model_mark_input.copy()
    
    if 'DSHARP' in para.dustopacname_1:
        model_mark = model_mark + '_DSHARPdust'
    if bool_dust_settling:
        model_mark = model_mark + '_dustset'
    # if bool_dust_fragmentation:
    #     model_mark = model_mark + '_dustfrag%i'%(vel_frag_t)
    # if bool_dust_radial_drifting:
    #     model_mark = model_mark + '_dustdrift' 
    # if bool_dust_inner_rim:
    #     model_mark = model_mark + '_dustrim'
    # if bool_same_rc_as_radmc3d: # Note: by default I am using the same rc as radmc3d
    #     model_mark = model_mark + '_samerc_in_chem'
    
    # define the name
    name_of_this_model = 'diskmintmodel%s_t%i_%s_mgas%.1ems_mdust%.1ems_gtd%.0f_rc%.1f_alphav%.1e'%(model_date, i_model, model_mark, mdisk_gas_t, para.get_parameter_value('mdiskd'), para.get_parameter_value('ratio_g2d_global'), para.get_parameter_value('Rtap'), para.get_parameter_value('visc_alpha'))
    name_of_this_model = name_of_this_model.replace('.', 'p')
    
    print('--------------------')
    time1 = time.time()
    print('start at %.2f [s]'%(time.time() - time0))
    print('at time %s'%(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))) 
    print('model #%i'%(i_model))
    print('--------------------')
    print('%s'%(name_of_this_model))
    print('--------------------')
    
    #
    # set up the DiskMINT mint model object
    #
    mint = model.Mint(para, file_dir=data_dir)

    #
    # read in the bools we set
    #
    mint.bool_MakeDustKappa           = bool_MakeDustKappa
    mint.bool_SED                     = bool_SED
    mint.bool_VHSE                    = bool_VHSE
    mint.n_vhse_loop                  = n_vhse_loop
    mint.bool_dust_settling           = bool_dust_settling
    mint.bool_chemistry               = bool_chemistry
    mint.bool_savemodel               = bool_savemodel
    mint.chem_code_dir                = chem_code_dir
    print(mint.chem_code_dir)
    
    # (new features under development) 
    mint.bool_temp_decouple           = bool_temp_decouple = False
    mint.bool_dust_fragmentation      = bool_dust_fragmentation = False
    mint.bool_dust_radial_drifting    = bool_dust_radial_drifting = False
    mint.bool_dust_inner_rim          = bool_dust_inner_rim = False
    mint.bool_same_rc_as_radmc3d      = bool_same_rc_as_radmc3d = True

    """
    #
    # Run model
    #
    """
    exe.runmodel(mint, test_alliteration=False)

    """save the CO data files"""
    data_save_dir = os.path.join(para.chemical_save_dir, para.chemical_save_name)
    # data_save_dir = save_dir # os.path.join(save_dir, 'model_outputs', 'data_files')
    # Check if the directythory exists
    if not os.path.exists(data_save_dir):
        # If it does not exist, create the directory
        os.makedirs(data_save_dir)
        print(f"Directory '{data_save_dir}' created.")
    else:
        print(f"Directory '{data_save_dir}' already exists.")

    file_move_list = ['COinitgrid-*.dat', 'COendgrid-*.chem']
    for file_t in file_move_list:
        command_t = 'mv ' + os.path.join(data_dir, file_t) + ' ' + data_save_dir
        print('executing command: %s'%(command_t))
        os.system(command_t)

    print('--------------------')
    print('model #%i'%(i_model))
    print('%s'%(para.chemical_save_name))
    print('END at %.2f [s]'%(time.time() - time0))
    print('at time %s'%(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))) 
    print('cost %.2f HOURS'%((time.time() - time0)/3600.0))
    print('cost %.2f CPU HOURS'%((time.time() - time0)/3600.0 * para.nthreads))
    print('--------------------')
