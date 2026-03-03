#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author: Dingshan Deng @ University of Arizona
# contact: dingshandeng@arizona.edu
# created: 04/02/2025
# updated: 02/23/2026

"""
model a random disk around 0.5 ms star that is similar to HH30.

to run this example code, use
screen -S $session_name # change $session_name to any screen session names you like
conda activate /tank/data/research/conda_env/diskmint_env
python -u 0-model_similar_to_HH30_example.py 2>&1 | tee -a output_20260223_t0.log
"""

import os, sys, copy, glob, shutil
import numpy as np
import time

time0 = time.time()

#
# Temporarily add the package position to the Python Path
# If it is already added to your Python path, then you can ignore this line
#
# import platform

# if platform.node() == 'garnet':
#     package_position = "/home/dingshandeng/github/DiskModeling/0-DiskMINT/"
# elif platform.node() == 'fld':
#     package_position = "/tank/data/research/software/DiskMINT/"
# #  package_position = "/Users/dingshandeng/github/DiskModeling/0-DiskMINT/"
# sys.path.append(os.path.join(package_position, 'src'))

import diskmint.constants as const
import diskmint.model as model
import diskmint.disk_density as dd
import diskmint.execute as exe

# import py_IMLup_utiles as utils
current_dir = os.getcwd()
utils_dir = os.path.join(current_dir, "..", "..", "example_utils") # "/tank/data/research/ice_in_disk/0-Utils"
sys.path.append(utils_dir)
import diskmint_utils as utils

#
# I. Set up the directories and the names
#
# (0) Working Direcotry: this is the place that you want to keep every files together
# working_dir = os.path.join(package_position, 'examples', 'example_RULup')
working_dir = current_dir
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
read_dir = os.path.join(working_dir, 'input')
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
# please check example_RULup_parameters.dat for detail
file_parameters = 'model_parameters_0p5ms_star'
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
#
bool_dust_settling           = True
#
# whether to run the chemical network
#
bool_chemistry               = True

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

# the directory to save a copy of the chemical input/output files
para.edit_parameter("chemical_save_dir", new_value = save_dir)

"""
# searching for the best-fit combination for the dust parameters
# example best-fit combination was
amax_set = 0.1
pla_set  = 3.5
"""
# # also change the plh to get a better fitting in initial Gaussian model.
# para.edit_parameter("plh", new_value = 1.20) # power-law index for h in Gaussian model. 1.20 is used for DIANA dust. 1.125 used for dsharp dust

"""
# Setting up the model parameters and the model name
"""

model_name_prefix = 'diskmint_similar_to_HH30_example'
model_date = '20260223'
ratio_g2d_t = 100.0 # the gas to dust mass ratio = Mgas / Mdust

i_start = 0
i_model = i_start-1
for mdisk_dust_t, \
    pla_set_t, \
    amax_set_t, \
    in zip([1.0e-4], 
              [3.5], 
              [0.1]
           ):
 
    bool_same_rc_as_radmc3d = True
    
    # record the model number
    i_model = int(i_model + 1)
     
    para.edit_parameter("mdiskd", new_value = mdisk_dust_t) # 2.3e-3 ms for the model change g2d ratio; 1.0e-3 ms for the model that changes the sigmad.
    
    print('-'*10)
    time1 = time.time()
    print('start at %.2f [s]'%(time.time() - time0))
    print('model #%i'%(i_model))
    print('the parameters changed from the csv file in this run:')
    print('mdiskd: %.2e ms'%(para.get_parameter_value('mdiskd')/const.ms))
    print('pla_dustsize: %.2f'%(para.get_parameter_value('pla_dustsize')))
    print('amax_all: %.2f cm'%(para.get_parameter_value('amax_all'))) 
    print('-'*10)
    
    # And apply ratio_g2d_t and dust params to para:
    para.edit_parameter("ratio_g2d_global", new_value=ratio_g2d_t)
    para.edit_parameter("pla_dustsize", new_value=pla_set_t)
    para.edit_parameter("amax_all", new_value=amax_set_t)
    
    """
    #
    # Setting up the path for save data and run RADMC-3D
    #
    """
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

    # Copy input files
    files_input_list = [
        para.fmodel_filename,
        f"dustkappa_{para.dustopacname_1}_*",
        f"dustkappa_{para.dustopacname_2}_*",
        "aave.inp", "fracs.inp", "fracs_numb.inp", "ndsd.inp",
    ]

    copied = 0
    for file_pattern in files_input_list:
        for source_path in glob.glob(os.path.join(read_dir, file_pattern)):
            if os.path.isfile(source_path):
                shutil.copy(source_path, data_dir)
                copied += 1
    print("Copied %d input files to data directory." % copied)

       
    """ 
    # change the used dust parameters
    #
    # You need to make your own dust opacity files in radmc3d format and bring them to the working directory
    # and specify opacity file names
    #
    """
    
    bool_MakeDustKappa_andFiles = True
    amin_set = 1e-6 # default of DSHARP is 1e-6 # cm # 1e-5 adopted for DIANA that is derived from Francischi # 5e-6 is the default amin in optool
    amax_set = amax_set_t # 0.07 # 0.1 # 0.5 # 1.0 # cm
    # 0.1 is used for DIANA dust
    pla_set = pla_set_t # 3.7
    
    """
    # WD01 + DSHARP dust
    """
    # amax_set = 0.2 # 0.1 # 0.5 # 1.0 # cm
    # amid = 1e-2
    # rhobulk_set = 3.224 # 2.08 for DIANA # 3.224 for DSHARP
    # nr_dust_list = [15, 5]
    # nr_dust_set = nr_dust_list[0] + nr_dust_list[1] # 3 # 20 # 20
    
    """
    # using the optool to make DIANA dust
    """
    porosity_t = 0.0 # 0.25 
    mass_fraction_carbon = 0.13 # 0.10 # default 0.13
    mass_fraction_silicate = 1-mass_fraction_carbon
    rho_carbon = 1.8
    rho_silicate = 3.01
    mass_fraction_carbon = 0.13
    mass_fraction_silicate = 1-mass_fraction_carbon
    porous_volumn = porosity_t # 0.0 # 0.25
    rhobulk = (1 - porous_volumn)/(mass_fraction_carbon/rho_carbon + mass_fraction_silicate/rho_silicate)
    # rhobulk_set = 2.70 # 2.08 for DIANA # 3.224 for DSHARP # 2.70 for solid DIANA dust
    rhobulk_set = rhobulk
    # the mass fractions below only work for the optool dust
    nr_dust_set = 20 
    
    """
    # change the dust parameters names
    """ 
    
    # para.dustopacname_1 = 'DSHARP_default_amax0p2'
    # # currently it supports two different compositons and can be put here
    # para.dustopacname_2 = 'DSHARP_default_amax0p2'
    
    # the optool DIANA dust
    dustopacname_base = 'optool_%s_amin1e-6_amax%.2f_pla%.2f_rhobulk%.1f'%(para.dustopacname_1, amax_set_t, pla_set_t, rhobulk_set)
    dustopacname_base = dustopacname_base.replace('.', 'p')
    para.dustopacname_1 = dustopacname_base
    para.dustopacname_2 = para.dustopacname_1
    
    # para.dustopacname_1 = 'DSHARP_WD01_n_default_amax0p1'
    # para.dustopacname_1 = 'DSHARP_default_amax0p1'
    
    """
    #
    # the new dust opac function
    #
    """
    para = utils.set_dust_parameters_to_nrdust1(para, amin_set=amin_set, amax_set=amax_set, nr_dust_set=nr_dust_set, pla_set=pla_set, rhobulk_set=rhobulk_set)
    
    """
    # the function below is for setting up the dust with a inner disk separately
    # or used to set up two different dust species
    """
    # nr_dust_set_2 = nr_dust_list[1]
    # para = utils.set_dust_parameters_to_nrdust_1n2(para, amin_set_list=[amin_set, amid], amax_set_list=[amid, amax_set], nr_dust_set_list=[nr_dust_set, nr_dust_set_2], pla_set=pla_set, rhobulk_set=rhobulk_set)
    
    if bool_MakeDustKappa_andFiles:
        os.system("rm " +os.path.join(data_dir, "dustkappa_*.inp"))
        
        print('-'*10)
        print('Making new dust kappa files with the new dust parameters...')
        print('-'*10)
        
        # utils.wrapper_dsharp_opac(para.amin_1, amid, amid, para.amax_1, nr_dust_list[0], nr_dust_list[1], para.pla_dustsize, para.dustopacname_1, datafile = '/home/dingshandeng/github/DiskModeling/3-opacities/opacities_DSHARP/notebooks/data/smallergrains_opacities.npz', bool_makefigures=False)
        
        # # save the time for making new opacities
        # # copy the DSHARP with WD01 with amid=1e-2
        # os.system("cp ./opacities_temp/dustkappa_*.inp ./")
        
        utils.wrapper_optool_opac(data_dir=data_dir, dust_opac_name=para.dustopacname_1,amin=para.amin_1, amax=para.amax_1, pla=para.pla_dustsize, nr_dust=para.dust_spec_nr, mass_frac_1=mass_fraction_silicate, mass_frac_2=mass_fraction_carbon, porosity=porosity_t)
    
    #
    # Have a look on all used parameters
    #
    # print(para.parameters)
    print('disk g2d: %s'%(str(para.ratio_g2d)))
    print('disk dust mass: %.3e'%(para.mdiskd/const.ms))
    print('stellar radius: %.3f'%(para.rstar/const.rs))
    print('stellar mass: %.3f'%(para.mstar/const.ms))

    # print(para.get_parameter_value('chemical_save_name'))

    """define the save name"""
    
    name_of_this_model = '%s_%s_t%i_mdust%.1ems_gtd%.0f_pla%.2f_amax%.2f'%(model_name_prefix, model_date, i_model, para.get_parameter_value('mdiskd'), para.get_parameter_value('ratio_g2d_global'), para.get_parameter_value('pla_dustsize'), para.get_parameter_value('amax_all'))
    name_of_this_model = name_of_this_model.replace('.', 'p')

    # the name of this current model
    para.edit_parameter("chemical_save_name", new_value = name_of_this_model)
    
    #
    # make a copy of the model parameters that is used (before running the model)
    #
    para.write_parameters_to_csv(para.chemical_save_name+'_parameters_setup', directory=data_dir, extension='.csv')

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
    mint.chem_code_dir                = chem_code_dir
    print(mint.chem_code_dir)

    """
    #
    # Run model
    #
    """
    exe.runmodel(mint, test_alliteration=False)

    """save the CO data files"""
    data_save_dir = os.path.join(working_dir, 'model_outputs', 'data_files')
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
    print('cost %.2f HOURS'%((time.time() - time1)/3600.0))
    print('cost %.2f CPU HOURS'%((time.time() - time1)/3600.0 * para.nthreads))
    print('--------------------')
    
