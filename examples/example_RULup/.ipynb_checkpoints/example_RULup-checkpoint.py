#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author: Dingshan Deng @ University of Arizona
# contact: dingshandeng@arizona.edu
# created: 05/22/2023

"""
The model example of RU Lup best-fit model
"""

import os, sys, copy
sys.path.append(os.getcwd())
package_position = "/home/dingshandeng/github/DiskModeling/0-DiskMINT/" 
sys.path.append(package_position)

# Get the absolute path to the directory
current_dir = os.path.dirname(os.path.abspath(__file__))
print('current working directory:\n'+current_dir)

# Change the current working directory to the path containing the radmc3d model files
work_dir = "/home/dingshandeng/github/DiskModeling/0-DiskMINT/examples/example_RULup/"
target_dir = work_dir 
os.chdir(target_dir)

current_dir = os.path.dirname(os.path.abspath(__file__))
print('current working directory (after changed):\n'+current_dir)

import diskmint.model as model
import diskmint.disk_density as dd
import diskmint.constants as const

import numpy as np

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
# para.read_parameters_from_csv(filename=os.path.join(work_dir,'test_parameters.dat'))

#
# Have a look on all used parameters
#
print(para.parameters)
# print('disk dust mass: %.3e'%(para.mdiskd/const.ms))
# print('stellar radius: %.3f'%(para.rstar/const.rs))

#
# model file name
#
# stellar model (the stellar model needs to have uv spectral)
# so that the UV field can be correctly computed
#
para.fmodel_filename = 'RULupSpecModel_wuv.inp'
print(para.fmodel_filename)

# the directory to save the chemical input/output files
para.chemical_save_dir  = '/home/dingshandeng/data/diskmass/DiskMINTv3p2_results/RULup_chemistry/' 

# the save name of this current model
para.chemical_save_name = 'testmodel_RULup_VHSE_bestfit_bestgtd_Jun11_t0'

# what to use to run the model
para.nphot = 1e7 # photon number in RADMC3D
para.nthreads = 24 # threads number in RADMC3D
para.ratio_g2d_global = 30.0 # the g2d ratio that will be used in the model

#
# make a copy of the model parameters that is used
#
para.write_parameters_to_csv(os.path.join(work_dir, para.chemical_save_name+'_parameters.dat'))

#
# You need to make your own dust opacity files in radmc3d format and bring them to the working directory
# and specify opacity file names
#
para.dustopacname_1 = 'DSHARP_repWD01_Dec6'
# currently it supports two different compositons and can be put here
#  para.dustopacname_2 = ''


# 
# besides, you need to make your own files for the info of those opacities
# the info that is needed including
# fracs: mass fractions of dust species
# a_ave: average dust radius of dust dpecies
# fracs_nb: the number fractions of dust species
# ndsd: nd(a) * pi * a**2
# they need to be the name of 'fracs.inp', 'a_ave.inp', 'fracs_nb.inp', 'ndsd.inp'
# if you don't want to calculate all these files by yourself, then
# you need to specify the dust size distribution info correctly in the model setup
# as the example below
#
bool_MakeDustKappa = True
if bool_MakeDustKappa:
    #
    # number of dust species
    # to setup our VHSE models
    # [NOTE] here we support two different dust species
    # which can be further divided into different subgroups
    # with different sizes.
    para.nr_dust_1         = 20
    para.nr_dust_2         = 0
    # automatically set the total nr of dust
    para.dust_spec_nr      = para.nr_dust_1 + para.nr_dust_2

    #
    # use following parameters to make dust kappa files
    #
    # here a is the radius of dust grain in the unit of cm.
    #
    # minimum size of the dust
    # use 100 AA as 1.0e-6 so that dust is not as small as molecule
    para.amin_1             = 1.0e-6
    para.amin_2             = 1.0e-6
    #
    # maximum size of the dust
    para.amax_1             = 0.3
    para.amax_2             = 0.3
    #
    # automatically set the overall smallest
    # and largest sizes among the two
    para.amin_all           = np.min([para.amin_1, para.amin_2])
    para.amax_all           = np.max([para.amax_1, para.amax_2])
    #
    # power law index of the dust size distribution
    para.pla_dustsize       = 3.5

    #
    # bulk density
    # the bulk density should be larger instead of this value
    # but it seems like it is going to be cancelled so double check the code
    # the rhobulk is used in chemical network
    # prev used 2.5 (too small)
    #
    para.rhobulk            = 3.224   # assume in the unit of [g cm^-3]
#
# where is the checmial network 
#
chem_code_dir = '/home/dingshandeng/github/DiskModeling/0-DiskMINT/chemistry/with_CO2_ice/'

#
# set up the DiskMINT mint model object
#
mint = model.Mint(para, file_dir=work_dir)

#
# defaultly the dust opacity info is pre set so read in the neccessary info
# from the .inp files
# but here in case you want to use the feature we provide to generate those
#
mint.bool_MakeDustKappa = bool_MakeDustKappa

# 
# defaultly the chemistry is not run (so bool_chemistry = False)
# here we want to run the chemistry so change to True
#
mint.bool_chemistry = True
mint.chem_code_dir  = chem_code_dir

#
# have a look on the set up g2d in this model on the screen to double check
# (or print out other parameters as needed)
#
print('the setup g2d global value is %.2f'%(mint.parameters.ratio_g2d_global))

#
# Run model
#
dd.runmodel(mint)


