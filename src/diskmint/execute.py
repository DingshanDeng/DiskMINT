#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author: Dingshan Deng @ University of Arizona
# contact: dingshandeng@arizona.edu
# created: 05/18/2023
# new version updated: 08/10/2025
# DiskMINT v1.5.0

import os, sys, copy
import datetime 
import glob
import shutil # , subprocess
import traceback

from . import constants as const
from . import model
from . import disk_density as dd
from . import modelgrid
from . import dustopac

# sys.path.append(os.getcwd())
# print('current working directory at: ', os.getcwd())

try:
    import numpy as np
    from scipy import interpolate
    import scipy.optimize
    from scipy.optimize import curve_fit
except ImportError:
    np = None
    interpolate = None
    scipy = None
    print(' Numpy or Scipy cannot be imported ')
    print(' To use this python module, you need to install Numpy, Scipy')
    print(traceback.format_exc())

# try:
#     import radmc3dPy
# except ImportError:
#     radmc3dPy = None
#     print(' radmc3dPy cannot be imported ')
#     print(' To use this python module, you need to install radmc3dPy')
#     print(traceback.format_exc())


def runmodel(mint, test_alliteration=False):
    """
    Run the DiskMINT model
    
    Args:
        mint (object): The DiskMINT model object containing the data and parameters to be saved.
        test_alliteration (bool): If True, force to run the all the iterations set up for VHSE to see whether it is truly converged at the last few iterations.
        test_io (bool): If True, this run model exe only tests the i/o from the datadir (assuming it was already executed before)
    """

    proceed_code = True
    
    para = mint.parameters
    
    if not os.path.isdir(mint.file_dir):
        os.mkdir(mint.file_dir)

    # Get the absolute path to the directory
    current_dir = os.path.dirname(os.path.abspath(__file__))
    print('current working directory:\n' + current_dir)
    code_start_dir = current_dir

    # Change the current working directory to the path containing the radmc3d model files
    target_dir = mint.file_dir
    os.chdir(target_dir)

    current_dir = os.path.dirname(os.path.abspath(__file__))
    print('current working directory (after changed):\n' + current_dir)

    if mint.bool_MakeDustKappa:
        print('the DiskMINT module does not support making the opacity files here')
        print('Please make your own Dust Kappa and feed into the Model')
        print('Here the setup amin and amax are used to make the neccessary information')
        #  os.system("rm dustkappa_*.inp")
        #  os.system("python makedustkappa.py")
        mint.setup_dust_info(para)
        #  elif len(para.ndsd)
    else:
        print('reading pre set a_ave.inp, fracs.inp, fracs_nb.inp, and ndsd.inp')
        #  then the dust info needs to in read from .inp files
        mint.readin_dust_info(file_dir=mint.file_dir)
        print('pre set a_ave.inp, fracs.inp, fracs_nb.inp, and ndsd.inp are read')

    #
    # Clean old files
    #
    # os.system("rm dust_density_*.dat")
    # os.system("rm dust_temperature_*.dat")
    # os.system("rm gas_temperature_*.new")
    file_list = ['dust_density.inp', 'dust_temperature.dat', 'gas_temperature.inp', 'gas_temperature.new', 'dust_density_*.dat', 'dust_temperature_*.dat', 'gas_temperature_*.new', 'gas_temperature_*.inp', 'dust_density_*.inp', 'dust_density_*.new', 'dust_density_settled.new', 'ratio_g2d_grid_*.dat', 'ratio_g2d_grid_*.new','ratio_g2d_grid.*', 'rhodust_all_probset.inp', 'rhogas_probset.inp', 'rhodust_probset.inp', 'dust_add_rim.new', 'taux.dat', 'tauz_dn.dat', 'tauz_up.dat', 'spectrum_*.out', 'spectrum.out']
    for file_name_i in file_list:
        os.system("rm %s"%(file_name_i))

    #
    # Call problem_setup
    #
    command_mctherm = 'radmc3d mctherm'
    
    # start over or read in the previous files
    if mint.bool_startover:
        #
        # The Module Below relies on reading in the new dust opacity, 
        # so import it after the dust opacity setup and old files removed.
        #
        #  import wrapper_scripts.wrapper_disk_density as dd

        dd.problem_setup(mint)
        print("Done with problem_setup setup")

        nr_redundant_itr = 3
        i_redundant_iter = 1
        
        # comp dust temperature
        os.system(command_mctherm)
        
        fname = 'dust_temperature.dat'
        key_out = 0
        while not key_out:
            if os.path.isfile(fname):
                print('successfully obtained the dust_temperature.dat')
                key_out = 1
            
            elif i_redundant_iter <= nr_redundant_itr:
                i_redundant_iter += 1
                print("run the initial radmc3d mctherm for #%i time"%(i_redundant_iter))
                
                # os.system(command_mctherm)
                dd.get_new_dust_temperature(mint, command=command_mctherm, bool_clean_dusttemp_midplane=mint.bool_clean_dusttemp_midplane)
                
            else:
                print('Something wrong with the initial iteration, check the code before proceeding')
                key_out = 1
                proceed_code = False
    
    # read in the previous files        
    else:
        print('Using previous setups as the initial Gaussian distribution')

        dd.problem_setup(mint)
        print("Done with problem_setup setup")
        print('Copy the files needed from previous runs')

        # file_dir_readin = '../result_radmcfiles/' + prev_setup_name + '/'
        file_dir_readin = mint.data_dir_readin + mint.prev_setup_name + '/'

        print(file_dir_readin)

        if mint.read_in_type == 'Gaussian':
            filein_list = ['stars.inp', 'dustopac.inp', 'fracs.inp', 'aave.inp', 'fracs_numb.inp', 'amr_grid.inp',
                           'wavelength_micron.inp', 'dust_density.inp', 'dust_temperature.dat', 'spectrum_init_GS.out']

        elif mint.read_in_type == 'VHSE':
            filein_list = ['stars.inp', 'dustopac.inp', 'fracs.inp', 'aave.inp', 'fracs_numb.inp', 'amr_grid.inp',
                           'wavelength_micron.inp', 'dust_density_mdust_VHSE_1.dat',
                           'dust_temperature_mdust_VHSE_1.dat', 'spectrum_init_GS.out']

        elif mint.read_in_type == 'VHSEend':
            filein_list = ['stars.inp', 'dustopac.inp', 'fracs.inp', 'aave.inp', 'fracs_numb.inp', 'amr_grid.inp',
                           'wavelength_micron.inp', 'dust_density.inp', 'dust_temperature.dat', 'spectrum_init_GS.out']

        file_list = ['stars.inp', 'dustopac.inp', 'fracs.inp', 'aave.inp', 'fracs_numb.inp', 'amr_grid.inp',
                     'wavelength_micron.inp', 'dust_density.inp', 'dust_temperature.dat', 'spectrum.out']

        # reset the rhogas_probset.inp
        #         rhogas0_prev        = np.loadtxt(file_dir_readin+'rhogas_probset.inp')
        #         ratio_g2d_prev_grid = place_ratio_g2d(mint.ratio_g2d_prev)
        #         rhogas0             = rhogas0_prev/mint.ratio_g2d_prev*mint.ratio_g2d

        #         np.savetxt('rhogas_probset.inp', rhogas0)
        #         mint.rhogas0 = rhogas0

        for ic, file_name_i in enumerate(file_list):
            filename = file_dir_readin + filein_list[ic]
            os.system("cp " + filename + " " + file_name_i)

        #
        # The Module Below relies on reading in the new dust opacity, 
        # so import it after the dust opacity setup and old files removed.
        #
        #  import wrapper_scripts.wrapper_disk_density as dd

        os.system('rm dustkappa_*.inp')
        os.system('cp ' + file_dir_readin + 'dustkappa_*.inp ./')

    command_sed = "radmc3d sed incl %.1f phi %.1f" % (para.incl, para.phi)  # setthreads 10"
    # Compute the sed for the Gaussian Model
    if mint.bool_SED and mint.bool_startover and proceed_code:
        if not mint.bool_VHSE:
            os.system(command_sed)
            filename = "spectrum_init_GS.out"
            os.system("cp spectrum.out " + filename)

    # Save the initial Gaussian disk before VHSE calculation    
    if mint.bool_chemistry and proceed_code:
        dd.chemistry_setup(mint, save_name='GSinit_' + para.chemical_save_name, bool_same_rc_as_radmc3d=mint.bool_same_rc_as_radmc3d)

    # set up the test values to see when to stop the vhse iterations
    difference_test_max_prev_itr = np.inf
    difference_test_mean_prev_itr = np.inf
    
    # start the main loop for solving VHSE
    if mint.bool_VHSE and proceed_code:
        # results above needs to be read in
        # when import the following package
        # so import here

        nloops = mint.n_vhse_loop  # set for testing now
        nr_redundant_itr = int(nloops / 1.0) # set the maximum roduntent iteration
        key_out = 0  # see whether to stop iterations

        # now we start from i_iter 0 for the initial Gaussian distribution
        i_iter = -1 
        i_redundant_iter = 0
        # for i_iter in range(1, nloops + 1):
        # while i_iter < nloops and key_out == 0:
        while key_out == 0:
            
            i_iter = int(i_iter + 1)
            print('start iteration #%i'%(i_iter))
            
            """
            previous version,
            solve VHSE and dust settling both based on 
            the previous density and thermal structure
            """
            
            # #
            # # solve vhse
            # #
            # if i_iter >= 1:
            #     bool_dust_settling_in_iteration_t = mint.bool_dust_settling
            # else:
            #     # we do not apply dust settling in the first iteration
            #     bool_dust_settling_in_iteration_t = False

            # #
            # # compute the dust density structure
            # #
            # """NOTE (developing) now we integrate the dust settling code inside the vhse loop"""
            # get_new_dust_density(mint, bool_dust_settling_in_iteration=bool_dust_settling_in_iteration_t)
            
            """
            01/08/2025
            now in each iteration after i_iter 1,
                first, we include dust settling,
                then, we run radmc3d mctherm to get new Tdust (and Tgas),
                fiannly, we solve vhse on the settled dust and new Tgas
            then go back to first step above
            """
            if i_iter >= 1:
                
                if mint.bool_dust_settling:
                    #
                    # solve the dust settling
                    #
                    dd.set_dust_settling(mint, verbose=True)
                    print('dust and gas are de-coupled, rhodust is updated according to the settling solver at itr #%i'%(i_iter))
                    
                    # move the new files as the input
                    os.system('mv dust_density_settled.new dust_density.inp')
                    os.system("mv ratio_g2d_grid.new ratio_g2d_grid.dat")
                
                # then we need to do the thermal calculation again to get the new dust Temperature
                # os.system("radmc3d mctherm")
                dd.get_new_dust_temperature(mint, command=command_mctherm, bool_clean_dusttemp_midplane=mint.bool_clean_dusttemp_midplane)
            
            # else: continue
            # we do not apply dust settling in the first iteration
            
            """
            # read in the dust density 1 (previous density from last iteration)
            # this is for comparison and see whether the density converges
            """
            # dustdens1 = radmc3dPy.analyze.readData(dtemp=False, ddens=True, gtemp=True)
            dustdens1 = modelgrid.readData(dtemp=False, ddens=True, gtemp=True)
            rhodust1  = dustdens1.rhodust.copy()
            rhogas1   = np.nansum(rhodust1, axis=3) * mint.ratio_g2d
            
            if i_iter >= 1:
                # in the first iteration before solving VHSE
                # Tgas is not available.
                Tgas1_0   = dustdens1.gastemp.copy()
                Tgas1_0   = Tgas1_0[:, :, :, 0]
            
            #
            # solve VHSE
            # in the case when dust settling is turned on (dust-gas de-coupled)
            # this function is getting the new gas density by solving VHSE
            # rather than really get the dust density
            dd.get_new_dust_density(mint, 
                                 bool_dust_settling_in_iteration=mint.bool_dust_settling) # ,
                                #  bool_dust_fragmentation_set=mint.bool_dust_fragmentation,
                                #  bool_dust_radial_drifting_set=mint.bool_dust_radial_drifting)
            
            #
            # Save old files and update dust density file
            #
            # the dust density for the previous iteration
            filename0 = "dust_density_mdust_VHSE_" + str(i_iter) + ".dat"
            os.system("cp dust_density.inp " + filename0)
            # the dust temperature for the previous iteration
            filename1 = "dust_temperature_mdust_VHSE_" + str(i_iter) + ".dat"
            os.system("cp dust_temperature.dat " + filename1)
            # the gas temperature for the previous iteration
            filename2 = "gas_temperature_mdust_VHSE_" + str(i_iter) + ".inp"
            os.system("cp gas_temperature.inp " + filename2)
            
            print('the input files for this iteration (= the files for previous iteration) are saved as\n%s\n%s\n%s'%(filename0, filename1, filename2))
            
            # the ratio g2d grid
            if os.path.exists('ratio_g2d_grid.dat'):
                filename3 = "ratio_g2d_grid_VHSE_" + str(i_iter) + ".dat"
                os.system("cp ratio_g2d_grid.dat " + filename3)
                
                print('%s is also saved'%(filename3))
                
            if mint.bool_clean_dusttemp_midplane:
                # the original dust temperature for the previous iteration
                filename4 = "dust_temperature_originalmidplane_mdust_VHSE_" + str(i_iter) + ".dat"
                os.system("cp dust_temperature_originalmidplane.dat " + filename4)

                print('%s is also saved'%(filename4))

            # 
            # # the NEW ratio g2d
            # filename4 = "ratio_g2d_grid_VHSE_" + str(i_iter) + ".new"
            # os.system("cp ratio_g2d_grid.new " + filename4)
            
            # put the newly calculated dust density and g2d grid to the inp files.
            os.system("mv dust_density.new dust_density.inp")
            os.system("mv ratio_g2d_grid.new ratio_g2d_grid.dat")
            
            """
            # read in the dust density 2 (new density this iteration)
            """
            # dustdens2 = radmc3dPy.analyze.readData(dtemp=False, ddens=True, gtemp=True)
            dustdens2 = modelgrid.readData(dtemp=False, ddens=True, gtemp=True)
            rhodust2  = dustdens2.rhodust.copy()
            rhogas2   = np.nansum(rhodust2, axis=3) * mint.ratio_g2d
            
            if i_iter >= 1:
                # in the first iteration before solving VHSE
                # Tgas is not available.
                Tgas2_0   = dustdens2.gastemp.copy()
                Tgas2_0   = Tgas2_0[:, :, :, 0]
            
            if i_iter >= 1:
                # in the first iteration before solving VHSE
                # Tgas is not available.
                difference_test = np.abs((Tgas2_0 - Tgas1_0)/Tgas1_0)
                difference_test[rhogas2 <= 1e-20] = 0
                difference_test[np.isnan(difference_test)] = 0
                difference_test_max_t = np.nanmax(difference_test)
                difference_test_mean_t = np.nanmean(difference_test)
                print('-------------')
                print('test the Tgas differences of rho between the two iterations')
                print('delta_T/T with mean %.3e and max %.3e' % (difference_test_mean_t, difference_test_max_t))

            difference_test = np.abs((rhogas2 - rhogas1)/rhogas1)
            difference_test[rhogas2 <= 1e-20] = 0
            difference_test[np.isnan(difference_test)] = 0
            difference_test_max_t = np.nanmax(difference_test)
            difference_test_mean_t = np.nanmean(difference_test)
            print('-------------')
            print('test the rhogas differences of rho between the two iterations')
            print('delta_rho/rho with mean %.3e and max %.3e' % (difference_test_mean_t, difference_test_max_t))
            print('-------------')
            
            if not test_alliteration:
                if i_iter <= 1:
                    print('finished the %i iteration for VHSE from an arbitrary gaussian distribution; now move on.'%(i_iter))

                    # here we make sure the max differences are still inf since the first iteration should not count.
                    difference_test_max_prev_itr = np.inf
                    difference_test_mean_prev_itr = np.inf
                
                elif np.nanmean(difference_test) < 1e-20:
                    print('the new density seem to be identical to the previous iteration')
                    
                    if i_redundant_iter < nr_redundant_itr:
                        print('something wrong in this iteration, repeat this iteration #%i'%(i_iter))
                        
                        i_iter = int(i_iter - 1)
                        i_redundant_iter = int(i_redundant_iter + 1)
                        
                        print('this will be the the #%i redundant iteration'%(i_redundant_iter))
                    
                    else:
                        print('a lot iterations have issues, meet the stop criteria, please check the code before proceeding')
                        key_out = 1
                        # proceed_code = False
                
                else:
                    
                    different_criteria = 0.05 # 0.1 and 0.05 tried in the past
                    if np.nanmax(difference_test) < different_criteria:
                        print('the delta_rho/rho is less than criteria of %.2f'%(different_criteria))
                        print('meet the criteria at itr %i' % (i_iter))
                        proceed_code = True
                        key_out = 1 
                        
                    elif (difference_test_mean_t >= difference_test_mean_prev_itr) & (difference_test_max_t >= difference_test_max_prev_itr):
                        print('the differences starts to increase from this iteration')
                        print('the delta_rho/rho in previous iterations are with\nmean %.3e and max %.3e' % (difference_test_mean_prev_itr, difference_test_max_prev_itr))
                        print('the current set up might not be able to decrease this difference and converge into the set up criteria')
                        print('please consider increase nr and nphot')
                        print('code still proceed; be careful with the results.')
                        proceed_code = True
                        key_out = 1
                        
                    else:
                        difference_test_max_prev_itr  = difference_test_max_t
                        difference_test_mean_prev_itr = difference_test_mean_t
            
            print('-------------')
            print("i_iter %i done" % (i_iter))
            print('-------------')
            
            if i_iter >= nloops:
                print('Went through (%i) iterations, larger than the requested number (%i); loop out'%(i_iter, nloops))
                key_out = 1
            
            if key_out == 1:
                print('this is the last iteration in VHSE')
                print('-------------')
                
                # # record the last ratio_g2d_grid as dat instead of new
                # filename = "ratio_g2d_grid.dat"
                # os.system("cp ratio_g2d_grid.new " + filename)
                if not proceed_code:
                    if i_iter > int(nloops / 2.0):
                        proceed_code = True
                        print('Went though (%i) iterations. At least half of the required iterations (requested nloop = %i) have finished, code proceed, but careful on the final results'%(i_iter, nloops))
                        
                    else:
                        proceed_code = False
                        print('The code has some issues, meet the stop criteria, please check the code before proceeding')
                    
                if proceed_code:
                    #
                    # calculate dust temperature
                    # solve the dust temperature in the last iteration
                    # to make sure it is the correct Tdust for this density
                    # os.system("radmc3d mctherm")
                    dd.get_new_dust_temperature(mint, command=command_mctherm, bool_clean_dusttemp_midplane=mint.bool_clean_dusttemp_midplane)
    
    else: 
        """For the case without using VHSE but want dust settling"""
        if mint.bool_dust_settling & proceed_code:
            #
            # solve the dust settling
            #
            dd.set_dust_settling(mint, verbose=True)
            
            # the dust density for the previous iteration
            filename0 = "dust_density_mdust_GS_" + "final_before_dustset.dat"
            os.system("cp dust_density.inp " + filename0)
            # the dust temperature for the previous iteration
            filename1 = "dust_temperature_mdust_GS_" + "final_before_dustset.dat"
            os.system("cp dust_temperature.dat " + filename1)
            # the gas temperature for the previous iteration
            filename2 = "gas_temperature_mdust_GS_" + "final_before_dustset.inp"
            os.system("cp gas_temperature.inp " + filename2)
            # the ratio g2d grid
            if os.path.exists('ratio_g2d_grid.dat'):
                filename3 = "ratio_g2d_grid_GS_" + "final_before_dustset.dat"
                os.system("cp ratio_g2d_grid.dat " + filename3)

            print('the files for this iteration are saved as\n%s\n%s\n%s'%(filename0, filename1, filename2))
            if os.path.exists('ratio_g2d_grid.dat'):
                print('%s is also saved'%(filename3))
            
            # move the new files as the input
            os.system('mv dust_density_settled.new dust_density.inp')
            os.system("mv ratio_g2d_grid.new ratio_g2d_grid.dat")
            
            # then we need to do the thermal calculation again to get the new dust Temperature
            # os.system("radmc3d mctherm")
            dd.get_new_dust_temperature(mint, command=command_mctherm, bool_clean_dusttemp_midplane=mint.bool_clean_dusttemp_midplane)
            
            # then if we need to get the new gas temperature based on this new dust thermal structure
            # dustdens = radmc3dPy.analyze.readData(dtemp=True, ddens=True, binary=False)
            # Tdust = dustdens.dusttemp.copy()
            # Tdust[Tdust < 2.7] = 2.7
            # solve_gastemp(mint, dustdens, Tdust, bool_output_gastemp=True)
            # os.system('mv gas_temperature_direct_output.inp gas_temperature.inp')
        
    bool_frag_andor_drift = mint.bool_dust_fragmentation or mint.bool_dust_radial_drifting
    if bool_frag_andor_drift & proceed_code:
        
        #
        # Save old files and update dust density file
        #
        # the dust density for the previous iteration
        filename0 = "dust_density_mdust_VHSE_" + str(i_iter) + "_before_fragordrift.dat"
        os.system("cp dust_density.inp " + filename0)
        # the dust temperature for the previous iteration
        filename1 = "dust_temperature_mdust_VHSE_" + str(i_iter) + "_before_fragordrift.dat"
        os.system("cp dust_temperature.dat " + filename1)
        # the gas temperature for the previous iteration
        filename2 = "gas_temperature_mdust_VHSE_" + str(i_iter) + "_before_fragordrift.inp"
        os.system("cp gas_temperature.inp " + filename2)
        # the ratio g2d grid
        if os.path.exists('ratio_g2d_grid.dat'):
            filename3 = "ratio_g2d_grid_VHSE_" + str(i_iter) + "_before_fragordrift.dat"
            os.system("cp ratio_g2d_grid.dat " + filename3)

        print('the files for this iteration are saved as\n%s\n%s\n%s'%(filename0, filename1, filename2))
        if os.path.exists('ratio_g2d_grid.dat'):
                print('%s is also saved'%(filename3))
                
        #
        # Put the fragmentation/radial drift
        #
        dd.get_new_dust_density(mint, 
                            bool_dust_vhse=False,
                            bool_dust_settling_in_iteration=False, # both vhse and settling is solved previously
                            bool_dust_fragmentation_set=mint.bool_dust_fragmentation,
                            bool_dust_radial_drifting_set=mint.bool_dust_radial_drifting)
        
        # move the new files as the input
        os.system("mv dust_density.new dust_density.inp")
        os.system("mv ratio_g2d_grid.new ratio_g2d_grid.dat")
        
        # then we need to do the thermal calculation again to get the new dust Temperature
        # os.system("radmc3d mctherm")
        dd.get_new_dust_temperature(mint, command=command_mctherm, bool_clean_dusttemp_midplane=mint.bool_clean_dusttemp_midplane)
        
    """add the dust inner rim"""
    if mint.bool_dust_inner_rim & proceed_code:
        """(devloping)"""
        dd.add_dust_rim(mint, verbose=True)
        
        # the dust density for the previous iteration
        filename0 = "dust_density_mdust_VHSE_" + "final_before_addrim.dat"
        os.system("cp dust_density.inp " + filename0)
        # the dust temperature for the previous iteration
        filename1 = "dust_temperature_mdust_VHSE_" + "final_before_addrim.dat"
        os.system("cp dust_temperature.dat " + filename1)
        # the gas temperature for the previous iteration
        filename2 = "gas_temperature_mdust_VHSE_" + "final_before_addrim.inp"
        os.system("cp gas_temperature.inp " + filename2)
        # the ratio g2d grid
        if os.path.exists('ratio_g2d_grid.dat'):
            filename3 = "ratio_g2d_grid_VHSE_" + "final_before_addrim.dat"
            os.system("cp ratio_g2d_grid.dat " + filename3)
        
        # move the new files as the input
        os.system('cp dust_add_rim.new dust_density.inp')
        os.system("mv ratio_g2d_grid.new ratio_g2d_grid.dat")
        
        # get the new thermal structure
        # os.system("radmc3d mctherm")
        dd.get_new_dust_temperature(mint, command=command_mctherm, bool_clean_dusttemp_midplane=mint.bool_clean_dusttemp_midplane)
    
    if mint.bool_SED & mint.bool_VHSE & proceed_code:
        #
        # Compute the sed
        #
        os.system(command_sed)
        #
        filename = "spectrum_VHSE_final.out"
        os.system("cp spectrum.out " + filename)
    
    elif mint.bool_SED & mint.bool_dust_settling & proceed_code:
        # Compute the sed
        os.system(command_sed)
        #
        filename = "spectrum_GS_settled.out"
        os.system("cp spectrum.out " + filename)

    #
    # After the above computation, save for chemical network
    #
    if mint.bool_chemistry & proceed_code:
        # os.system("python chemical_initial_setup.py")
        dd.chemistry_setup(mint, save_name=para.chemical_save_name, bool_same_rc_as_radmc3d=mint.bool_same_rc_as_radmc3d)
        # os.system("./run_chem.sh")

        # for IM Lup
        # G0Hab_set = 3.439e+10 # this was what was used
        # G0Hab_set = 1.0e11 # this is the one that could give G0Hab_1AU = 1e7 (actually G0Hab_set at 2.5Rsolar is 7.4e10)
        
        dd.runchemistry(mint, chem_code_dir=mint.chem_code_dir, G0Hab_set=para.G0Hab_set) 
        # , bool_save_allab=True) # NOTE:now the bool_save_allab is abandonded. All the chemistry is saving all ab files inside the mint.file_dir

        """
        After the chemical network, the .chem file will be created
        and then use the write_LIME_grid.py file to create the LIME grids
        for running LIME and get the synthetic line flux from LIME
        """

    if mint.bool_savemodel & proceed_code: 
        #
        # by default, the bool_savemodel is on
        # if not true, the models will only stay in the ./data file where the code is executed
        # if true, and the para.chemical_save_dir and para.chemical_save_name are specified, 
        # then all of the results will be saved inside the directory of
        # {chemical_save_dir}/{chemical_save_name}/
        #
        save_dir_t = os.path.join(para.chemical_save_dir, para.chemical_save_name)
        if not os.path.exists(save_dir_t):
        # If it does not exist, create the directory
            os.makedirs(save_dir_t)
            print(f"Directory '{save_dir_t}' created.")
        else:
            print(f"Directory '{save_dir_t}' already exists.")
        
        #
        # save the radmc3d model to the data directory
        #
        dd.savemodel(mint, save_dir=os.path.join(para.chemical_save_dir, para.chemical_save_name), model_dir=mint.file_dir)
        
        print('also saving the parameter files to')
        print('%s'%(os.path.join(para.chemical_save_dir, para.chemical_save_name)))
        
        #
        # save the chemistry files
        #
        if mint.bool_savemodel & mint.bool_chemistry:
            # """save the CO data files and the main output"""
            # data_save_dir = os.path.join(para.chemical_save_dir, para.chemical_save_name) # save_dir # os.path.join(save_dir, 'model_outputs', 'data_files')
            # # Check if the directythory exists
            # if not os.path.exists(data_save_dir):
            #     # If it does not exist, create the directory
            #     os.makedirs(data_save_dir)
            #     print(f"Directory '{data_save_dir}' created.")
            # else:
            #     print(f"Directory '{data_save_dir}' already exists.")

            # file_move_list = ['COinitgrid-*.dat', 'COendgrid-*.chem']
            # for file_t in file_move_list:
            #     command_t = 'mv ' + os.path.join(mint.file_dir, file_t) + ' ' + data_save_dir
            #     print('executing command: %s'%(command_t))
            #     os.system(command_t)
            
            # Save the CO data files and the main output
            data_save_dir = os.path.join(para.chemical_save_dir, para.chemical_save_name)

            # Create the directory if it doesn't exist
            os.makedirs(data_save_dir, exist_ok=True)
            print(f"Directory '{data_save_dir}' is ready.")

            # Patterns of files to move
            file_patterns = ['COinitgrid-*.dat', 'COendgrid-*.chem']
            for pattern in file_patterns:
                matching_files = glob.glob(os.path.join(mint.file_dir, pattern))
                if not matching_files:
                    print(f"No files found for pattern: {pattern}")
                for file_path in matching_files:
                    dst_path = os.path.join(data_save_dir, os.path.basename(file_path))
                    print(f"Moving: {file_path} → {dst_path}")
                    shutil.move(file_path, dst_path)
                
            """save all ab files for all the molecules that are included in the current chemical network -- this is most to help specify the name of the file"""
            # but now we are also going to record the name of the chemical_network we just ran so that we do not mess up with the anything
            src = os.path.join(mint.file_dir, 'chemistry')
            # dst = os.path.join(mint.file_dir, f'chemistry-{os.path.basename(mint.chem_code_dir)}')
            dst_base = os.path.join(mint.file_dir, f'chemistry-{os.path.basename(mint.chem_code_dir)}')
            dst = dst_base

            # If destination exists, append timestamp until unique
            if os.path.exists(dst):
                timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
                dst = f"{dst_base}_copy_{timestamp}"

            print(f"Moving: {src} → {dst}")
            shutil.move(src, dst)
        
        # NOTE (Developing)
        # save the parameter files this is now saved in the savemodel function above
        # mint.parameters.write_parameters_to_csv(para.chemical_save_name+'_parameters', directory=save_dir_t, extension='.csv')
        
        # save all of the chemistry outputs
        
    print('----------------------')
    print('finished the model for')
    print('%s'%(para.chemical_save_name))
    print('----------------------')
    
    return proceed_code

