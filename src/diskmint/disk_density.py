#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author: Dingshan Deng @ University of Arizona
# contact: dingshandeng@arizona.edu
# created: 05/18/2023
# new version updated: 07/11/2025
# DiskMINT v1.5.0

"""
Here I will define the key functions to calculate the disk density
in VHSE and also set up the initial inputs for chemical network

This will be based on the disk density file in v3.1alpha
as well as the initial setups that was called as
wrapper_model_parameters.py
"""

import os, sys, copy, time
import shutil, subprocess
import traceback

from . import constants as const
from . import model
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

"""
functions
"""
def save_3dgrid_to_dat(array, filename):
    """
    Saves a 3D NumPy array to a .dat file, including the shape as the first line.
    
    Args:
        array (numpy.ndarray): The 3D NumPy array to save.
        filename (str): The name of the .dat file to save the array in.
    """
    # Save the shape of the array
    original_shape = array.shape

    # Reshape the 3D array to 2D
    array_2d = array.reshape(-1, array.shape[-1])

    # Open the file for writing
    with open(filename, 'w') as f:
        # Write the shape as the first line
        f.write(' '.join(map(str, original_shape)) + '\n')
        
        # Save the 2D array to the file
        np.savetxt(f, array_2d)
        
    return 1


def load_3dgrid_from_dat(filename):
    """
    Loads a 3D NumPy array from a .dat file, assuming the first line contains the shape.
    
    Args:
        filename (str): The name of the .dat file to load the array from.
    
    Returns:
        numpy.ndarray: The loaded 3D NumPy array.
    """
    # Open the file for reading
    with open(filename, 'r') as f:
        # Read the first line to get the shape
        shape_str = f.readline().strip()
        original_shape = tuple(map(int, shape_str.split()))
        
        # Load the rest of the data
        array_2d = np.loadtxt(f)

    # Reshape the 2D array back to the original 3D shape
    array_3d = array_2d.reshape(original_shape)
    
    return array_3d

    

def grid_refine_inner_edge(x_orig, nlev, nspan):
    """
    Refines the grid by adding points to the inner edges of a given grid.

    Args:
        x_orig (ndarray): Original grid points as a 1D numpy array.
        nlev (int): Number of refinement levels.
        nspan (int): Number of spans or intervals to refine.

    Returns:
        ndarray: Refined grid points as a 1D numpy array.

    Description:
    This function takes an original grid, specified by the `x_orig` array, and refines it by adding
    points to the inner edges of the grid. The refinement is performed `nlev` times, with `nspan`
    specifying the number of spans or intervals to refine. The function uses a midpoint interpolation
    scheme to calculate the new refined points.

    Note:
    - The function modifies the `x_orig` array internally but returns a new array with the refined grid points.
    - The function requires the numpy module to be imported as `np`.

    Example:
    x = np.array([0.0, 1.0, 2.0, 3.0])
    refined_x = grid_refine_inner_edge(x, 3, 2)
    print(refined_x)
    Output: [0.    0.125 0.25  0.375 0.5   0.75  1.    1.5   2.    3.   ]
    """
    x = x_orig.copy()
    rev = x[0] > x[1]
    for ilev in range(nlev):
        x_new = 0.5 * (x[1:nspan + 1] + x[:nspan])
        x_ref = np.hstack((x, x_new))
        x_ref.sort()
        x = x_ref
        if rev:
            x = x[::-1]
    return x


def grid_radial_profile(rin, rout, nr, r_step_max=0.0):
    """
    Generates a radial grid based on the desired surface density profile.

    Args:
        rin (float): Inner radius of the radial grid.
        rout (float): Outer radius of the radial grid.
        nr (int): Number of grid points in the radial direction.
        r_step_max (float): The maximum step lengths for r grid, 
                            if the step lengths is larger in the logspace,
                            it will change to linear space to fullfill this
                            requirement.

    Returns:
        ndarray: The radial grid points.

    Description:
    This function generates a radial grid based on the desired surface density profile. The grid points are computed
    using logarithmically spaced values to capture the range of radial column densities. The grid is defined from the
    inner radius (`rin`) to the outer radius (`rout`) with `nr` number of grid points.
    """

    # collog = np.arange(19, 26, 0.1)  # radial column goes from 1e20-1e26
    # del_rim = 10 ** collog / 1.0e16  # approx n ~ 1e16 at inner ridge to get delta r
    # rinner = del_rim + rin  #
    # rmin = rinner[-1] + del_rim[-1]
    
    # simply use the log space for ri if we don't 
    # want to capture the outer rim
    rmin = rin
    ri = np.logspace(np.log10(rmin), np.log10(rout), nr + 1)
        
    if r_step_max > 0:
        ri_t_log = ri
        
        # Find the index where the difference exceeds r_step_max
        diffs = ri_t_log[1:] - ri_t_log[:-1]
        switch_idx = np.argmax(diffs > r_step_max)

        # If there's a point where the difference exceeds r_step_max, switch to linear
        if switch_idx > 0:
            # Create linear grid from the switch point to the end
            # linear_part = np.linspace(ri_t_log[switch_idx], rout, nr - switch_idx + 1)
            linear_part = np.arange(ri_t_log[switch_idx], rout + r_step_max, r_step_max)
            # Combine the logarithmic part and the linear part
            ri_t = np.concatenate((ri_t_log[:switch_idx + 1], linear_part[1:]))
        else:
            # If no switching is needed, just use the logarithmic grid
            ri_t = ri_t_log

        # Calculate the final number of points
        final_nr = len(ri_t) - 1

        # Print or store the final grid and number of points
        ri = ri_t
        nr = final_nr
        
    return ri, nr


def grid_vertical_layer(thetaup, ntheta, ncoarse):
    """
    Set up the grid in a vertical layer to correctly trace the temperature
    for the molecular vertical emitting layers.

    Args:
        thetaup(float): The upper limit of the refined layer in theta coordinates.
        ntheta(int): The number of cells in total
        ncoarse(int): The number of cells in the coarse regions.

    Returns:
        ndarray: The refined grid points in theta coordinates.

    Description:
    This function refines the grid in a vertical layer by generating a set of grid points in theta coordinates.
    The vertical layer is divided into three regions: coarse regions at the top and bottom, and a fine region in between.
    The number of cells in the coarse regions is specified by the `ncoarse` argument, and the grid points are generated
    accordingly using linearly spaced values.
    """
    ncoarse_bottom = int(0.25*ncoarse)
    ths = np.linspace(thetaup, 1.35, int(ncoarse-ncoarse_bottom))
    thh = np.linspace(1.351, 1.5, ntheta + 1 - ncoarse)
    thm = np.linspace(1.510, np.pi / 2.0, ncoarse_bottom)
    thetai = np.concatenate([ths, thh])
    thetai = np.concatenate([thetai, thm])

    return thetai


def place_ratio_g2d(ratio_g2d, shape=None, 
                    rc_cyl=None, zrc_cyl=None, 
                    rc=None, thetac=None, phic=None, 
                    nr_dense=None, ntheta_dense=None, 
                    rc_cyl_set = None,
                    method_polate1='linear', method_polate2='linear', 
                    coord='cyl', ratio_g2d_outgrid=100.):
    """
    **NOTE** (Developing)
    Places the g2d ratio values based on the provided input.

    Args:
        ratio_g2d (float or ndarray): The g2d ratio values to be placed.
        it can be either constant (float or int)
        or an array contains ratio_g2d[:, 0] as the radial grid (in cm)
        and the ratio_g2d[:, 1] as the corresponding gtd value at each radii.
        
        shape (list): the shape of the target coordinate.
        
        **NOTE**(Developing)
        rc_cyl_set (ndarry, optional): the cylindrical coordinate we want, this is used for sph coordinate, where we also need to do the conversion.
        
        TBA

    Returns:
        ndarray: The placed g2d ratio values.

    Description:
    This function places the g2d ratio values based on the provided input. The placement can be a constant value
    (float or int) or a function of radius. The function checks the type and shape of the input to determine the
    appropriate placement approach.

    If the g2d ratio is a constant value, the function prints a message indicating that the ratio is set up as a
    constant. If the g2d ratio is a function of radius, the function checks if the shape of the input matches the
    desired shape. If it matches, the function prints a message indicating that the pre-set function will be used.
    Otherwise, if the shape does not match or is not 2D, the function treats the input as a reference value and
    performs interpolation to generate the estimated ratio_g2d.
    """
    #  print(type(ratio_g2d))
    #  print(ratio_g2d.shape)
    
    if type(ratio_g2d) in [int, float, np.int32, np.int64, np.float32, np.float64]:
        print('the g2d ratio set up is a constant (float or float64 or integer), with ratio_g2d=%.2f' % (ratio_g2d))

        return ratio_g2d

    elif coord == 'cyl':
        print('sample the g2d to cylindrical coordinate acoording to cyl and the shape')

        if shape == None:
            nr, nz, nphi = len(rc_cyl), len(zrc_cyl), len(phic)
            shape = tuple(np.array([nr, ntheta, nphi]))

        if ratio_g2d.shape == shape:
            print('the shape of set up ratio_g2d matches with the sigmad,\nthus using the pre-set ratio g2d grid')

            return ratio_g2d

        else:
            assert len(ratio_g2d.shape) == 2,\
                'the shape of ratio_g2d does not match with the grid,'+\
                'nor a 2d reference function.\n'+\
                'It has to be a 3d grid that matches with the spatial grid in radmc3d (in spherical coordinate) or '+\
                'a reference 2d array (refer to the example).'
            print(
                'the set up ratio_g2d is a reference value, '+\
                'interpolate will be made to generate the estimated ratio_g2d at different spatial grid')

            # the reference ratio g2d is the one setup in the parameters.py
            # ratio_g2d_reference = np.loadtxt('Ratio_of_Luminosity.txt')
            ratio_g2d_reference = ratio_g2d.copy()

            # note that here the xdata and ydata is setup as the [:,0] and [:,1] in the dat
            x_data = ratio_g2d_reference[:, 0]  # df_read['radial[cm]'].values * au
            y_data = ratio_g2d_reference[:, 1]  # df_read['g2d'].values

            # extending the y value using the inner and outer edge of the y values
            x_fit = rc_cyl.copy()
            y_fit = np.interp(x_fit, x_data, y_data, left=None, right=None, period=None)
            y_fit[x_fit < 0.5 * (x_data[0] + x_data[1])] = y_data[0]
            y_fit[x_fit > np.max(x_data)] = y_data[-1]

            ratio_g2d_grid_cyl = np.ones(list(shape))
            # set up the whole vertical layer (same r in cyl coordinate) with the same gtd ratio
            for j_height in np.arange(shape[1]):
                ratio_g2d_grid_cyl[:, j_height, 0] = y_fit

            return ratio_g2d_grid_cyl

    else:
        # here coord == 'sph'
        print('the g2d ratio is setup as a function of radius and it needs to be sampled to spherical coordinate')

        if shape == None:
            nr, ntheta, nphi = len(rc), len(thetac), len(phic)
            shape = tuple(np.array([nr, ntheta, nphi]))
        else:
            nr, ntheta, nphi = shape

        if ratio_g2d.shape == shape:
            print('the shape of set up ratio_g2d matches with the sigmad,\nthus using the pre-set ratio g2d grid')

            return ratio_g2d

        else:
            assert len(ratio_g2d.shape) == 2, \
                'the shape of ratio_g2d does not match with the grid,'+\
                ' nor a 2d reference function.\nIt has to be a 3d grid that matches with the spatial grid in radmc3d'+\
                '(in spherical coordinate) or a reference 2d array (refer to the example).'
            print('the set up ratio_g2d is a reference value,'+\
                  'interpolate will be made to generate the estimated ratio_g2d at different spatial grid')

            # DENSE grid
            if nr_dense is None:
                nr_dense = 5 * nr
            if ntheta_dense is None:
                ntheta_dense = 5 * ntheta
            
            # the reference ratio g2d is the one setup in the parameters.py
            # ratio_g2d_reference = np.loadtxt('Ratio_of_Luminosity.txt')
            ratio_g2d_reference = ratio_g2d.copy()

            # note that here the xdata and ydata is setup as the [:,0] and [:,1] in the dat
            x_data = ratio_g2d_reference[:, 0]  # df_read['radial[AU]'].values * au
            y_data = ratio_g2d_reference[:, 1]  # df_read['g2d'].values

            # the mesh grid for the desired spherical coordinate
            qq = np.meshgrid(rc, thetac, phic, indexing='ij')
            rr = qq[0]  # Spherical R 
            tt = qq[1]  # Angle with z axis ( = pi/2 for the z=0 axis)
            # the cyl grid made to better sample the g2d
            
            #
            # New Grid for Sampling in the scipy griddata regridding
            # The new grid will exactly in r and z/r grids
            # 
            if rc_cyl_set is None: 
                rc_new_cyl = np.logspace(np.log10(0.95 * np.min(rc)), np.log10(1.05 * np.max(rc)), int(nr_dense / 3))
            else:
                rc_new_cyl = rc_cyl_set
            
            zrc_new_cyl = np.linspace(np.min(1 / np.tan(tt)), np.max(1 / np.tan(tt)), int(ntheta_dense / 3))

            x_fit = rc_new_cyl.copy()
            y_fit = np.interp(x_fit, x_data, y_data, left=None, right=None, period=None)
            y_fit[x_fit < 0.5 * (x_data[0] + x_data[1])] = y_data[0]
            y_fit[x_fit > np.max(x_data)] = y_data[-1]

            # put the ratio in the new cylindrical coordinate
            ratio_g2d_grid_cyl = np.ones([len(rc_new_cyl), len(zrc_new_cyl)])
            # set up the whole vertical layer (same r in cyl coordinate) with the same gtd ratio
            for j_height in np.arange(len(zrc_new_cyl)):
                ratio_g2d_grid_cyl[:, j_height] = y_fit

            """
            # now we get the ratio_g2d_grid in the cylindrical grid,
            # we then need to sample it to the spherical grid that the radmc3d requires.
            """
            
            qq_cyl = np.meshgrid(rc_new_cyl, zrc_new_cyl, phic, indexing='ij')
            rr_cyl = qq_cyl[0]
            zrr_cyl = qq_cyl[1]

            # flatten the coordinate
            R_flat_cyl = np.sqrt(rr_cyl.flatten() ** 2 + (zrr_cyl.flatten() * rr_cyl.flatten()) ** 2)
            tt_flat_cyl = np.pi / 2 - np.arctan(zrr_cyl.flatten())

            ratio_g2d_grid_flat_cyl = ratio_g2d_grid_cyl.flatten()

            #
            # DENSE grid for first sampling back
            #
            # Make Meshgrid
            rc_dense = rc.copy() # np.logspace(np.log10(np.min(rc)), np.log10(np.max(rc)), nr_dense)
            while len(rc_dense) <= nr_dense:
                rc_dense = np.hstack([rc_dense, 0.5 * (rc_dense[:-1] + rc_dense[1:])])
            rc_dense = np.sort(rc_dense)
            
            thetac_dense = thetac.copy()
            while len(thetac_dense) <= ntheta_dense:
                thetac_dense = np.hstack([thetac_dense, 0.5 * (thetac_dense[:-1] + thetac_dense[1:])])
            thetac_dense = np.sort(thetac_dense)

            qq_dense = np.meshgrid(rc_dense, thetac_dense, phic, indexing='ij')
            rr_dense = qq_dense[0]  # Spherical R 
            tt_dense = qq_dense[1]  # Angle with z axis ( = pi/2 for the z=0 axis)
            # zr_dense = np.pi/2.e0 - qq_dense[1]

            nr_dense = len(rc_dense)
            ntheta_dense = len(thetac_dense)

            # sampling to the DENSE spherical coordinate
            ratio_g2d_grid_sph_dense_new = np.zeros([nr_dense, ntheta_dense, nphi])
            ratio_g2d_grid_sph_dense_new = interpolate.griddata((np.log10(R_flat_cyl / const.au), tt_flat_cyl),
                                                                ratio_g2d_grid_flat_cyl,
                                                                (np.log10(rr_dense / const.au), tt_dense),
                                                                method=method_polate1)

            # Make Flat data for further sampling back to COARSE grid
            rr_flat_dense_new = rr_dense[:, :, 0].flatten()
            tt_flat_dense_new = tt_dense[:, :, 0].flatten()
            ratio_g2d_grid_sph_flat_dense_new = ratio_g2d_grid_sph_dense_new[:, :].flatten()

            #
            # New Grid for Sampling in the scipy griddata regridding
            # Now the grid is swiched back to the original grids
            # 
            qq = np.meshgrid(rc, thetac, phic, indexing='ij')
            rr = qq[0]  # Spherical R 
            tt = qq[1]  # Angle with z axis ( = pi/2 for the z=0 axis)
            zr = np.pi / 2.e0 - qq[1]
            zz = np.sin(zr) * rr

            #
            # Regridding to a New Grid
            # Adopt method of 'linear' because of undersampling
            #
            ratio_g2d_grid_sph = np.zeros([nr, ntheta, nphi])
            ratio_g2d_grid_sph = interpolate.griddata((np.log10(rr_flat_dense_new / const.au), tt_flat_dense_new),
                                                      ratio_g2d_grid_sph_flat_dense_new, (np.log10(rr / const.au), tt),
                                                      method=method_polate2)

            # get the output
            ratio_g2d = ratio_g2d_grid_sph.copy()
            # set the nan values (out of grid) points to the reference g2d
            # those regions don't matter because they should be very low density regions
            ratio_g2d[np.isnan(ratio_g2d)] = ratio_g2d_outgrid

            return ratio_g2d


def clean_rhodust_array(mint, rhodust):
    """
    Cleans a 4D dust density array by replacing non-finite values with 0.0 and 
    applying a minimum threshold for each dust species.

    Args:
        rhodust (ndarray): A 4D array of dust density values with shape (dim1, dim2, dim3, nrspec).
        mint (object): An object containing simulation parameters. Must have an attribute `fracs`
                       that provides the dust fractions for each species.

    Returns:
        ndarray: The cleaned 4D dust density array with the same shape as the input.
        
    Description:
        This function performs the following steps:
            1. Replaces any NaN or infinite values in the input array with 0.0.
            2. Computes the minimum dust density threshold for each dust species using the formula:
               rhodust_min = const.mu * 0.01 * mint.fracs
            3. Reshapes the computed thresholds for broadcasting across the first three dimensions.
            4. Applies the threshold across all dust species simultaneously using element-wise maximum.
    """
    # Copy the input array to avoid modifying it in place.
    rhodust_clean = rhodust.copy()

    # Clean up the array by replacing NaN and infinite values with 0.0
    rhodust_clean[~np.isfinite(rhodust_clean)] = 0.0

    # Calculate the minimum dust density for each dust species.
    # Ensure mint.fracs is converted to a NumPy array.
    rhodust_min = const.mu * 0.01 * np.array(mint.fracs)

    # Reshape the threshold array for broadcasting: shape becomes (1, 1, 1, nrspec)
    rhodust_min = rhodust_min.reshape(1, 1, 1, -1)

    # Apply the threshold by taking the element-wise maximum between the array and the threshold.
    rhodust_clean = np.maximum(rhodust_clean, rhodust_min)

    return rhodust_clean


def place_rhodust_all(mint, rhod_total_input,
                      rhodust_a_input = None, 
                      bool_dust_fragmentation=False,
                      bool_dust_radial_drifting=False,
                      tgas_input=[], 
                      verbose=True):
    """
    **NOTE** (Developing)
    
    Place the dust densities for each sub species.
    
    Args:
        mint (object): The DiskMINT model object containing the parameters and data.
        rhod_total_input (ndarray): the input rhodust for all.
                                    rhod_total = rhogas / ratio_g2d.
        rhodust_a_input (ndarray): the input rhodust for each size. 
                                   Default is none, and the rhodust_a will be set according to
                                   the rhod_total_input. 
        bool_dust_fragmentation (bool): if true, the rhodust will be distributed according to
                                 the size of each grains for each size bins.
                                 the fragmentation is based on the surface density of the gas,
                                 and depends on the fragmentation threshold (ufrag), 
                                 viscous alpha (alphat), and the grain material density (rhogr).
        bool_dust_radial_drifting (bool): if true, the rhodust will be drifted in. We move all the grains that need to drift (a > a_drift) inside the radius_drift (in the unit of au), and keep the total masses for each grains be exactly the same.
        tgas_input (ndarray): the gas temperature, used to relocate the dust grains per sizes
        verbose (bool): whether to print on screen for debugging
        
        TBA

    Returns:
        rhodust_a: The placed rhodust for each species.
        
    Descriptions:
        Currently we can set up two different dust surface density profiles, 
        and each has a sigma_i.
        
        Here we set up not only the sigma_i for each dust components, but also
        the subspecies (distribut to each dust size bins) rho_i_j.
    
    **NOTE**
    What is new:
        Adding the function of distributing the grains with different sizes to different sigmad,
        so now the sigmad would also be a function of size (a).
    """
    
    para = mint.parameters
    
    ############
    # this section below (defining sigmad_1 and sigmad_2) 
    # enables us to set up two different dust sigmad distributions
    # This is for TWO completely different dust species, so we can have two species with two 
    # different sigmad. This is useful when we want to add an inner disk for example.
    ############
    
    #
    # Set up the rhodust grid from the rhodust all
    #
    
    if rhodust_a_input is None:
        # set up the fraction grid for each dust species
        fracgrid_total_dust_1 = mint.sigmad_1/mint.sigmad
        fracgrid_total_dust_2 = mint.sigmad_2/mint.sigmad

        if verbose:
            print('-------------------------')
            index_t = (fracgrid_total_dust_1 + fracgrid_total_dust_2 - 1 < 1e-90)
            print('check: where the total of sigma is not 1:', mint.sigmad[index_t])

        # set up the fraction of the dust for each size bin
        fracs_sizebin_dust_1 = np.array(mint.fracs[0:para.nr_dust_1])
        fracs_sizebin_dust_2 = np.array(mint.fracs[para.nr_dust_1:para.nr_dust_1+para.nr_dust_2])
        if verbose:
            print('frac for dust1, dust2:', np.average(fracgrid_total_dust_1), np.average(fracgrid_total_dust_2))
            print('frac_dust_1:', fracs_sizebin_dust_1)
            print('frac_dust_2:', fracs_sizebin_dust_2)
            print('sum of frac1, frac2:', np.sum(fracs_sizebin_dust_1), np.sum(fracs_sizebin_dust_2))

        if verbose:
            print('-------------------------')
            index_t = (mint.rhod * fracgrid_total_dust_1 - mint.rhod_1 < 1e-90)
            print('check: where the rhod_1 is not correctly cast to fraction_1:', mint.sigmad[index_t])

        rhodust_a = np.zeros(list(mint.rhod.shape) + [para.dust_spec_nr])
        # mdust_orig_eachsize = np.zeros(list(mint.rhod.shape) + [para.dust_spec_nr])
        i_dust = -1
        for i_sizebin in np.arange(para.nr_dust_1):
            i_dust += 1
            rhodust_a[:, :, :, i_dust] = rhod_total_input * fracgrid_total_dust_1 * fracs_sizebin_dust_1[i_sizebin]
            # mdust_orig_eachsize[i_dust] = (rhodust_a[:, :, :, i_dust] * mint.vol).sum(0).sum(0).sum(0)

        for i_sizebin in np.arange(para.nr_dust_2):
            i_dust += 1
            rhodust_a[:, :, :, i_dust] = rhod_total_input * fracgrid_total_dust_2 * fracs_sizebin_dust_2[i_sizebin]
            # mdust_orig_eachsize[i_dust] = (rhodust_a[:, :, :, i_dust] * mint.vol).sum(0).sum(0).sum(0)
    else:
        # take the input rhodust_a
        rhodust_a = rhodust_a_input
        
    if bool_dust_fragmentation:
        """NOTE under development"""
        if verbose:
            print('start redistribute dust based on the fragmentation threshold')
        
        #-------------  NEED THESE INPUTS ------------------
        # these input values will be put in the parameter array
        ufrag   = para.vel_frag   # [cm/s] Fragmentation threshold (test value: 400.0, plan to adopt: 100.0)
        alphat  = para.visc_alpha # [unitless] Viscous alpha (test value: 0.001)
        rhogr   = para.rhobulk    # [g cm^-3] Grain material density being used (test value: 2.3)
        
        assert tgas_input.shape == rhod_total_input.shape, 'the input tgas shape does not match the input rhod shape, check the input!'
        
        # midplane temperature
        tgas_mid = np.ones(list(tgas_input.shape))
        for j_height in np.arange(mint.ntheta):
            tgas_mid[:, j_height, 0] = tgas_input[:, -1, 0]
                
        # nr = nr
        na = np.size(mint.a_ave) # = para.nr_dust_1 + para.nr_dust_2

        # Thermal speed squared
        cs2  = const.kk * tgas_mid / (const.mu) # H2 thermal speed

        # Max garain size at this radius
        # this is for defining what would be the max grain sizes at such radii
        # and then also define the radial profile of the dust emission extents
        afrag = (2.0/(3.0*np.pi)) * (ufrag**2.0/(alphat*cs2)) * (1.0/rhogr) * (mint.sigmagas_setup) # (mint.sigmad * mint.ratio_g2d)
        afrag[np.isnan(afrag)] = 1e-199
        afrag[afrag < 1e-199] = 1e-199

        # Redistribute mass in grains with size > afrag
        rhodust_a_orig = copy.deepcopy(rhodust_a)
        rhodust_a_t = copy.deepcopy(rhodust_a)
        for ir in range(mint.nr):
            for iz in range(mint.ntheta):
                adum = np.ones(na)
                for ia in range(na):
                    if mint.a_ave[ia] > afrag[ir, iz, 0]:
                        adum[ia] = 0
                
                # mass in dust grains below this value
                rhod_tot_below_afrag = np.sum(rhodust_a_t[ir,iz,0,:]*adum[:])
                # print(rhod_tot_below_afrag)
                
                if rhod_tot_below_afrag > 0.0:
                    mass_increase = rhod_total_input[ir,iz,0]/rhod_tot_below_afrag
                    # print(mass_increase)

                    # Re-normalize sigmas to add mass from size > afrag
                    rhodust_a_t[ir,iz,0,:] = rhodust_a_t[ir,iz,0,:] * adum[:] * mass_increase
        
        # set to the output   
        rhodust_a = rhodust_a_t
        
    if bool_dust_radial_drifting:
        """
        NOTE: recommend only use radial drift or fragmentation at a time for now
        (developing) use fragmentation and drifting together
        """
        if verbose:
            print('start redistribute dust based on the radial drift threshold')
              
        #-------------  NEED THESE INPUTS ------------------
        # these input values will be put in the parameter array
        a_drift = para.a_drift # cm. The grains that are larger than this size will be drifted in
        radius_drift = para.radius_drift # au. The grains that needs to be drifted will be dirfted inside this radius.
        
        # Redistribute mass in grains with size > a_drift into the radius_drift
        rhodust_a_orig = copy.deepcopy(rhodust_a)
        rhodust_a_t = copy.deepcopy(rhodust_a)
        
        # find the index for rc
        select_outside_radii = mint.rc > radius_drift
        select_inside_radii = mint.rc <= radius_drift
        
        # nr = nr
        na = np.size(mint.a_ave) # = para.nr_dust_1 + para.nr_dust_2
        # iterate through the disk
        mass_drift_in = 0.0
        for ia in range(na):
            # for each size, we drift the grains that are larger than a_drift
            # by putting all the masses of these grains ths is outside of radius_drift
            # inside. We add the masses equally to all the grains inside the radius_drift
            # because the grain collisional time scale is very short
            # and they will become equalibrium states in a short time
            if mint.a_ave[ia] >= a_drift:
                if verbose:
                    print('the %.2f cm dust grain drifted in'%(mint.a_ave[ia]))
                    print('this is dust grain #%i'%(ia))
                    
                # mass_increase = Mtotal / old_Minside 
                # 2 * (rhod_total_input * mint.vol).sum(0).sum(0).sum(0)
                # mass_increase = (rhodust_a_t[:, :, 0, ia] * mint.vol[:, :, 0]).sum(0).sum(0)/(rhodust_a_t[select_inside_radii, :, 0, ia] * mint.vol[select_inside_radii, :, 0]).sum(0).sum(0)
                
                # if verbose:
                #     print(mass_increase)
            
                mass_drift_in += (rhodust_a_t[select_outside_radii, :, 0, ia] * mint.vol[select_outside_radii, :, 0]).sum(0).sum(0)
                
                # set the smallest boundary of the rhodust
                # it is rhogas_minimum (1 molecule/cm^3) and gtd=100
                rhodust_a_t[select_outside_radii, :, 0, ia] = const.mu * 0.01 * mint.fracs[ia] # 1e-199
                
        if verbose:
            print('the mass drifted in is')
            print(mass_drift_in/const.ms, '[ms]')

        # sum up the mass inside, which was the previous mass inside before the drifting
        rhodust_all_t = rhodust_a_t.sum(3)
        mass_was_inside = (rhodust_all_t[select_inside_radii, :, 0] * mint.vol[select_inside_radii, :, 0]).sum(0).sum(0)
        
        # calculate the mass increase, this should be the value applied to all the grains
        mass_increase = mass_drift_in/mass_was_inside
        
        if verbose:
            print('the mass that was previously inside was')
            print(mass_was_inside/const.ms, '[ms]')
            print('the mass increase fraction is')
            print(mass_increase)
            
        rhodust_a_t[select_inside_radii, :, 0, :] = rhodust_a_t[select_inside_radii, :, 0, :] * (1 + mass_increase)
                
        # # redistribute these mass to the radii inside
        # rhodust_all_t = rhodust_a_t.sum(3)
        # total_mass_inside_for_alla = (rhodust_all_t[select_inside_radii, :, 0] * mint.vol[select_inside_radii, :, 0]).sum(0).sum(0)
        # for j_a in range(na):
        #     fraction_inside_ja_to_alla = \
        #         (rhodust_a_t[select_inside_radii, :, 0, j_a] * mint.vol[select_inside_radii, :, 0]).sum(0).sum(0)/total_mass_inside_for_alla
        #     rhodust_a_t[select_inside_radii, :, 0, j_a] = rhodust_a_t[select_inside_radii, :, 0, ia] * mass_increase * fraction_inside_ja_to_alla
            
        # rhodust_a_t[select_outside_radii, :, 0, ia] = 0.0
                    
        rhodust_a = rhodust_a_t
    
    if verbose:
        print('-------------------------')
        print('place dust density for each subspecies and sub sizes')
        
        print('the size of rhodust_all: ', rhodust_a.shape)
        
        mass_t = 2 * (rhod_total_input * mint.vol).sum(0).sum(0).sum(0)
        print('input dust disk mass = %.3e [ms]'%(mass_t/const.ms))

        if bool_dust_fragmentation:
            print('output disk mass before dust grain fragmentation')
            rhodust_a_total = rhodust_a_orig.sum(3)
            mass_t = 2 * (rhodust_a_total * mint.vol).sum(0).sum(0).sum(0)
            print('output dust disk mass = %.3e [ms]'%(mass_t/const.ms)) 
        
            print('-------')
            print('after dust grain relocating')
        
        rhodust_a_total = np.nansum(rhodust_a, axis=3)
        mass_t = 2 * (rhodust_a_total * mint.vol).sum(0).sum(0).sum(0)
        print('output dust disk mass = %.3e [ms]'%(mass_t/const.ms))

        rhodust_a_total = np.nansum(rhodust_a[:, :, :, 0:para.nr_dust_1], axis=3)
        mass_t = 2 * (rhodust_a_total * mint.vol).sum(0).sum(0).sum(0)
        print('output primary dust disk mass = %.3e [ms]'%(mass_t/const.ms))

        rhodust_a_total = np.nansum(rhodust_a[:, :, :, para.nr_dust_1:], axis=3)
        mass_t = 2 * (rhodust_a_total * mint.vol).sum(0).sum(0).sum(0)
        print('output secondary dust disk mass = %.3e [ms]'%(mass_t/const.ms))
    
    """
    do not convert the mass back because there might be some grains should have smaller mass now;
    we should do this properly above by checking their masses
    """
    # # make sure the mass is conservative to the dust mass that was set up
    # rhodust_a_total = rhodust_a.sum(3)
    # mass_t = 2 * (rhodust_a_total * mint.vol).sum(0).sum(0).sum(0)
    # rhodust_a = rhodust_a * (para.)
    
    """
    Put some minimum values for the rhodust to avoid numerical issues
    """
    rhodust_a = clean_rhodust_array(mint, rhodust_a)
    
    return rhodust_a


def place_surface_density(sigmad_input, shape=None, 
                          rc_cyl=None, zrc_cyl=None, 
                          rc=None, thetac=None, phic=None, 
                          nr_dense=500, ntheta_dense=500, 
                          method_polate1='linear', method_polate2='linear', 
                          coord='cyl'):
    """ 
    **NOTE** (Developing)
    Places the surface density values based on the provided input.

    **NOTE** call the sigma to sigmad to avoid confusion.
    
    Args:
        sigmad_input (ndarray): The sigmad_input. It needs to be in the format of 2d array where the first col is radius, and the 2nd col is sigmad.
        
        shape (list): the shape of the target coordinate.
        
        TBA

    Returns:
        ndarray: The placed surface density ratio values.

    Description:
    This function places the surface density distributions ratio values based on the provided input. 
    The placement needs to be a function of radius. It first checks the type and shape of the input to determine the appropriate placement approach. The function checks if the shape of the input matches the desired shape. If it matches, the function prints a message indicating that the pre-set function will be used. Otherwise, if the shape does not match or is not 2D, the function treats the input as a reference value and performs interpolation to generate the estimated surface density distribution (sigmad_output).
    """
    #  print(type(sigmad_input))
    #  print(sigmad_input.shape)

    if coord == 'cyl':
        print('sample the sigma to cylindrical coordinate acoording to cyl and the shape')

        if shape == None:
            nr, nz, nphi = len(rc_cyl), len(zrc_cyl), len(phic)
            shape = np.array([nr, nz, nphi])

        if sigmad_input.shape == shape:
            print('the shape of set up sigmad_input matches with the sigmad,\nthus using the pre-set ratio sigma grid')

            return sigmad_input

        else:
            # here coord == 'sph'
            assert len(sigmad_input.shape) == 2,\
                'the shape of sigmad_input does not match with the grid,'+\
                'nor a 2d reference function.\n'+\
                'It has to be a 3d grid that matches with the spatial grid in radmc3d (in spherical coordinate) or '+\
                'a reference 2d array (refer to the example).'
            
            print('the set up sigmad_input is a reference value, '+\
                'interpolate will be made to generate the estimated sigmad_input at different spatial grid')

            # the reference ratio sigma is the one setup in the parameters.py
            # sigmad_input_reference = np.loadtxt('Ratio_of_Luminosity.txt')
            sigmad_input_reference = sigmad_input.copy()

            # note that here the xdata and ydata is setup as the [:,0] and [:,1] in the dat
            x_data = sigmad_input_reference[:, 0]  # df_read['radial[cm]'].values * au
            y_data = sigmad_input_reference[:, 1]  # df_read['sigma'].values

            x_fit = rc_cyl.copy()
            y_fit = np.interp(x_fit, x_data, y_data, left=None, right=None, period=None)
            y_fit[x_fit < 0.5 * (x_data[0] + x_data[1])] = y_data[0]
            y_fit[x_fit > np.max(x_data)] = y_data[-1]

            sigmad_input_grid_cyl = np.ones(list(shape))
            # set up the whole vertical layer (same r in cyl coordinate) with the same gtd ratio
            for j_height in np.arange(shape[1]):
                sigmad_input_grid_cyl[:, j_height, 0] = y_fit

            return sigmad_input_grid_cyl

    else:
        print('the sigma ratio is setup as a function of radius and it needs to be sampled to spherical coordinate')

        if shape == None:
            nr, ntheta, nphi = len(rc), len(thetac), len(phic)
            shape = tuple(np.array([nr, ntheta, nphi]))
        else:
            nr, ntheta, nphi = shape
        
        # print(sigmad_input.shape, shape)
        # print(sigmad_input.shape == shape)
        
        if sigmad_input.shape == shape:
            print('the shape of set up sigmad_input matches with the sigmad,\nthus using the pre-set ratio sigma grid')

            return sigmad_input

        else:
            assert len(sigmad_input.shape) == 2, \
                'the shape of sigmad_input does not match with the grid,'+\
                ' nor a 2d reference function.\nIt has to be a 3d grid that matches with the spatial grid in radmc3d'+\
                '(in spherical coordinate) or a reference 2d array (refer to the example).'
            
            print('the set up sigmad_input is a reference value,'+\
                  'interpolate will be made to generate the estimated sigmad_input at different spatial grid')

            # the reference ratio sigma is the one setup in the parameters.py
            # sigmad_input_reference = np.loadtxt('Ratio_of_Luminosity.txt')
            sigmad_input_reference = sigmad_input.copy()

            # note that here the xdata and ydata is setup as the [:,0] and [:,1] in the dat
            x_data = sigmad_input_reference[:, 0]  # df_read['radial[AU]'].values * au
            y_data = sigmad_input_reference[:, 1]  # df_read['sigma'].values

            # the mesh grid for the desired spherical coordinate
            qq = np.meshgrid(rc, thetac, phic, indexing='ij')
            rr = qq[0]  # Spherical R 
            tt = qq[1]  # Angle with z axis ( = pi/2 for the z=0 axis)
            # the cyl grid made to better sample the sigma
            rc_new_cyl = np.logspace(np.log10(0.95 * np.min(rc)), np.log10(1.05 * np.max(rc)), int(nr_dense / 3))
            zrc_new_cyl = np.linspace(np.min(1 / np.tan(tt)), np.max(1 / np.tan(tt)), int(ntheta_dense / 3))

            x_fit = rc_new_cyl.copy()
            y_fit = np.interp(x_fit, x_data, y_data, left=None, right=None, period=None)
            y_fit[x_fit < 0.5 * (x_data[0] + x_data[1])] = y_data[0]
            y_fit[x_fit > np.max(x_data)] = y_data[-1]

            # put the ratio in the new cylindrical coordinate
            sigmad_input_grid_cyl = np.ones([len(rc_new_cyl), len(zrc_new_cyl)])
            # set up the whole vertical layer (same r in cyl coordinate) with the same gtd ratio
            for j_height in np.arange(len(zrc_new_cyl)):
                sigmad_input_grid_cyl[:, j_height] = y_fit

            # now we get the sigmad_input_grid in the cylindrical grid,
            # we then need to sample it to the spherical grid that the radmc3d requires.

            qq_cyl = np.meshgrid(rc_new_cyl, zrc_new_cyl, phic, indexing='ij')
            rr_cyl = qq_cyl[0]
            zrr_cyl = qq_cyl[1]

            # flatten the coordinate
            R_flat_cyl = np.sqrt(rr_cyl.flatten() ** 2 + (zrr_cyl.flatten() * rr_cyl.flatten()) ** 2)
            tt_flat_cyl = np.pi / 2 - np.arctan(zrr_cyl.flatten())

            sigmad_input_grid_flat_cyl = sigmad_input_grid_cyl.flatten()

            #
            # DENSE grid for first sampling back
            #
            # Make Meshgrid
            rc_dense = np.logspace(np.log10(np.min(rc)), np.log10(np.max(rc)), nr_dense)
            thetac_dense = thetac.copy()
            while len(thetac_dense) <= ntheta_dense:
                thetac_dense = np.hstack([thetac_dense, 0.5 * (thetac_dense[:-1] + thetac_dense[1:])])
            thetac_dense = np.sort(thetac_dense)

            qq_dense = np.meshgrid(rc_dense, thetac_dense, phic, indexing='ij')
            rr_dense = qq_dense[0]  # Spherical R 
            tt_dense = qq_dense[1]  # Angle with z axis ( = pi/2 for the z=0 axis)
            # zr_dense = np.pi/2.e0 - qq_dense[1]

            nr_dense = len(rc_dense)
            ntheta_dense = len(thetac_dense)

            # sampling to the DENSE spherical coordinate
            sigmad_input_grid_sph_dense_new = np.zeros([nr_dense, ntheta_dense, nphi])
            sigmad_input_grid_sph_dense_new = interpolate.griddata((np.log10(R_flat_cyl / const.au), tt_flat_cyl),
                                                                sigmad_input_grid_flat_cyl,
                                                                (np.log10(rr_dense / const.au), tt_dense),
                                                                method=method_polate1)

            # Make Flat data for further sampling back to COARSE grid
            rr_flat_dense_new = rr_dense[:, :, 0].flatten()
            tt_flat_dense_new = tt_dense[:, :, 0].flatten()
            sigmad_input_grid_sph_flat_dense_new = sigmad_input_grid_sph_dense_new[:, :].flatten()

            #
            # New Grid for Sampling in the scipy griddata regridding
            # Now the grid is swiched back to the original grids
            # 
            qq = np.meshgrid(rc, thetac, phic, indexing='ij')
            rr = qq[0]  # Spherical R 
            tt = qq[1]  # Angle with z axis ( = pi/2 for the z=0 axis)
            zr = np.pi / 2.e0 - qq[1]
            zz = np.sin(zr) * rr

            #
            # Regridding to a New Grid
            # Adopt method of 'linear' because of undersampling
            #
            sigmad_input_grid_sph = np.zeros([nr, ntheta, nphi])
            sigmad_input_grid_sph = interpolate.griddata((np.log10(rr_flat_dense_new / const.au), tt_flat_dense_new),
                                                      sigmad_input_grid_sph_flat_dense_new, (np.log10(rr / const.au), tt),
                                                      method=method_polate2)

            # get the output
            sigmad_output = sigmad_input_grid_sph.copy()
            # set the nan values (out of grid) points to the 1e-199 to avoid numerical errors
            # those regions don't matter because they should be very low density regions
            sigmad_output[np.isnan(sigmad_output)] = 1e-199

            return sigmad_output


########################
# the model procedures # 
########################

def problem_setup(mint):
    """
    Sets up the initial dust density distribution and the required files for running radmc3d.

    Args:
        mint: The DiskMINT model object containing the necessary parameters and data.

    Description:
    This function performs the initial setup for running radmc3d. It creates the necessary files and directories, including files for dust density, stars, grid, dust opacity, and the radmc3d control file. The function saves the gas density, writes the wavelength file, writes the stars.inp file with information about the central star, writes the grid file with the coordinates, writes the dust mass density file, and writes the dust opacity control file. Finally, it writes the radmc3d.inp control file with various simulation parameters.

    Note:
    - The function assumes the presence of the `os` module and the `np` module from numpy.
    - The function modifies the `model` object internally to set up the required files and directories.
    - The function utilizes the parameters and data stored in the `model` object.
    - The function saves all files to the mint.file_dir.
    """
    para = mint.parameters

    if not os.path.isdir(mint.file_dir):
        os.mkdir(mint.file_dir)

    ########################################
    # If Calculate Dust Settling is needed #
    ########################################
    #
    # Save rhogas
    #
    np.savetxt(os.path.join(mint.file_dir, 'rhogas_probset.inp'), mint.rhogas[:, :, 0])
    mint.rhogas_f = mint.rhogas.copy()

    ############################
    #### Write Input Files #####
    ############################

    #
    # Write the wavelength file
    #
    with open(os.path.join(mint.file_dir, 'wavelength_micron.inp'), 'w+') as f:
        f.write('%d\n' % (mint.nlam))
        for value in mint.lam:
            f.write('%13.6e\n' % (value))

    #
    # Write the stars.inp file
    #
    print('writing stars.inp, ')
    print('rstar: ', para.rstar / const.rs, ' rsolar')
    print('mstar: ', para.mstar / const.ms, ' msolar')
    # the position of the star is always set at the center
    pstar = np.array([0., 0., 0.])
    if para.fmodel_filename != '':
        with open(os.path.join(mint.file_dir, 'stars.inp'), 'w+') as f:
            f.write('2\n')
            f.write('1 %d\n\n' % (mint.nlam))
            f.write(
                '%13.13e %13.13e %13.13e %13.13e %13.13e\n\n' % (para.rstar, para.mstar, pstar[0], pstar[1], pstar[2]))
            for value in mint.lam:
                f.write('%13.6e\n' % (value))
            for value in mint.Fnu:
                f.write('%13.6e\n' % (value))
        #     f.write('\n%13.6e\n'%(-tstar))
    else:
        with open(os.path.join(mint.file_dir, 'stars.inp'), 'w+') as f:
            f.write('2\n')
            f.write('1 %d\n\n' % (mint.nlam))
            f.write(
                '%13.13e %13.13e %13.13e %13.13e %13.13e\n\n' % (para.rstar, para.mstar, pstar[0], pstar[1], pstar[2]))
            for value in mint.lam:
                f.write('%13.6e\n' % (value))
            f.write('\n%13.6e\n' % (-para.tstar))

    #
    # Write the grid file
    #
    with open(os.path.join(mint.file_dir, 'amr_grid.inp'), 'w+') as f:
        f.write('1\n')  # iformat
        f.write('0\n')  # AMR grid style  (0=regular grid, no AMR)
        f.write('100\n')  # Coordinate system: spherical
        f.write('0\n')  # gridinfo
        f.write('1 1 0\n')  # Include r,theta coordinates
        f.write('%d %d %d\n' % (mint.nr, mint.ntheta, 1))  # Size of grid
        for value in mint.ri:
            f.write('%13.13e\n' % (value))  # X coordinates (cell walls)
        for value in mint.thetai:
            f.write('%13.13e\n' % (value))  # Y coordinates (cell walls)
        for value in mint.phii:
            f.write('%13.13e\n' % (value))  # Z coordinates (cell walls)

    #
    # Write the dust mass density file
    #
    # fracs = np.loadtxt('fracs.inp')
    dust_nr = para.dust_spec_nr # len(mint.fracs)
    rhodust_input = mint.rhod.copy()
    # save rhodust input
    np.savetxt(os.path.join(mint.file_dir, 'rhodust_all_probset.inp'), mint.rhod[:, :, 0])
    # place the rhodust to different sizes
    rhodust_a = place_rhodust_all(mint, rhodust_input)
    # make a copy of the problem setup in my style
    rhodust_probset_flat = rhodust_a.ravel(order='F')
    with open('rhodust_probset.inp', 'w') as f:
        f.write('1\n')                      # Format number
        f.write(f'{mint.nr}\n')             # e.g., number of r
        f.write(f'{mint.ntheta}\n')         # e.g., number of theta
        f.write(f'{mint.nphi}\n')           # e.g., number of phi
        f.write(f'{dust_nr}\n')             # e.g., number of species or additional header info
        np.savetxt(f, rhodust_probset_flat, fmt='%13.6e')
    # save as the radmc3d format
    with open(os.path.join(mint.file_dir, 'dust_density.inp'), 'w+') as f:
        f.write('1\n')  # Format number
        f.write('%d\n' % (mint.nr * mint.ntheta * mint.nphi))  # Nr of cells
        f.write('%i\n' % (dust_nr))  # Nr of dust species
        for i in range(dust_nr):
            # data = mint.fracs[i] * mint.rhod.ravel(order='F')  # Create a 1-D view, fortran-style indexing
            data = rhodust_a[:,:,:,i].ravel(order='F')         # Create a 1-D view, fortran-style indexing
            data.tofile(f, sep='\n', format="%13.6e")
            f.write('\n')

    #
    # Dust opacity control file
    #
    """
    **NOTE**
    Jun-17, 2023
    add the option of setup multiple dust species instead of just two of them
    
    Apr-29, 2024
    but actually more than 2 dust species would quickly turn impractical in a sense
    Because each dust would need at least 10 dust size bins, then it adds up to the 
    total calculation time in RADMC3D and could quickly become impractical.
    
    So here, a combination of 2 sub-species make sense, and I would recommend try not 
    to use more than 2.
    """
    dustopacname_list = [para.dustopacname_1, para.dustopacname_2]
    with open(os.path.join(mint.file_dir, 'dustopac.inp'), 'w+') as f:
        f.write('2               Format number of this file\n')
        f.write('%i              Nr of dust species\n' % (dust_nr))
        f.write('============================================================================\n')
        for i_index, i_name in enumerate(dustopacname_list):
            if i_name == para.dustopacname_1:
                for j in range(para.nr_dust_1):
                    f.write('1               Way in which this dust species is read\n')
                    f.write('0               0=Thermal grain\n')
                    f.write('%s_%s           Extension of name of dustkappa_***.inp file\n' % (i_name, j))
                    f.write('----------------------------------------------------------------------------\n')
            elif i_name == para.dustopacname_2:
                for j in range(para.nr_dust_2):
                    f.write('1               Way in which this dust species is read\n')
                    f.write('0               0=Thermal grain\n')
                    f.write('%s_%s           Extension of name of dustkappa_***.inp file\n' % (i_name, j))
                    f.write('----------------------------------------------------------------------------\n')

    #
    # Write the radmc3d.inp control file
    #
    with open(os.path.join(mint.file_dir, 'radmc3d.inp'), 'w+') as f:
        f.write('nphot = %d\n' % (para.nphot))
        #                         f.write('scattering_mode_max = 1\n') # treat scattering in an isotropic way
        f.write('scattering_mode_max = 0\n')  # don't include scattering effect
        f.write('iranfreqmode = %d\n' % (para.scattering_mode))  # in example as 1; but 1 gives prob in midplane
        f.write('setthreads = %d\n' % (para.nthreads))  # 10
        f.write('modified_random_walk = 1\n')  # accelerate by MRW method
        f.write('istar_sphere = 1\n')  # recommended by output warning  
        # Cause central star size is more than 10% of inner grid radius.
        f.write('itempdecoup  = 1\n')  # Enable for different dust components to have different temperatures

    return 1


def func_Tgas(Tgas, nbin, Tdust, ndsd):
    """
    The function that would be used to calculate gas temperature from collisional thermal equilibrium
    based on the given gas temperature, dust temperature, and dust surface area distribution.

    Args:
        Tgas (float): Gas temperature.
        nbin (int): Number of bins in the dust size distribution.
        Tdust (ndarray): Array of dust temperatures as a 1D numpy array.
        ndsd (ndarray): Array of dust surface area distribution (n_d(a) pi a^2) values as a 1D numpy array.

    Returns:
        float: Average gas temperature.

    Description:
    This function describes the gas-dust thermal equilibrium.
    """
    func = 0.0
    for i in range(nbin):
        func += ndsd[i] * (Tgas - Tdust[i])
    return func


def func_Tgas_visc(Tgas, nbin, Tdust, ndsd, Omega, alpha_h=0.01):
    """
    **NOTE** This function is currently not used. Put here for reference.
    The function that would be used to calculate gas temperature from collisional thermal equilibrium
    this is the function that considers the viscous heating term

    Args:
        Tgas (float): Gas temperature.
        nbin (int): Number of bins in the dust size distribution.
        Tdust (ndarray): Array of dust temperatures as a 1D numpy array.
        ndsd (ndarray): Array of dust surface area distribution (n_d(a) pi a^2) values as a 1D numpy array.
        Omega (float): Angular velocity.
        alpha_h (float, optional): Viscosity parameter. Defaults to 0.01.

    Returns:
        float: gas temperature with viscous heating.

    Description:
    This function describes the gas-dust thermal equilibrium including considering viscous heating.

    Note:
    - The function requires the numpy module to be imported as `np`.
    - The default value of `alpha_h` is set to 0.01.
    """

    gamma_visc = 2.25 * alpha_h * Tgas * Omega
    v_H2 = np.sqrt(8 * const.kk * Tgas / (np.pi * const.mu))  # sound speed of H molecule
    # here we use the mu as the mean molecule as 3.8e-24 (overall mean molecular mass)
    func = gamma_visc / (0.3 * v_H2 * 2 * np.pi)  # add the normalized viscuos heating term
    #     print(func, v_H2, const.mu)
    for i in range(nbin):
        func = func - (ndsd[i] * (Tgas - Tdust[i]))

    return func


def func_ext(r, a, b):
    """
    A simple linear function based on the given parameters.

    Args:
        r (float): Input value.
        a (float): Coefficient for the linear term.
        b (float): Coefficient for the constant term.

    Returns:
        float: Result of the linear function.
    """
    return a * r + b


def extend_insphr(rc, thetac, phic, rhogas, Tgas):
    """
    Extends the inner and outer grids to avoid potential issues in spherical and cylindrical coordinate conversion.

    Args:
        rc (ndarray): Original radial grid points as a 1D numpy array.
        thetac (ndarray): Original polar angle grid points as a 1D numpy array.
        phic (ndarray): Original azimuthal angle grid points as a 1D numpy array.
        rhogas (ndarray): Original gas density distribution as a 3D numpy array.
        Tgas (ndarray): Original gas temperature distribution as a 3D numpy array.

    Returns:
        tuple: Tuple containing the extended gas density distribution as a 2D numpy array, the extended gas temperature distribution as a 2D numpy array, and the extended meshgrid.

    Description:
    This function extends the inner and outer grids of the original radial grid points (`rc`) to avoid potential issues during spherical and cylindrical coordinate conversion. It creates a meshgrid using the input grid points and then extends the `rc` grid by adding 10 inner and 10 outer grid points. The gas density distribution (`rhogas`) and gas temperature distribution (`Tgas`) are extended by linearly extrapolating the values from the original grid to the extended grid. The extended grid and the extrapolated gas density and temperature distributions are returned as a tuple.

    Note:
    - The function requires the numpy module to be imported as `np`.
    - The function assumes that the `curve_fit` function from the scipy.optimize module is available.
    - The function depends on the `func_ext` function, func_ext(r, a, b) = a*r + b.
    """
    #
    # Make the meshgrid
    #
    qq = np.meshgrid(rc, thetac, phic, indexing='ij')
    rr = qq[0]  # Spherical R 
    tt = qq[1]  # Angle with z axis ( = pi/2 for the z=0 axis)
    zr = np.pi / 2.e0 - qq[1]
    zz = np.sin(zr) * rr

    #
    # Extend the rc gird by adding 10 inner and 10 outer the origianl grids
    #
    rc_ext1 = np.logspace(np.log10(0.8 * np.min(np.min(rr) * np.sin(np.min(tt)))), np.log10(np.min(rc)), 10)
    rc_ext2 = np.logspace(np.log10(np.max(rc)), np.log10(1.2 * np.max(rr) / np.sin(np.min(tt))), 11)
    rc_ext2 = rc_ext2[1:]
    rc_ext_a = np.hstack([rc_ext1, rc, rc_ext2])
    #     print(len(rc_ext_a))

    qq_ext = np.meshgrid(rc_ext_a, thetac, phic, indexing='ij')

    #
    # extrapolate with a linear fit according to each radial column
    #
    # Set up the blanck array
    rhogas_sph_ext = np.zeros([len(rc_ext_a), len(thetac)])
    Tgas_sph_ext = np.zeros([len(rc_ext_a), len(thetac)])
    #  rhodust_sph_ext = np.zeros([len(rc_ext_a), len(thetac), 1, dust_spec_nr])

    # fill the middle part
    assert rhogas_sph_ext[10:-10, :].shape == rhogas[:, :, 0].shape, 'rhogas shape seems wrong [%s]' % (
        str(rhogas[:, :, 0].shape))
    assert Tgas_sph_ext[10:-10, :].shape == Tgas[:, :, 0].shape, 'Tgas shape seems wrong [%s]' % (
        str(Tgas[:, :, 0].shape))

    rhogas_sph_ext[10:-10, :] = rhogas[:, :, 0]
    Tgas_sph_ext[10:-10, :] = Tgas[:, :, 0]

    # Ext rhogas
    print('extending grid for rhogas')
    for i_t in range(len(thetac)):
        xdata = rr[:, i_t, 0]
        ydata = rhogas[:, i_t, 0]

        xdata = xdata[np.isnan(ydata) == 0]
        ydata = ydata[np.isnan(ydata) == 0]

        # extended inner
        xdata_inn = xdata[:10]
        ydata_inn = ydata[:10]

        popt, pcov = curve_fit(func_ext, np.log10(xdata_inn), np.log10(ydata_inn))

        xfit = rc_ext1
        lg10_yfit = func_ext(np.log10(xfit), *popt)
        yfit = 10 ** lg10_yfit

        rhogas_sph_ext[:10, i_t] = yfit

        # extended outer
        xdata_out = xdata[-10:]
        ydata_out = ydata[-10:]

        popt, pcov = curve_fit(func_ext, np.log10(xdata_out), np.log10(ydata_out))

        xfit = rc_ext2
        lg10_yfit = func_ext(np.log10(xfit), *popt)
        yfit = 10 ** lg10_yfit

        rhogas_sph_ext[-10:, i_t] = yfit

    # Ext Tgas
    print('extending grid for Tgas')
    for i_t in range(len(thetac)):
        xdata = rr[:, i_t, 0]
        ydata = Tgas[:, i_t, 0]

        xdata = xdata[np.isnan(ydata) == 0]
        ydata = ydata[np.isnan(ydata) == 0]

        # extended inner
        xdata_inn = xdata[:10]
        ydata_inn = ydata[:10]

        popt, pcov = curve_fit(func_ext, np.log10(xdata_inn), np.log10(ydata_inn))

        xfit = rc_ext1
        lg10_yfit = func_ext(np.log10(xfit), *popt)
        yfit = 10 ** lg10_yfit

        Tgas_sph_ext[:10, i_t] = yfit

        # extended outer
        xdata_out = xdata[-10:]
        ydata_out = ydata[-10:]

        popt, pcov = curve_fit(func_ext, np.log10(xdata_out), np.log10(ydata_out))

        xfit = rc_ext2
        lg10_yfit = func_ext(np.log10(xfit), *popt)
        yfit = 10 ** lg10_yfit

        Tgas_sph_ext[-10:, i_t] = yfit

    return rhogas_sph_ext, Tgas_sph_ext, qq_ext


def sph2cyl(rc, thetac, phic, rhogas, Tgas, 
            nr_cyl=500, nz_cyl=500, 
            nr_dense=None, ntheta_dense=None,
            rc_cyl_set = None, 
            method_polate1='linear', method_polate2='linear'):
    """
    Performs the conversion from spherical coordinates to cylindrical coordinates.

    Args:
        rc (ndarray): Original radial grid points as a 1D numpy array.
        thetac (ndarray): Original polar angle grid points as a 1D numpy array.
        phic (ndarray): Original azimuthal angle grid points as a 1D numpy array.
        rhogas (ndarray): Original gas density distribution as a 3D numpy array.
        Tgas (ndarray): Original gas temperature distribution as a 3D numpy array.
        nr_cyl (int, optional): Number of radial grid points in the cylindrical coordinate. Defaults to 500.
        nz_cyl (int, optional): Number of axial grid points in the cylindrical coordinate. Defaults to 500.
        nr_dense (int, optional): Number of radial grid points in the dense spherical coordinate. Defaults to 500.
        ntheta_dense (int, optional): Number of polar angle grid points in the dense spherical coordinate. Defaults to 500.
        method_polate1 (str, optional): Interpolation method for extending the grid in spherical coordinates. Defaults to 'linear'.
        method_polate2 (str, optional): Interpolation method for regridding in cylindrical coordinates. Defaults to 'nearest'. (changed to linear on 01/28)

        **NOTE**(Developing)
        rc_cyl_set (ndarry, optional): the cylindrical coordinate we want
        
    Returns:
        tuple: Tuple containing the new radial grid points in cylindrical coordinates as a 1D numpy array, the new axial grid points in cylindrical coordinates as a 1D numpy array, the new meshgrid in cylindrical coordinates, the extended and regridded gas density distribution as a 3D numpy array, and the extended and regridded gas temperature distribution as a 3D numpy array.

    Description:
    This function converts the input gas density distribution (`rhogas`) and gas temperature distribution (`Tgas`) from spherical coordinates to cylindrical coordinates. First, it extends and samples the grids in the spherical coordinates by performing interpolation based on the specified interpolation method (`method_polate1` default `linear`). This step samples the grid from a COARSE spherical to a DENSE spherical coordinate. Then, it converts the DENSE grids to cylindrical coordinates, and it regrids the gas density and temperature distributions from the DENSE spherical grids to the new cylindrical grids using the specified interpolation method (`method_polate2`).

    Note:
    - The function requires the numpy module to be imported as `np`.
    - The function assumes that the `extend_insphr` function and the `interpolate.griddata` function from the scipy.interpolate module are available.
    - The constants (for example `const.au`) need to be defined before calling this function.
    """

    rhogas_sph_ext, Tgas_sph_ext, qq_ext = extend_insphr(rc, thetac, phic, rhogas, Tgas)

    #
    # Make the meshgrid
    # 
    # For original grid
    qq = np.meshgrid(rc, thetac, phic, indexing='ij')
    rr = qq[0]  # Spherical R 
    tt = qq[1]  # Angle with z axis ( = pi/2 for the z=0 axis)
    zr = np.pi / 2.e0 - qq[1]
    zz = np.sin(zr) * rr
    nphi = len(phic)
    # For extended grid
    rr_ext = qq_ext[0]  # Spherical R 
    tt_ext = qq_ext[1]  # Angle with z axis ( = pi/2 for the z=0 axis)
    zr_ext = np.pi / 2.e0 - qq_ext[1]
    zz_ext = np.sin(zr_ext) * rr_ext
    #
    # Name all the data in the extended grid
    #
    rr_flat_ext_old = rr_ext[:, :, 0].flatten()
    tt_flat_ext_old = tt_ext[:, :, 0].flatten()
    # zr_flat_ext_old = z_flat_ext_old/x_flat_ext_old
    Tgas_flat_ext_old = Tgas_sph_ext[:, :].flatten()
    rhogas_flat_ext_old = rhogas_sph_ext[:, :].flatten()

    #
    # DENSE grid
    #
    if nr_dense is None:
        nr_dense = nr_cyl
    if ntheta_dense is None:
        ntheta_dense = nz_cyl
        
    # Make Meshgrid
    # rc_dense = np.logspace(np.log10(np.min(rr_ext)), np.log10(np.max(rr_ext)), nr_dense)
    rc_dense = rr_ext[:, -1, 0].copy()
    while len(rc_dense) <= nr_dense:
        rc_dense = np.hstack([rc_dense, 0.5 * (rc_dense[:-1] + rc_dense[1:])])
    rc_dense = np.sort(rc_dense)
    
    thetac_dense = thetac.copy()
    while len(thetac_dense) <= ntheta_dense:
        thetac_dense = np.hstack([thetac_dense, 0.5 * (thetac_dense[:-1] + thetac_dense[1:])])
    thetac_dense = np.sort(thetac_dense)

    qq_dense = np.meshgrid(rc_dense, thetac_dense, phic, indexing='ij')
    rr_dense = qq_dense[0]  # Spherical R 
    tt_dense = qq_dense[1]  # Angle with z axis ( = pi/2 for the z=0 axis)
    # zr_dense = np.pi/2.e0 - qq_dense[1]

    nr_dense = len(rc_dense)
    ntheta_dense = len(thetac_dense)

    print('the DENSE grid: ', nr_dense, ntheta_dense)

    #
    # linear interpolate to setup all the cells in spherical coordinate
    #
    grid_old = (np.log10(rr_flat_ext_old / const.au), tt_flat_ext_old)
    grid_new = (np.log10(rr_dense / const.au), tt_dense)

    # Tgas_sph_dense_old = np.zeros([nr_dense, ntheta_dense, nphi])
    logvalue = np.log10(Tgas_flat_ext_old)
    grid_dense_Tg = interpolate.griddata(grid_old, logvalue, grid_new, method=method_polate1)
    Tgas_sph_dense_old = 10 ** grid_dense_Tg.copy()

    # rhogas_sph_dense_old = np.zeros([nr_dense, ntheta_dense, nphi])
    logvalue = np.log10(rhogas_flat_ext_old)
    grid_sph_rg = interpolate.griddata(grid_old, logvalue, grid_new, method=method_polate1)
    rhogas_sph_dense_old = 10 ** grid_sph_rg.copy()

    #
    # Name all the data in extended DENSE grid
    #
    r_flat_cyl_old = rr_dense[:, :, 0].flatten() * np.sin(tt_dense[:, :, 0].flatten())
    z_flat_cyl_old = rr_dense[:, :, 0].flatten() * np.cos(tt_dense[:, :, 0].flatten())
    zr_flat_cyl_old = z_flat_cyl_old / r_flat_cyl_old
    Tgas_flat_sph_dense_old = Tgas_sph_dense_old[:, :].flatten()
    rhogas_flat_sph_dense_old = rhogas_sph_dense_old[:, :].flatten()

    #
    # New Grid for Sampling in the scipy griddata regridding
    # The new grid will exactly in r and z/r grids
    # 
    if rc_cyl_set is None: 
        rc_new = np.logspace(np.log10(0.95 * np.min(rc)), np.log10(1.05 * np.max(rc)), nr_cyl)
    else:
        rc_new = rc_cyl_set
        
    zrc_new = np.linspace(np.min(1 / np.tan(tt)), np.max(1 / np.tan(tt)), nz_cyl)
    zrc_new = zrc_new[::-1]
    # zrc_new = 1/np.tan(tt)
    qq_new = np.meshgrid(rc_new, zrc_new, phic, indexing='ij')
    rq_new = qq_new[0]
    zrq_new = qq_new[1]

    #
    # Regridding to a New Grid
    # The Method is 'cubic' because it is an Oversampling
    # But the 'cubic' method will bring nan value and further destroy the results
    # so now adopt 'linear'
    # method_polate = 'cubic' # now it is defined in the input
    grid_old = (np.log10(r_flat_cyl_old / const.au), zr_flat_cyl_old)
    grid_new = (np.log10(rq_new / const.au), zrq_new)

    # Tgas_cyl_old = np.zeros([nr_cyl, nz_cyl, nphi])
    logvalue = np.log10(Tgas_flat_sph_dense_old)
    grid_rz_Tg = interpolate.griddata(grid_old, logvalue, grid_new, method=method_polate2)
    Tgas_cyl_old = 10 ** grid_rz_Tg.copy()

    # rhogas_cyl_old = np.zeros([nr_cyl, nz_cyl, nphi])
    logvalue = np.log10(rhogas_flat_sph_dense_old)
    grid_rz_rg = interpolate.griddata(grid_old, logvalue, grid_new, method=method_polate2)
    rhogas_cyl_old = 10 ** grid_rz_rg.copy()
    
    # there are a few grid point at the edge might not be sampled (at the top/bottom in the outer disk)
    Tgas_cyl_old[:, :, 0] = fill_nan_in_2d_array_inlog(Tgas_cyl_old[:, :, 0])
    rhogas_cyl_old[:, :, 0] = fill_nan_in_2d_array_inlog(rhogas_cyl_old[:, :, 0])

    return rc_new, zrc_new, qq_new, rhogas_cyl_old, Tgas_cyl_old


def cyl2sph(rc, thetac, phic, rr_cyl, zrr_cyl, rhogas_cyl_new, Tgas_cyl_old, nr_dense=None, ntheta_dense=None,
            method_polate1='linear', method_polate2='linear'):
    """
    Performs the conversion from cylindrical coordinates back to spherical coordinates.

    Args:
        rc (ndarray): Original radial grid points as a 1D numpy array.
        thetac (ndarray): Original polar angle grid points as a 1D numpy array.
        phic (ndarray): Original azimuthal angle grid points as a 1D numpy array.
        rr_cyl (ndarray): Radial grid points in the cylindrical coordinate as a 3D numpy array.
        zrr_cyl (ndarray): Vertical grid points (z/r ratio) in the cylindrical coordinate as a 3D numpy array.
        rhogas_cyl_new (ndarray): Gas density distribution in the cylindrical coordinate as a 3D numpy array.
        Tgas_cyl_old (ndarray): Gas temperature distribution in the cylindrical coordinate as a 3D numpy array.
        nr_dense (int, optional): Number of radial grid points in the dense spherical coordinate. Defaults to 500.
        ntheta_dense (int, optional): Number of polar angle grid points in the dense spherical coordinate. Defaults to 500.
        method_polate1 (str, optional): Interpolation method for sampling back to the dense spherical coordinate. Defaults to 'nearest'. (changed to linear on 01/28)
        method_polate2 (str, optional): Interpolation method for regridding to the original spherical coordinate. Defaults to 'linear'.

    Returns:
        tuple: Tuple containing the regridded gas density distribution in the spherical coordinate as a 3D numpy array and the regridded gas temperature distribution in the spherical coordinate as a 3D numpy array.

    Description:
    This function converts the input gas density distribution (`rhogas_cyl_new`) and gas temperature distribution (`Tgas_cyl_old`) from cylindrical coordinates back to spherical coordinates. It first samples the input distributions to a DENSE spherical coordinate using the specified interpolation method (`method_polate1`). Then, it regrids the DENSE spherical grids to the original (COARSE) spherical grids using the specified interpolation method (`method_polate2`).

    Note:
    - The function requires the numpy module to be imported as `np`.
    - The function assumes that the `interpolate.griddata` function from the scipy.interpolate module is available.
    - The constants (for example `const.au`) need to be defined before calling this function.
    """

    nr, ntheta, nphi = len(rc), len(thetac), len(phic)

    R_flat_cyl = np.sqrt(rr_cyl.flatten() ** 2 + (zrr_cyl.flatten() * rr_cyl.flatten()) ** 2)
    tt_flat_cyl = np.pi / 2 - np.arctan(zrr_cyl.flatten())

    rhogas_flat_cyl_new = rhogas_cyl_new.flatten()
    Tgas_flat_cyl = Tgas_cyl_old.flatten()

    #
    # DENSE grid for first sampling back
    #
    if nr_dense is None:
        nr_dense = rr_cyl.shape[0]
    if ntheta_dense is None:
        ntheta_dense = rr_cyl.shape[1]
        
    # Make Meshgrid
    rc_dense = rc.copy() # np.logspace(np.log10(np.min(rc)), np.log10(np.max(rc)), nr_dense)
    while len(rc_dense) <= nr_dense:
        rc_dense = np.hstack([rc_dense, 0.5 * (rc_dense[:-1] + rc_dense[1:])])
    rc_dense = np.sort(rc_dense)
    
    thetac_dense = thetac.copy()
    while len(thetac_dense) <= ntheta_dense:
        thetac_dense = np.hstack([thetac_dense, 0.5 * (thetac_dense[:-1] + thetac_dense[1:])])
    thetac_dense = np.sort(thetac_dense)

    qq_dense = np.meshgrid(rc_dense, thetac_dense, phic, indexing='ij')
    rr_dense = qq_dense[0]  # Spherical R 
    tt_dense = qq_dense[1]  # Angle with z axis ( = pi/2 for the z=0 axis)
    # zr_dense = np.pi/2.e0 - qq_dense[1]

    nr_dense = len(rc_dense)
    ntheta_dense = len(thetac_dense)

    # sampling to the DENSE spherical coordinate
    grid_old = (np.log10(R_flat_cyl / const.au), tt_flat_cyl)
    grid_new = (np.log10(rr_dense / const.au), tt_dense)
    
    # Tgas_sph_dense_new = np.zeros([nr_dense, ntheta_dense, nphi])
    logvalue = np.log10(Tgas_flat_cyl)
    grid_dense_Tg = interpolate.griddata(grid_old, logvalue, grid_new, method=method_polate1)
    Tgas_sph_dense_new = 10 ** grid_dense_Tg.copy()

    # rhogas_sph_dense_new = np.zeros([nr_dense, ntheta_dense, nphi])
    # rhogas_flat_cyl_new[rhogas_flat_cyl_new <= const.mu] = const.mu # assume there is 1 molecule in 1 cm^-3 in the ISM # this is done for the input rhogas now.
    logvalue = np.log10(rhogas_flat_cyl_new)
    grid_sph_rg = interpolate.griddata(grid_old, logvalue, grid_new, method=method_polate1)
    rhogas_sph_dense_new = 10 ** grid_sph_rg.copy()

    # clean out the NaNs at the edge
    # and assign those edge grids with the density in ISM
    rhogas_sph_dense_new[np.isnan(rhogas_sph_dense_new)] = 0.0
    rhogas_sph_dense_new[rhogas_sph_dense_new <= const.mu] = const.mu
    
    # Make Flat data for further sampling back to COARSE grid
    rr_flat_dense_new = rr_dense[:, :, 0].flatten()
    tt_flat_dense_new = tt_dense[:, :, 0].flatten()
    # zr_flat_ext_old = z_flat_ext_old/x_flat_ext_old
    Tgas_flat_dense_new = Tgas_sph_dense_new[:, :].flatten()
    rhogas_flat_dense_new = rhogas_sph_dense_new[:, :].flatten()

    #
    # New Grid for Sampling in the scipy griddata regridding
    # Now the grid is swiched back to the original grids
    # 
    qq = np.meshgrid(rc, thetac, phic, indexing='ij')
    rr = qq[0]  # Spherical R 
    tt = qq[1]  # Angle with z axis ( = pi/2 for the z=0 axis)
    zr = np.pi / 2.e0 - qq[1]
    zz = np.sin(zr) * rr

    #
    # Regridding to a New Grid
    # Adopt method of 'linear' because of undersampling
    # 
    grid_old = (np.log10(rr_flat_dense_new / const.au), tt_flat_dense_new)
    grid_new = (np.log10(rr / const.au), tt)

    # Tg_sphr_new = np.zeros([nr, ntheta, nphi])
    logvalue = np.log10(Tgas_flat_dense_new)
    grid_rz_Tg = interpolate.griddata(grid_old, logvalue, grid_new, method=method_polate2)
    Tg_sphr_new = 10 ** grid_rz_Tg.copy()

    # rhogas_sphr_new = np.zeros([nr, ntheta, nphi])
    logvalue = np.log10(rhogas_flat_dense_new)
    grid_rz_rg = interpolate.griddata(grid_old, logvalue, grid_new, method=method_polate2)
    rhogas_sphr_new = 10 ** grid_rz_rg.copy()
    
    # there are a few grid point at the edge might not be sampled (at the top/bottom in the outer disk)
    Tg_sphr_new[:, :, 0] = fill_nan_in_2d_array_inlog(Tg_sphr_new[:, :, 0])
    rhogas_sphr_new[:, :, 0] = fill_nan_in_2d_array_inlog(rhogas_sphr_new[:, :, 0])

    return rhogas_sphr_new, Tg_sphr_new


def extend_insphr_dust(rc, thetac, phic, rhodust):
    """
    Extends the inner and outer grids to avoid potential issues in spherical and cylindrical coordinate conversion.

    Args:
        rc (ndarray): Original radial grid points as a 1D numpy array.
        thetac (ndarray): Original polar angle grid points as a 1D numpy array.
        phic (ndarray): Original azimuthal angle grid points as a 1D numpy array.
        rhodust (ndarray): Original gas density distribution as a 3D numpy array.
        Tgas (ndarray): Original gas temperature distribution as a 3D numpy array.

    Returns:
        tuple: Tuple containing the extended gas density distribution as a 2D numpy array, the extended gas temperature distribution as a 2D numpy array, and the extended meshgrid.

    Description:
    This function extends the inner and outer grids of the original radial grid points (`rc`) to avoid potential issues during spherical and cylindrical coordinate conversion. It creates a meshgrid using the input grid points and then extends the `rc` grid by adding 10 inner and 10 outer grid points. The dust density distribution (`rhodust`) are extended by linearly extrapolating the values from the original grid to the extended grid. The extended grid and the extrapolated gas density and temperature distributions are returned as a tuple.

    Note:
    - The function requires the numpy module to be imported as `np`.
    - The function assumes that the `curve_fit` function from the scipy.optimize module is available.
    - The function depends on the `func_ext` function, func_ext(r, a, b) = a*r + b.
    """
    #
    # Make the meshgrid
    #
    qq = np.meshgrid(rc, thetac, phic, indexing='ij')
    rr = qq[0]  # Spherical R 
    tt = qq[1]  # Angle with z axis ( = pi/2 for the z=0 axis)
    zr = np.pi / 2.e0 - qq[1]
    zz = np.sin(zr) * rr

    #
    # Extend the rc gird by adding 10 inner and 10 outer the origianl grids
    #
    rc_ext1 = np.logspace(np.log10(0.8 * np.min(np.min(rr) * np.sin(np.min(tt)))), np.log10(np.min(rc)), 10)
    rc_ext2 = np.logspace(np.log10(np.max(rc)), np.log10(1.2 * np.max(rr) / np.sin(np.min(tt))), 11)
    rc_ext2 = rc_ext2[1:]
    rc_ext_a = np.hstack([rc_ext1, rc, rc_ext2])
    #     print(len(rc_ext_a))

    qq_ext = np.meshgrid(rc_ext_a, thetac, phic, indexing='ij')

    #
    # extrapolate with a linear fit according to each radial column
    #
    # Set up the blanck array
    nr_dust_spec = int(rhodust.shape[3])
    rhodust_sph_ext = np.zeros([len(rc_ext_a), len(thetac), 1, nr_dust_spec])
    print('the shape of rhodust_sph_ext %s'%(str(rhodust_sph_ext.shape)))
    
    # fill the middle part
    for i_dust in range(nr_dust_spec):
        assert rhodust_sph_ext[10:-10, :, 0, i_dust].shape == rhodust[:, :, 0, i_dust].shape, 'rhodust shape seems wrong [extended %s] vs. [original %s] for species [%s]' % (
            str(rhodust_sph_ext[10:-10, :, 0, i_dust].shape),
            str(rhodust[:, :, 0, i_dust].shape),
            str(i_dust),
            )

    rhodust_sph_ext[10:-10, :, 0, :] = rhodust[:, :, 0, :]

    # Ext rhodust
    for i_dust in range(nr_dust_spec):
        for i_t in range(len(thetac)):
            xdata = rr[:, i_t, 0]
            ydata = rhodust[:, i_t, 0, i_dust]
            ydata[ydata < 1.e-300] = 1.e-300
            
            xdata = xdata[np.isnan(ydata) == 0]
            ydata = ydata[np.isnan(ydata) == 0]

            # extended inner
            xdata_inn = xdata[:10]
            ydata_inn = ydata[:10]
            
            popt, pcov = curve_fit(func_ext, np.log10(xdata_inn), np.log10(ydata_inn))

            xfit = rc_ext1
            lg10_yfit = func_ext(np.log10(xfit), *popt)
            yfit = 10 ** lg10_yfit

            rhodust_sph_ext[:10, i_t, 0, i_dust] = yfit

            # extended outer
            xdata_out = xdata[-10:]
            ydata_out = ydata[-10:]

            popt, pcov = curve_fit(func_ext, np.log10(xdata_out), np.log10(ydata_out))

            xfit = rc_ext2
            lg10_yfit = func_ext(np.log10(xfit), *popt)
            yfit = 10 ** lg10_yfit

            rhodust_sph_ext[-10:, i_t, 0, i_dust] = yfit

    return rhodust_sph_ext, qq_ext


def sph2cyl_dust(rc, thetac, phic, rhodust, 
                 nr_cyl=500, nz_cyl=500, 
                 nr_dense=None, ntheta_dense=None, 
                 rc_cyl_set = None, 
                 method_polate1='linear', method_polate2='linear'):
    """
    Performs the conversion from spherical coordinates to cylindrical coordinates.

    Args:
        rc (ndarray): Original radial grid points as a 1D numpy array.
        thetac (ndarray): Original polar angle grid points as a 1D numpy array.
        phic (ndarray): Original azimuthal angle grid points as a 1D numpy array.
        rhodust (ndarray): Original dust density distribution as a 3D+1 numpy array.
        nr_cyl (int, optional): Number of radial grid points in the cylindrical coordinate. Defaults to 500.
        nz_cyl (int, optional): Number of axial grid points in the cylindrical coordinate. Defaults to 500.
        nr_dense (int, optional): Number of radial grid points in the dense spherical coordinate. Defaults to 500.
        ntheta_dense (int, optional): Number of polar angle grid points in the dense spherical coordinate. Defaults to 500.
        method_polate1 (str, optional): Interpolation method for extending the grid in spherical coordinates. Defaults to 'linear'.
        method_polate2 (str, optional): Interpolation method for regridding in cylindrical coordinates. Defaults to 'nearest'. (changed to linear on 01/28)
        
        **NOTE**(Developing)
        rc_cyl_set (ndarry, optional): the cylindrical coordinate we want
        
    Returns:
        tuple: Tuple containing the new radial grid points in cylindrical coordinates as a 1D numpy array, the new axial grid points in cylindrical coordinates as a 1D numpy array, the new meshgrid in cylindrical coordinates, the extended and regridded gas density distribution as a 3D numpy array, and the extended and regridded gas temperature distribution as a 3D numpy array.

    Description:
    This function converts the input dust density distribution (`rhodust`) from spherical coordinates to cylindrical coordinates. First, it extends and samples the grids in the spherical coordinates by performing interpolation based on the specified interpolation method (`method_polate1` default `linear`). This step samples the grid from a COARSE spherical to a DENSE spherical coordinate. Then, it converts the DENSE grids to cylindrical coordinates, and it regrids the gas density and temperature distributions from the DENSE spherical grids to the new cylindrical grids using the specified interpolation method (`method_polate2`).

    Note:
    - The function requires the numpy module to be imported as `np`.
    - The function assumes that the `extend_insphr_dust` function and the `interpolate.griddata` function from the scipy.interpolate module are available.
    - The constants (for example `const.au`) need to be defined before calling this function.
    """
    #
    # extend the rhodust first inorder to cover enough r grid in cylindrical coordinate
    #
    rhodust_sph_ext, qq_ext = extend_insphr_dust(rc, thetac, phic, rhodust)

    #
    # Make the meshgrid
    # 
    # For original grid
    qq = np.meshgrid(rc, thetac, phic, indexing='ij')
    rr = qq[0]  # Spherical R 
    tt = qq[1]  # Angle with z axis ( = pi/2 for the z=0 axis)
    zr = np.pi / 2.e0 - qq[1]
    zz = np.sin(zr) * rr
    nphi = len(phic)
    # For extended grid
    rr_ext = qq_ext[0]  # Spherical R 
    tt_ext = qq_ext[1]  # Angle with z axis ( = pi/2 for the z=0 axis)
    zr_ext = np.pi / 2.e0 - qq_ext[1]
    zz_ext = np.sin(zr_ext) * rr_ext
    #
    # Name all the data in the extended grid
    #
    rr_flat_ext_old = rr_ext[:, :, 0].flatten()
    tt_flat_ext_old = tt_ext[:, :, 0].flatten()
    # zr_flat_ext_old = z_flat_ext_old/x_flat_ext_old
    
    #
    # DENSE grid
    #
    if nr_dense is None:
        nr_dense = nr_cyl
    if ntheta_dense is None:
        ntheta_dense = nz_cyl
        
    # Make Meshgrid
    # rc_dense = np.logspace(np.log10(np.min(rr_ext)), np.log10(np.max(rr_ext)), nr_dense)
    rc_dense = rr_ext[:, -1, 0].copy()
    while len(rc_dense) <= nr_dense:
        rc_dense = np.hstack([rc_dense, 0.5 * (rc_dense[:-1] + rc_dense[1:])])
    rc_dense = np.sort(rc_dense)
    
    thetac_dense = thetac.copy()
    while len(thetac_dense) <= ntheta_dense:
        thetac_dense = np.hstack([thetac_dense, 0.5 * (thetac_dense[:-1] + thetac_dense[1:])])
    thetac_dense = np.sort(thetac_dense)

    qq_dense = np.meshgrid(rc_dense, thetac_dense, phic, indexing='ij')
    rr_dense = qq_dense[0]  # Spherical R 
    tt_dense = qq_dense[1]  # Angle with z axis ( = pi/2 for the z=0 axis)
    # zr_dense = np.pi/2.e0 - qq_dense[1]

    nr_dense = len(rc_dense)
    ntheta_dense = len(thetac_dense)

    print('the DENSE grid: ', nr_dense, ntheta_dense)
    
    #
    # New Grid for Sampling in the scipy griddata regridding
    # The new grid will exactly in r and z/r grids
    # 
    # rc_new = np.logspace(np.log10(0.95 * np.min(rc)), np.log10(1.05 * np.max(rc)), nr_cyl)
    if rc_cyl_set is None: 
        rc_new = np.logspace(np.log10(0.95 * np.min(rc)), np.log10(1.05 * np.max(rc)), nr_cyl)
    else:
        rc_new = rc_cyl_set
        
    zrc_new = np.linspace(np.min(1 / np.tan(tt)), np.max(1 / np.tan(tt)), nz_cyl)
    zrc_new = zrc_new[::-1]
    # zrc_new = 1/np.tan(tt)
    qq_new = np.meshgrid(rc_new, zrc_new, phic, indexing='ij')
    rq_new = qq_new[0]
    zrq_new = qq_new[1]
    
    nr_dust_spec = int(rhodust.shape[3])
    rhodust_cyl_old = np.zeros([nr_cyl, nz_cyl, nphi, nr_dust_spec])
    
    # define the coordinates
    # from coarse sph to dense sph
    grid_dense_old = (np.log10(rr_flat_ext_old / const.au), tt_flat_ext_old)
    grid_dense_new = (np.log10(rr_dense / const.au), tt_dense)
    
    # define the coordinates
    # from dense sph to dense cyl
    r_flat_cyl_old = rr_dense[:, :, 0].flatten() * np.sin(tt_dense[:, :, 0].flatten())
    z_flat_cyl_old = rr_dense[:, :, 0].flatten() * np.cos(tt_dense[:, :, 0].flatten())
    zr_flat_cyl_old = z_flat_cyl_old / r_flat_cyl_old
    grid_cyl_old = (np.log10(r_flat_cyl_old / const.au), zr_flat_cyl_old)
    grid_cyl_new = (np.log10(rq_new / const.au), zrq_new)
    
    for i_dust in range(nr_dust_spec):
        
        rhodust_sph_ext_idust_t = rhodust_sph_ext[:, :, 0, i_dust]
        rhodust_sph_ext_idust_t[rhodust_sph_ext_idust_t <= 0.0] = np.nan
        rhodust_sph_ext_idust_t = fill_nan_in_2d_array_inlog(rhodust_sph_ext_idust_t)
        rhodust_flat_ext_old_idust = rhodust_sph_ext_idust_t.flatten()
        
        #
        # linear interpolate to setup all the cells in spherical coordinate
        #
        # rhodust_sph_dense_old_idust = np.zeros([nr_dense, ntheta_dense, nphi])
        logvalue = np.log10(rhodust_flat_ext_old_idust)
        grid_rhodust_sph = interpolate.griddata(grid_dense_old, logvalue, grid_dense_new, method=method_polate1)
        rhodust_sph_dense_old_idust = 10 ** grid_rhodust_sph.copy()

        #
        # Name the data in extended DENSE grid
        #
        rhodust_sph_dense_old_idust[:, :, 0] = fill_nan_in_2d_array_inlog(rhodust_sph_dense_old_idust[:, :, 0])
        rhodust_flat_sph_dense_old_idust = rhodust_sph_dense_old_idust[:, :, 0].flatten()
        
        #
        # Regridding to a New Grid
        # The Method is 'cubic' because it is an Oversampling
        # But the 'cubic' method will bring nan value and further destroy the results
        # so now adopt 'linear'
        #
        # method_polate = 'cubic' # now it is defined in the input
        # rhodust_cyl_old_idust = np.zeros([nr_cyl, nz_cyl, nphi])
        logvalue = np.log10(rhodust_flat_sph_dense_old_idust)
        grid_rhodust_cyl = interpolate.griddata(grid_cyl_old, logvalue, grid_cyl_new, method=method_polate2)
        rhodust_cyl_old_idust = 10 ** grid_rhodust_cyl.copy()
        
        rhodust_cyl_old[:, :, :, i_dust] = rhodust_cyl_old_idust
    
    return rc_new, zrc_new, qq_new, rhodust_cyl_old


def cyl2sph_dust(rc, thetac, phic, rr_cyl, zrr_cyl, rhodust_cyl_new, nr_dense=None, ntheta_dense=None,
            method_polate1='linear', method_polate2='linear'):
    """
    Performs the conversion from cylindrical coordinates back to spherical coordinates.

    Args:
        rc (ndarray): Original radial grid points as a 1D numpy array.
        thetac (ndarray): Original polar angle grid points as a 1D numpy array.
        phic (ndarray): Original azimuthal angle grid points as a 1D numpy array.
        rr_cyl (ndarray): Radial grid points in the cylindrical coordinate as a 3D numpy array.
        zrr_cyl (ndarray): Vertical grid points (z/r ratio) in the cylindrical coordinate as a 3D numpy array.
        rhodust_cyl_new (ndarray): Dust density distribution in the cylindrical coordinate as a 3D+1D numpy array.
        nr_dense (int, optional): Number of radial grid points in the dense spherical coordinate. Defaults to 500.
        ntheta_dense (int, optional): Number of polar angle grid points in the dense spherical coordinate. Defaults to 500.
        method_polate1 (str, optional): Interpolation method for sampling back to the dense spherical coordinate. Defaults to 'nearest'. (changed to linear on 01/28)
        method_polate2 (str, optional): Interpolation method for regridding to the original spherical coordinate. Defaults to 'linear'.

    Returns:
        tuple: Tuple containing the regridded gas density distribution in the spherical coordinate as a 3D numpy array and the regridded gas temperature distribution in the spherical coordinate as a 3D numpy array.

    Description:
    This function converts the input dust density distribution (`rhodust_cyl_new`) from cylindrical coordinates back to spherical coordinates. It first samples the input distributions to a DENSE spherical coordinate using the specified interpolation method (`method_polate1`). Then, it regrids the DENSE spherical grids to the original (COARSE) spherical grids using the specified interpolation method (`method_polate2`).

    Note:
    - The function requires the numpy module to be imported as `np`.
    - The function assumes that the `interpolate.griddata` function from the scipy.interpolate module is available.
    - The constants (for example `const.au`) need to be defined before calling this function.
    """

    nr, ntheta, nphi = len(rc), len(thetac), len(phic)

    R_flat_cyl = np.sqrt(rr_cyl.flatten() ** 2 + (zrr_cyl.flatten() * rr_cyl.flatten()) ** 2)
    tt_flat_cyl = np.pi / 2 - np.arctan(zrr_cyl.flatten())
    
    #
    # DENSE grid for first sampling back
    #
    if nr_dense is None:
        nr_dense = rr_cyl.shape[0]
    if ntheta_dense is None:
        ntheta_dense = rr_cyl.shape[1]
    
    # Make Meshgrid
    rc_dense = np.logspace(np.log10(np.min(rc)), np.log10(np.max(rc)), nr_dense)
    thetac_dense = thetac.copy()
    while len(thetac_dense) <= ntheta_dense:
        thetac_dense = np.hstack([thetac_dense, 0.5 * (thetac_dense[:-1] + thetac_dense[1:])])
    thetac_dense = np.sort(thetac_dense)

    qq_dense = np.meshgrid(rc_dense, thetac_dense, phic, indexing='ij')
    rr_dense = qq_dense[0]  # Spherical R 
    tt_dense = qq_dense[1]  # Angle with z axis ( = pi/2 for the z=0 axis)
    # zr_dense = np.pi/2.e0 - qq_dense[1]

    nr_dense = len(rc_dense)
    ntheta_dense = len(thetac_dense)
        
    # Make Flat data for further sampling back to COARSE grid
    rr_flat_dense_new = rr_dense[:, :, 0].flatten()
    tt_flat_dense_new = tt_dense[:, :, 0].flatten()
    
    #
    # New Grid for Sampling in the scipy griddata regridding
    # Now the grid is swiched back to the original grids
    # 
    qq = np.meshgrid(rc, thetac, phic, indexing='ij')
    rr = qq[0]  # Spherical R 
    tt = qq[1]  # Angle with z axis ( = pi/2 for the z=0 axis)
    zr = np.pi / 2.e0 - qq[1]
    zz = np.sin(zr) * rr
    
    # define the coordinates
    # from coarse sph to dense sph
    grid_dense_old = (np.log10(R_flat_cyl / const.au), tt_flat_cyl)
    grid_dense_new = (np.log10(rr_dense / const.au), tt_dense)
    
    # define the coordinates
    # from dense sph to dense cyl
    r_flat_cyl_old = rr_dense[:, :, 0].flatten() * np.sin(tt_dense[:, :, 0].flatten())
    z_flat_cyl_old = rr_dense[:, :, 0].flatten() * np.cos(tt_dense[:, :, 0].flatten())
    zr_flat_cyl_old = z_flat_cyl_old / r_flat_cyl_old
    grid_sph_old = (np.log10(rr_flat_dense_new / const.au), tt_flat_dense_new)
    grid_sph_new = (np.log10(rr / const.au), tt)
    
    rhodust_sphr_new = np.zeros([nr, ntheta, nphi, rhodust_cyl_new.shape[-1]])
    for i_dust in range(rhodust_cyl_new.shape[-1]):
        
        rhodust_cyl_idust_t = rhodust_cyl_new[:,:,0,i_dust]
        rhodust_cyl_idust_t = fill_nan_in_2d_array_inlog(rhodust_cyl_idust_t)
        rhodust_i_flat_cyl_new = rhodust_cyl_idust_t.flatten()

        # sampling to the DENSE spherical coordinate
        logvalue = np.log10(rhodust_i_flat_cyl_new)
        grid_dense_rhodust_i = interpolate.griddata(grid_dense_old, logvalue, grid_dense_new, method=method_polate1)
        rhodust_i_sph_dense_new = 10 ** grid_dense_rhodust_i.copy()

        # zr_flat_ext_old = z_flat_ext_old/x_flat_ext_old
        rhodust_i_sph_dense_new[:, :, 0] = fill_nan_in_2d_array_inlog(rhodust_i_sph_dense_new[:, :, 0])
        rhodust_i_dense_new = rhodust_i_sph_dense_new[:, :, 0].flatten()
        
        #
        # Regridding to a New Grid
        # Adopt method of 'linear' because of undersampling
        # 
        # rhodust_i_sphr_new = np.zeros([nr, ntheta, nphi])
        logvalue = np.log10(rhodust_i_dense_new)
        grid_rz_rhodust_i = interpolate.griddata(grid_sph_old, logvalue, grid_sph_new,method=method_polate2)
        rhodust_i_sphr_new = 10 ** grid_rz_rhodust_i.copy()
        # rhodust_i_sphr_new[np.isnan(rhodust_i_sphr_new)] = 0
        rhodust_i_sphr_new[:, :, 0] = fill_nan_in_2d_array_inlog(rhodust_i_sphr_new[:, :, 0])

        rhodust_sphr_new[:,:,:,i_dust] = rhodust_i_sphr_new[:,:,:]
    
    return rhodust_sphr_new


def calc_tau_new(dataT, rr_cyl, zr_cyl, rhodust_cyl, dust_spec_nr, wavelength=0.55,
                 savenames=['taux.dat', 'tauz_up.dat', 'tauz_dn.dat', 'dtauz.dat']):
    """
    **NOTE** currently the scattering coefficients are ignored, because they are neglegible compared to the absorption coeffcients in most cases.
    Calculates the optical depth in cylindrical coordinates. The function calculates the optical depth in three different directions: R direction (taux), upward z direction (tauz_up), and downward z direction (tauz_dn). The calculations are done in cylindrical coordinates.

    Args:
        dataT: The data file read in by radmc3dPy.
        rr_cyl (ndarray): Radial grid points in cylindrical coordinates.
        zr_cyl (ndarray): Vertical grid points (z/r) in cylindrical coordinates.
        rhodust_cyl (ndarray): Dust density distribution in cylindrical coordinates.
        dust_spec_nr (int): Number of dust species.
        wavelength (float, optional): Wavelength of the observation in micrometers (um). Defaults to 0.55.
        savenames (list, optional): List of filenames to save the calculated optical depth values. Defaults to ['taux.dat', 'tauz_up.dat', 'tauz_dn.dat', 'dtauz.dat'].

    Returns:
        tuple: Tuple containing the calculated optical depth values in cylindrical coordinates: taux, tauz_up, tauz_dn, and dtauz.

    Description:
    This function calculates the optical depth in cylindrical coordinates based on the input data. It first calculates the optical depth in the R direction (taux) using the `getTau` function from the `dataT` object by radmc3dPy. Then, it calculates the optical depth in the upward z direction (tauz_up) and downward z direction (tauz_dn) in cylindrical coordinates. Finally, it saves the calculated optical depth values to the specified filenames.
    
    In short,
    The New way to calculate Tau is done by following:
    
    1. taux (used radmc3dPy methods) -- which is the tau in the R direction to the star, in spherical coordinate
    2. tauz up and down -- tau in z direction is now calculated in cylindrical coordinate

    Note:
    - The function assumes that the `radmc3dPy` module is available and that the `radmc3dDustOpac` class from this module is imported.
    - The function requires the numpy module to be imported as `np`.
    - The constant `const.au` needs to be defined before calling this function.
    """

    ## dataT.getTau(wav=0.55)

    dataT.getTau(idust=[i for i in range(dust_spec_nr)], wav=wavelength)
    taux = dataT.taux[:, :, 0]

    ##################################################
    # Calculate Everything in Cylindrical Coordinate #
    ##################################################

    wav = wavelength # 0.55 #u.um
    tauz_up_cyl = np.zeros([len(rr_cyl), len(zr_cyl)], dtype=np.float64)
    tauz_dn_cyl = np.zeros([len(rr_cyl), len(zr_cyl)], dtype=np.float64)
    dtauz_cyl = np.zeros([len(rr_cyl), len(zr_cyl)], dtype=np.float64)

    #
    # Read Kappa Value at V-band for each dust specie
    #
    for i_dust in range(dust_spec_nr):
        # for i_dust in [6]:
        # opac = radmc3dPy.analyze.radmc3dDustOpac()
        opac = dustopac.radmc3dDustOpac()
        mo = opac.readMasterOpac()
        ext = mo['ext']
        scatmat = mo['scatmat']
        opac.readOpac(ext=ext[i_dust], scatmat=None)

        #     taur_idust = np.zeros([dataT.grid.nx, dataT.grid.ny, dataT.grid.nz], dtype = np.float64)
        #     tauz_idust = np.zeros([dataT.grid.nx, dataT.grid.ny, dataT.grid.nz], dtype = np.float64)

        if opac.wav[0][-1] > opac.wav[0][0]:
            kabs = 10. ** np.interp(np.log10(np.array(wav)), np.log10(opac.wav[0]), np.log10(opac.kabs[0]))
        else:
            kabs = 10. ** np.interp(np.log10(np.array(wav)), np.log10(opac.wav[0][::-1]),
                                    np.log10(opac.kabs[0][::-1]))

        if opac.ksca[0][0] > 0:
            if opac.wav[0][-1] > opac.wav[0][0]:
                ksca = 10. ** np.interp(np.log10(np.array(wav)), np.log10(opac.wav[0]), np.log10(opac.ksca[0]))
            else:
                ksca = 10. ** np.interp(np.log10(np.array(wav)), np.log10(opac.wav[0][::-1]),
                                        np.log10(opac.ksca[0][::-1]))

        print('Opacity at ' + ("%.2f" % wav) + 'um : ', kabs + ksca)

        kappa = kabs + ksca

        # here z is really the z in the cylindrical coordinate

        dum_x = rr_cyl.copy()
        tauz_idust_up = np.zeros([len(rr_cyl), len(zr_cyl)], dtype=np.float64)
        tauz_idust_dn = np.zeros([len(rr_cyl), len(zr_cyl)], dtype=np.float64)
        dtauz_idust = np.zeros([len(rr_cyl), len(zr_cyl)], dtype=np.float64)
        diff_zr = np.abs(zr_cyl[1:] - zr_cyl[:-1])

        # calculate the dtau, not considering the dz and the final tauz now.
        # for better comparison, provide both dtau and dz in the end
        for iz in range(0, len(zr_cyl)):
            dtauz_idust[:, iz] = rhodust_cyl[:, iz, 0, i_dust] * kappa * diff_zr[iz - 1] * dum_x

        tauz_idust_up[:, 0] = rhodust_cyl[:, 0, 0, i_dust] * kappa * diff_zr[0] * dum_x
        for iz in range(1, len(zr_cyl)):
            tauz_idust_up[:, iz] = tauz_idust_up[:, iz - 1] + rhodust_cyl[:, iz, 0, i_dust] * kappa * diff_zr[
                iz - 1] * dum_x

        tauz_idust_dn[:, -1] = rhodust_cyl[:, -1, 0, i_dust] * kappa * diff_zr[-1] * dum_x
        for iz in range(1, len(zr_cyl)):
            tauz_idust_dn[:, -iz - 1] = tauz_idust_dn[:, -iz] + rhodust_cyl[:, -iz - 1, 0, i_dust] * kappa * diff_zr[
                -iz] * dum_x

        # summing up all of the dust subspecies
        dtauz_cyl = dtauz_cyl + dtauz_idust
        tauz_up_cyl = tauz_up_cyl + tauz_idust_up
        tauz_dn_cyl = tauz_dn_cyl + tauz_idust_dn

    # as there is another half of the disk, 
    # so the tauz_dn needs to sum the tauz_up at the bottom
    tauz_up_cyl[np.isnan(tauz_up_cyl)] = 1e-300
    for i_r in range(len(rr_cyl)):
        tauz_dn_cyl[i_r, :] = tauz_dn_cyl[i_r, :] + tauz_up_cyl[i_r, -1]

    np.savetxt(savenames[0], taux)
    np.savetxt(savenames[1], tauz_up_cyl)
    np.savetxt(savenames[2], tauz_dn_cyl)
    np.savetxt(savenames[3], dtauz_cyl)

    return taux, tauz_up_cyl, tauz_dn_cyl, dtauz_cyl


def set_dusttemp_at_midplane(mint, dustdens, Tdust, wav=1330, tau_threshold=1000):
    """
    Sets the dust temperature at the midplane based on the optical depth criteria.

    Args:
        mint (object): An object containing simulation parameters, including the number of dust species.
        dustdens (object): the radmc3dPy object that handles dust density and optical depth calculations. Must have a `getTau` method and a `taux` attribute.
        Tdust (ndarray): 4D array of dust temperatures with shape (nr, ntheta, nz, dust_spec_nr).
        wav (float, optional): Wavelength of the observation in micrometers (um). Defaults to 1330 (for mm C18O wavelength); updated to optical wavelength (0.55) now.
        tau_threshold (float, optional): The optical depth threshold to determine the dust temperature edge. Defaults to 10 (updated to 100 (then 1000) with optical wave).

    Returns:
        ndarray: Updated 4D array of dust temperatures with the same shape as `Tdust`.

    Description:
        This function updates the dust temperature (`Tdust`) at the midplane based on the optical depth (`taux_mm`). For each dust species, it identifies regions where the optical depth is below a specified threshold (`tau_threshold`). In these regions, the dust temperature is retained; otherwise, it is set to the last valid temperature value where the optical depth was below the threshold. This effectively propagates the temperature from the edge of the optical depth threshold inward.

    Raises:
        AttributeError: If `dustdens` does not have the required `getTau` method or `taux` attribute.
        ValueError: If the shapes of input arrays do not match expected dimensions.
    """
    
    # Validate input dimensions
    if Tdust.ndim != 4:
        raise ValueError(f"Expected Tdust to be a 4D array, got {Tdust.ndim}D array instead.")
    
    nr, ntheta, nz, dust_spec_nr = Tdust.shape
    
    # Make a copy of Tdust to preserve the original data
    Tdust_new = Tdust.copy()

    # Retrieve the number of dust species from `mint`
    dust_spec_nr_mint = getattr(mint, 'dust_spec_nr', None)
    if dust_spec_nr_mint is None:
        raise AttributeError("The 'mint' object must have a 'dust_spec_nr' attribute.")
    if dust_spec_nr_mint != dust_spec_nr:
        raise ValueError(f"Mismatch in number of dust species: mint.dust_spec_nr={dust_spec_nr_mint} vs Tdust.shape[3]={dust_spec_nr}.")

    # Calculate optical depth using the provided `dustdens` object
    if not hasattr(dustdens, 'getTau') or not hasattr(dustdens, 'taux'):
        raise AttributeError("The 'dustdens' object must have 'getTau' method and 'taux' attribute.")
    
    # Assuming `getTau` updates the `taux` attribute
    dustdens.getTau(idust=list(range(dust_spec_nr)), wav=wav)
    
    taux_mm = dustdens.taux[:, :, 0]  # Extracting midplane (z=0) optical depth
    if taux_mm.shape != (nr, ntheta):
        raise ValueError(f"Expected taux_mm shape to be ({nr}, {ntheta}), got {taux_mm.shape} instead.")

    # Create a mask where taux_mm is less than the threshold
    mask = taux_mm < tau_threshold  # Shape: (nr, ntheta)

    for i_dust in range(dust_spec_nr):
        # Extract the Tdust values for the current dust species at midplane
        Tdust_species = Tdust[:, :, 0, i_dust]  # Shape: (nr, ntheta)
        
        if Tdust_species.shape != (nr, ntheta):
            raise ValueError(f"Expected Tdust_species shape to be ({nr}, {ntheta}), got {Tdust_species.shape} instead.")

        # Use where mask is True, else NaN
        Tdust_filled = np.where(mask, Tdust_species, np.nan)
        
        # Forward fill NaN values along the theta axis
        mask_valid = ~np.isnan(Tdust_filled)
        
        # Create an array of indices along theta
        theta_indices = np.arange(ntheta)
        
        # Replace invalid entries with the last valid index
        idx = np.where(mask_valid, theta_indices, 0)
        idx = np.maximum.accumulate(idx, axis=1)
        
        # Create row indices for advanced indexing
        row_indices = np.arange(nr)[:, np.newaxis]
        
        # Assign the forward-filled values to Tdust_new
        Tdust_new[:, :, 0, i_dust] = Tdust_filled[row_indices, idx]
        
    return Tdust_new


def get_new_dust_temperature(mint, bool_clean_dusttemp_midplane=True, command='radmc3d mctherm'):
    """
    Calculates and updates the dust temperature by running RADMC-3D's Monte Carlo thermal simulation and optionally cleans up noise in the midplane.

    Args:
        mint (object): An object containing simulation parameters, including dust species information and fractions.
        bool_clean_dusttemp_midplane (bool, optional): Flag to determine whether to clean noise in the midplane of the dust temperature. Defaults to True.
        command (str, optional): Shell command to execute RADMC-3D's Monte Carlo thermal simulation. Defaults to 'radmc3d mctherm'.

    Returns:
        int: Returns 1 upon successful completion.

    Description:
        This function performs the following steps:
        
        1. **Run RADMC-3D Monte Carlo Thermal Simulation**:
            - Executes the provided shell command (`command`) to run RADMC-3D's `mctherm` mode, which calculates the dust temperature distribution.
        
        2. **Optional Midplane Noise Cleanup**:
            - If `bool_clean_dusttemp_midplane` is set to `True`, the function proceeds to clean up noise in the midplane of the dust temperature distribution.
            - **Steps for Cleaning Midplane Noise**:
                a. **Read Dust Density and Temperature Data**:
                    - Utilizes `radmc3dPy.analyze.readData` to read the dust density (`ddens`) and dust temperature (`dtemp`) data from RADMC-3D output files.
                
                b. **Update Dust Temperature at Midplane**:
                    - Copies the original dust temperature data and passes it to the `set_dusttemp_at_midplane` function along with the `mint` and `dustdens` objects to obtain a cleaned version of the dust temperature (`Tdust_new`).
                
                c. **Backup Original Dust Temperature File**:
                    - Renames the original `dust_temperature.dat` to `dust_temperature_originalmidplane.dat` to preserve the uncleaned data.
                
                d. **Write Cleaned Dust Temperature Data**:
                    - Opens a new `dust_temperature.dat` file and writes the cleaned dust temperature data in the required format for RADMC-3D.
                    - The data is written in Fortran-style (column-major order) to ensure compatibility with RADMC-3D's expectations.
        
        **Note**: Ensure that the `set_dusttemp_at_midplane` function is defined and accessible within the scope where this function is used.

    Raises:
        FileNotFoundError: If the RADMC-3D output files (`dust_density.inp`, `dust_temperature.dat`) are not found.
        RuntimeError: If the RADMC-3D command execution fails.
        AttributeError: If the `mint` object lacks required attributes such as `dust_spec_nr` or `fracs`.
        ValueError: If the shapes of the temperature arrays do not match expected dimensions.

    Example:
        ```python
        # Assuming `mint` and `dustdens` objects are properly initialized
        result = get_new_dust_temperature(mint)
        if result == 1:
            print("Dust temperature updated successfully.")
        else:
            print("An error occurred while updating dust temperature.")
        ```
    """
    try:
        # Run RADMC-3D Monte Carlo Thermal Simulation
        if isinstance(command, str):
            command_status = os.system(command)
            if command_status != 0:
                raise RuntimeError(f"Command '{command}' failed with status {command_status}.")
        else:
            raise TypeError("The 'command' parameter must be a string.")

        # Optional: Clean noise in the midplane of the dust temperature
        if bool_clean_dusttemp_midplane:
            # Read dust temperature and density data from RADMC-3D output
            try:
                dustdens = modelgrid.readData(dtemp=True, ddens=True) # , binary=False)
            except FileNotFoundError as e:
                raise FileNotFoundError("Required RADMC-3D output files not found.") from e

            # Copy original dust temperature data
            Tdust = dustdens.dusttemp.copy()

            # Clean the dust temperature at the midplane
            Tdust_new = set_dusttemp_at_midplane(mint, dustdens, Tdust)

            # Backup the original dust_temperature.dat file
            backup_command = 'mv dust_temperature.dat dust_temperature_originalmidplane.dat'
            backup_status = os.system(backup_command)
            if backup_status != 0:
                raise RuntimeError(f"Failed to backup dust_temperature.dat with status {backup_status}.")

            # Retrieve dimensions from the dust temperature array
            nr, ntheta, nphi, nrspec = Tdust.shape

            # Write the cleaned dust temperature data to a new dust_temperature.dat file
            try:
                with open('dust_temperature.dat', 'w') as f:
                    f.write('1\n')  # Format number
                    f.write(f"{nr * ntheta * nphi}\n")  # Number of cells
                    f.write(f"{nrspec}\n")  # Number of dust species

                    for i in range(nrspec):
                        # Flatten the 3D dust temperature array for the current species in Fortran order
                        data = Tdust_new[:, :, :, i].ravel(order='F')
                        # Write the flattened data to the file with scientific notation
                        np.savetxt(f, data, fmt="%13.6e")
                        f.write('\n')  # Ensure separation between species
            except IOError as e:
                raise IOError("Failed to write the cleaned dust_temperature.dat file.") from e

        return 1

    except Exception as e:
        # Handle unexpected exceptions and provide informative error messages
        print(f"Error in get_new_dust_temperature: {e}")
        return 0


def smooth_gas_temp(rc, thetac, Tg_init, nr_regridr=0.2, verbose=False):
    """
    Smooths the gas temperature distribution to obtain a noise-free result.

    **IMPORTANT**
    Becarefull using the smooth_gas_temp function!
    This will only be usefull when we are using in-sufficient nphot number so that there are some noise in mid-plane.
    We always need to compare the temperature smoothed inner disk with the case setting up large enough nphot.
    In the final result, if sufficient nphot is setup, then we don't need to do smooth.

    Args:
        rc (ndarray): Radial grid points.
        thetac (ndarray): Polar grid points.
        Tg_init (ndarray): Initial gas temperature distribution (2D array).
        nr_regridr (float, optional): Number specifying the regridding factor for the radial grid. Defaults to 0.2.
        verbose (bool): whether to print on screen for debugging

    Returns:
        ndarray: Smoothed gas temperature distribution.

    Description:
    This function is used to smooth the gas temperature distribution to remove noise. It takes as input the initial gas temperature distribution (Tg_init) and the radial (rc) and polar (thetac) grid points. The gas temperature distribution is first regridded to a new grid with a higher resolution in the radial direction (specified by nr_regridr). The regridded temperature distribution is then smoothed using linear interpolation. Finally, the smoothed temperature distribution is regridded back to the original grid.

    Note:
    - The function requires the numpy module to be imported as `np`.
    - The function assumes that the `interpolate.griddata` function from the scipy module is available.
    """

    #
    # Old Grid
    #
    qq = np.meshgrid(rc, thetac, indexing='ij')
    r_q = qq[0]
    t_q = qq[1]
    zr_q = np.pi / 2.e0 - qq[1]

    #
    # New Grid for Sampling in the scipy griddata regridding
    #
    rc_new = np.logspace(np.log10(np.min(rc)), np.log10(np.max(rc)), int(len(rc) * nr_regridr))
    zrc = np.pi / 2.e0 - thetac
    #     zrc_new = np.sort(np.hstack([zrc, 0.5*(zrc[1:] + zrc[:-1])]))
    zrc_new = np.sort(zrc)
    zrc_new = zrc_new[::-1]
    qq_new = np.meshgrid(rc_new, zrc_new, indexing='ij')
    rq_new, zr_q_new = qq_new

    #
    # Regridding to a New Grid
    #
    Tg_f_all = Tg_init.copy()
    grid_zt = interpolate.griddata((r_q.flatten(), zr_q.flatten()), Tg_f_all.flatten(), (rq_new, zr_q_new),
                                   method='linear')

    # remove the nan values at the edge of the grid and put them to the nearest values
    grid_zt_nearest = interpolate.griddata((r_q.flatten(), zr_q.flatten()), Tg_f_all.flatten(), (rq_new, zr_q_new),
                                   method='nearest')
    nanplaces = np.isnan(grid_zt)
    grid_zt[nanplaces] = grid_zt_nearest[nanplaces]
    
    if verbose:
        print('for Tgas smoothing.')
        print('checking the lowest layer at z -> 0')
        print(rq_new[:, -1])
        print(r_q[:, -1])
        print('original T')
        print(Tg_f_all[:, -1])
        print('resampled T')
        print(grid_zt[:, -1])
    
    #
    # Regridding Back to the OLD Grid
    #  
    # Tg_f_all_nRG = Tg_f_all.copy()
    Tg_f_all = interpolate.griddata((rq_new.flatten(), zr_q_new.flatten()), grid_zt.flatten(), (r_q, zr_q),
                                    method='linear')

    if verbose:
        print('smoothed T')
        print(Tg_f_all[:, -1])
        
    return Tg_f_all


def solve_gastemp(mint, dustdens, Tdust, 
                  bool_smoothing_Tgas=True, 
                  bool_visc_heating=True, 
                  bool_output_gastemp=False,
                  bool_temp_decouple=False, 
                  verbose=False):
    """
    calculate the gas temperature based on the dust-gas thermal equilibrium together with the viscous heating
    
    Args:
        mint (object): The DiskMINT model object containing the parameters and data.
        Tdust: (3D array) dust thermal structure
        dustdens: the rdmc3dPy dust density object that read by `radmc3dPy.analyze.readData()`
        bool_smoothing_Tgas (bool): whether to smooth the gas temperature so that we can get a more noise-free Tgas 
                                    structure especially in the inner disk close to the mid-plane. Default True.
        bool_visc_heating (bool): whether to add viscous heating
        bool_output_gastemp (bool): whether to directly output gas temperature, the saved file would be named as `gas_temperature_direct_output.inp`
        bool_temp_decouple (bool): whether to increase the Tgas with certain amount from the ave(Tdust)
                                   the way we make it is to have a Transiaction radius R_temp_trans (in the unit of au), 
                                   and the factor to which the gas temperature is larger than the dust temperature fact_Tgas_2_Tdust
        verbose (bool): if true, output the prints to screen for debugging
        
    Return:
        Tgas: (3D array) the solved gas temperature structure
    """
    
    rhodust = dustdens.rhodust.copy()
    
    nrspec = int(dustdens.rhodust.shape[-1])
    # Use this new defination
    # rhogas = mint.ratio_g2d * np.sum(dustdens.rhodust, axis=3)
    # rhogas = dustdens.rhodust[:,:,:,0]/fracs[0] * ratio_g2d
    # ngas = rhogas / const.mu
    nr = dustdens.grid.nx
    ntheta = dustdens.grid.ny
    nphi = dustdens.grid.nz
    rc = dustdens.grid.x
    thetac = dustdens.grid.y
    phic = dustdens.grid.z
    
    para = mint.parameters
    
    if verbose:
        # print the current working dir
        print(os.getcwd())
    
    #
    # Read ndsd [call to fortran]
    #
    #     ndsd = np.loadtxt('ndsd.inp') # n_d(a) pi a^2
    #     if(len(ndsd)!=nrspec):
    #         print('Problem..')
    #         quit()

    # mstar = mstar
    # Tdust[Tdust < 2.7] = 2.7
    Tgas = np.zeros([nr, ntheta, nphi])

    #
    # Compute the gas temperature
    #
    if verbose:
        print('aave:')
        print(mint.a_ave)
        # print('ndsd:')
        # print(mint.ndsd)
        print('size of Tdust: ', Tdust.shape)
        print('size of Tgas (prepared)', Tgas.shape) 
    
    #
    # Compute the ndsd
    #
    # ndsd = nd * ng * a_ave**2
    # ndsd = fracs_nb * a_ave**2 # number fractions for fracs
    # ndsd_t = mint.ndsd 
    # Compute ndust in a vectorized way
    ndust = rhodust / (para.rhobulk * (4/3) * np.pi * (np.array(mint.a_ave) / 2) ** 3)

    # Compute fracs_nb_eachsize only for the third dimension index 0
    fracs_nb_eachsize = np.zeros_like(ndust)
    ndsd_eachsize = np.zeros_like(ndust)

    # Compute fraction and ndsd only for the third dimension index 0
    sum_ndust = np.nansum(ndust[:, :, 0, :], axis=-1, keepdims=True)  # Sum over last axis
    fracs_nb_eachsize[:, :, 0, :] = ndust[:, :, 0, :] / sum_ndust  # Broadcasting applies
    ndsd_eachsize[:, :, 0, :] = fracs_nb_eachsize[:, :, 0, :] * np.array(mint.a_ave)**2  # Broadcasting applies
    
    #
    # solve the Tgas at different location
    #
    for i in range(0, nr):
        for j in range(0, ntheta):

            Tdmax_here = np.nanmax(Tdust[i, j, 0, :])
            ndsd_t = ndsd_eachsize[i, j, 0, :]
            
            if verbose:
                try: 
                    Tgas[i, j, 0] = scipy.optimize.bisect(func_Tgas, 2.7, Tdmax_here, args=(nrspec, Tdust[i, j, 0, :], ndsd_t))
                    
                except:
                    print(Tdmax_here)
                    print(Tdust[i,j,0,:])
                    print(ndsd_t)
                    
            else:
                Tgas[i, j, 0] = scipy.optimize.bisect(func_Tgas, 2.7, Tdmax_here, args=(nrspec, Tdust[i, j, 0, :], ndsd_t))
            
    if verbose:
        print('solved the Tgas')

    #
    # Smooth the Tg
    #
    if bool_smoothing_Tgas:
        Tgas_new = np.zeros([nr, ntheta, 1])
        Tgas_new[:, :, 0] = smooth_gas_temp(rc / const.au, thetac, Tgas[:, :, 0], nr_regridr=0.2)

        Tgas = Tgas_new.copy()
        
        if verbose:
            print('smoothed Tgas')

    """
    Modifications to the solved Tgas
    """
    
    qq = np.meshgrid(rc, thetac, phic, indexing='ij')
    rr = qq[0]  # Spherical R 
    tt = qq[1]  # Angle with z axis ( = pi/2 for the z=0 axis)
    zr = np.pi / 2.e0 - qq[1]
    
    r_cyl_t = rr * np.cos(zr)

    #
    # (Devleoping)
    # set up a transaction temperature
    # we find that the Tgas and Tdust would decoupled within certain radius
    # therefore we need to increase the Tgas with certain amount from the ave(Tdust)
    # the way we make it is to have a Transiaction radius R_temp_trans (in the unit of au), 
    # and the factor to which the gas temperature is larger than the dust temperature fact_Tgas_2_Tdust
    #
    if bool_temp_decouple:
        R_temp_trans = para.R_temp_trans
        fact_Tgas_2_Tdust = para.fact_Tgas_2_Tdust
        
        # # Tgas_new[:, :, 0] = smooth_gas_temp(rc / const.au, thetac, Tgas[:, :, 0], nr_regridr=0.2)
        # index_radius = rc < R_temp_trans
        # Tgas[index_radius, :, 0] = Tgas[index_radius, :, 0] * fact_Tgas_2_Tdust
        
        """
        **NOTE**
        now the vertical cut of 0.12 is a hard cut for testing
        in the future I will update the number according to the Av values
        
        Below is
        flavor ZERO, simply set the jump for transition
        """
        # index_zr     = (np.pi/2.0 - thetac) > 0.12
        # index_radius = rc/const.au < R_temp_trans

        # Tgas[np.ix_(index_radius, index_zr, [0])] *= fact_Tgas_2_Tdust
        
        """
        **NOTE**
        flavor ONE, set up the transition smoother
        """
        # Calculate the condition based on thetac
        index_zr = (np.pi/2.0 - thetac) > 0.02
        
        # Assuming rc, const, Tgas, thetac, and R_temp_trans are defined
        # Define the transition range around R_temp_trans
        transition_start = R_temp_trans - 0.9*R_temp_trans
        transition_end = R_temp_trans

        # Convert rc to the same units as R_temp_trans for comparison (assuming au)
        r_in_au = rc / const.au

        # Initialize the factor array with ones
        smooth_factor = np.ones_like(r_in_au)

        # Indices for the smooth transition based on radius
        transition_indices = (r_in_au >= transition_start) & (r_in_au <= transition_end)

        # Calculate the smooth factor for the transition range (linear for simplicity)
        # Linear transition: (r - start) / (end - start) * (fact - 1) + 1
        smooth_factor[transition_indices] = ((transition_end - r_in_au[transition_indices]) /
                                            (transition_end - transition_start) *
                                            (fact_Tgas_2_Tdust - 1)) + 1

        # Apply the full factor beyond the transition end
        smooth_factor[r_in_au < transition_start] = fact_Tgas_2_Tdust

        # Now, apply this factor to Tgas, considering both the radius and thetac conditions
        for i in range(len(r_in_au)):
            for j in range(len(index_zr)):
                if index_zr[j]:
                    # Multiply by the smooth factor if index_zr[j] is True
                    Tgas[i, j, 0] *= smooth_factor[i]

    """
    Add viscous heating to the Tgas computation
    
    Find Tgas the way you do now with the balance between heating and cooling by grains: let us call this Tgas_d for now.

    Then, estimate a viscous term as follows (this assumes all the heat is radiated away as a black body): 
    $Tgas_{visc} = ( fac1 * fac2 ) ** (1./4.)$ where $fac1 = 3 G Mstar Mdotacc /(8 \pi \sigma_{SB} r**3)$ and $fac2 = 1 - sqrt(rstar/r)$. Here $G=6.67e-8$, use the mstar, rstar and mdot for RU Lup. sigma_SB is the $stefan boltzmann constant = 5.67e-5$ and $r$ is the cylindrical radius. Everything is in cgs units.

    Then set $Tgas**4 = ( Tgas_{visc}**4 + Tgas_d **4)$
    """
    
    #
    # Compute the Viscous Heating Tgas
    #
    if bool_visc_heating:
        Tgas_d = Tgas.copy()

        fac1 = 3 * const.gg * para.mstar * para.mdotacc / (8 * np.pi * const.sigma_SB * r_cyl_t ** 3)
        fac2 = 1 - np.sqrt(para.rstar / r_cyl_t)
        Tgas_visc = (fac1 * fac2) ** (1 / 4)

        Tgas = (Tgas_d ** 4 + Tgas_visc ** 4) ** (1 / 4)

    print('DONE! setup Tgas')
    
    if bool_output_gastemp:
        # write the gastemp file as the inp file
        with open('gas_temperature_direct_output.new', 'w+') as f:
            f.write('1\n')  # Format number
            f.write('%d\n' % (nr * ntheta * nphi))  # Nr of cells
            data = Tgas.ravel(order='F')  # Create a 1-D view, fortran-style indexing
            data.tofile(f, sep='\n', format="%13.6e")
            f.write('\n')
    
    return Tgas


def get_new_dust_density(mint, bool_smoothing_Tgas=True,
                         bool_dust_vhse=True,
                         bool_dust_fragmentation_set=False,
                         bool_dust_radial_drifting_set=False, 
                         bool_dust_settling_in_iteration=False,
                         verbose=True):
    """
    Calculates the new density distribution by solving the VHSE (Vertical Hydrostatic Equilibrium).
    This is just a single interation of solving VHSE, the stable/final VHSE solution requires multiple iterations.
    This function requires the radmc3dPy module.

    Args:
        model (object): The DiskMINT model object containing the parameters and data.
        bool_smoothing_Tgas (bool): whether to smooth the gas temperature so that we can get a more noise-free Tgas structure especially in the inner disk close to the mid-plane. Default True.
        bool_dust_vhse (bool): whether to solve the vertical hydrastatic equilibrium (VHSE). Default True.
        bool_dust_fragmentation_set (bool): whether to turn on the dust fragmentation feature. Default False.
        bool_dust_radial_drifting_set (bool): whether to turn on the dust radial drift feature. Default False.
        
        **NOTE: Developing**
        bool_dust_settling_in_iteration (bool): whether to solve the dust settling in each iteration when solving the VHSE.
        verbose (bool): whether to output the prints.
        
    Returns:
        None.
        
        The new density distribution in RADMC3D format will be saved to the file_dir.

    Description:
    This function calculates the new density distribution by solving the VHSE equation, following the steps:
    
    0. setting up: It reads the density and temperature data from the dust_density.inp and dust_temperature.inp files, respectively. It then computes the gas temperature (Tgas). (optional) It smooths Tgas using the `smooth_gas_temp` function. It also computes the viscous heating term and adds it to the gas temperature.
    1. converting from spherical coordinate to cylindrical: The density distribution is then transformed to cylindrical coordinates using the `sph2cyl` function.
    2. solving VHSE: The new density distribution is solved in VHSE, and the density structure is updated based on the gas temperature. Then the new density structure is re-normalized in each vertical column to ensure the surface density distribution is not changed.
    3. converting from cylindrical coordinate to spherical: The new density distribution is converted back to spherical coordinates.
    4. saving files: The new density distribution saved to the dust_density.new and gas_temperature.new files.

    Note:
    - The function assumes that the `interpolate.griddata` function from the scipy module is available.
    - The function assumes that the `scipy.optimize.bisect` function is available.
    - The function assumes that the constants `const.gg`, `const.mu`, `const.kk`, `const.sigma_SB`, and `const.au` are defined and imported.
    - The function assumes that the `smooth_gas_temp`, `sph2cyl`, `cyl2sph`, and `func_Tgas` functions are defined and imported.
    - The function assumes that the `np` module is imported.
    - The function assumes that the `os` module is imported.
    - The function assumes that the `scipy.optimize` module is imported.
    - The function assumes that the `radmc3dPy.analyze.readData` function is available.
    """

    para = mint.parameters

    # print the current working dir
    print(os.getcwd())

    # 
    # Read dust_density.inp
    #
    # dustdens = read_dustdens()
    # dustdens = radmc3dPy.analyze.readData(dtemp=True, ddens=True, binary=False)
    dustdens = modelgrid.readData(dtemp=True, ddens=True) # , binary=False)
    # nrspec = int(dustdens.nrspec)
    nrspec = int(dustdens.rhodust.shape[-1])
    # Use this new defination
    rhogas = mint.ratio_g2d * np.nansum(dustdens.rhodust, axis=3)
    rhogas[np.isnan(rhogas)] = 1e-199
    rhogas[rhogas <= const.mu] = const.mu
    # rhogas = dustdens.rhodust[:,:,:,0]/fracs[0] * ratio_g2d
    ngas = rhogas / const.mu
    nr = dustdens.grid.nx
    ntheta = dustdens.grid.ny
    nphi = dustdens.grid.nz
    rc = dustdens.grid.x
    thetac = dustdens.grid.y
    phic = dustdens.grid.z
    qq = np.meshgrid(rc, thetac, phic, indexing='ij')
    rr = qq[0]  # Spherical R 
    tt = qq[1]  # Angle with z axis ( = pi/2 for the z=0 axis)
    zr = np.pi / 2.e0 - qq[1]
    zz = np.sin(zr) * rr
    # 
    # Read dust_density.inp
    #
    # dusttemp = read_dusttemp()
    Tdust = dustdens.dusttemp.copy()
    #
    # Read ndsd [call to fortran]
    #
    #     ndsd = np.loadtxt('ndsd.inp') # n_d(a) pi a^2
    #     if(len(ndsd)!=nrspec):
    #         print('Problem..')
    #         quit()

    # mstar = mstar
    Tdust[Tdust < 2.7] = 2.7
    print('DONE! READ in dust thermal structure from previous iteration')
    
    mass_dust = (np.sum(dustdens.rhodust, axis=3) * mint.vol).sum(0).sum(0).sum(0)
    print('dust mass read in = %.3e'%(mass_dust * 2 / const.ms))
    
    mass_gas = (rhogas * mint.vol).sum(0).sum(0).sum(0)
    print('gas mass read in = %.3e'%(mass_gas * 2 / const.ms))
    
    #
    # Compute the gas temperature
    #
    Tgas = solve_gastemp(mint, dustdens, Tdust, 
                         bool_smoothing_Tgas=bool_smoothing_Tgas,
                         bool_temp_decouple=mint.bool_temp_decouple)
    
    # Tgas = np.zeros([nr, ntheta, nphi])
    # #  print(para.ndsd)
    # print(mint.ndsd)
    # for i in range(0, nr):
    #     for j in range(0, ntheta):
    #         #             ndsd = nd * ng * a_ave**2
    #         # ndsd = fracs_nb * a_ave**2 # number fractions for fracs
    #         ndsd_t = mint.ndsd

    #         Tdmax_here = np.max(Tdust[i, j, 0, :])
    #         #  try:
    #         Tgas[i, j, 0] = scipy.optimize.bisect(func_Tgas, 2.7, Tdmax_here, args=(nrspec, Tdust[i, j, 0, :], ndsd_t))
    #         #  except:
    #         #  print(Tdmax_here)
    #         #  print(Tdust[i,j,0,:])
    #         #  print(ndsd_t)

    # #
    # # Smooth the Tg
    # #
    # if bool_smoothing_Tgas:
    #     Tgas_new = np.zeros([nr, ntheta, 1])
    #     Tgas_new[:, :, 0] = smooth_gas_temp(rc / const.au, thetac, Tgas[:, :, 0], nr_regridr=0.2)

    #     Tgas = Tgas_new.copy()

    # """
    # Add viscous heating to the Tgas computation
    
    # Find Tgas the way you do now with the balance between heating and cooling by grains: let us call this Tgas_d for now.

    # Then, estimate a viscous term as follows (this assumes all the heat is radiated away as a black body): 
    # $Tgas_{visc} = ( fac1 * fac2 ) ** (1./4.)$ where $fac1 = 3 G Mstar Mdotacc /(8 \pi \sigma_{SB} r**3)$ and $fac2 = 1 - sqrt(rstar/r)$. Here $G=6.67e-8$, use the mstar, rstar and mdot for RU Lup. sigma_SB is the $stefan boltzmann constant = 5.67e-5$ and $r$ is the cylindrical radius. Everything is in cgs units.

    # Then set $Tgas**4 = ( Tgas_{visc}**4 + Tgas_d **4)$
    # """

    # r_cyl_t = rr * np.cos(zr)

    # #
    # # Compute the Viscous Heating Tgas
    # #

    # Tgas_d = Tgas.copy()

    # fac1 = 3 * const.gg * para.mstar * para.mdotacc / (8 * np.pi * const.sigma_SB * r_cyl_t ** 3)
    # fac2 = 1 - np.sqrt(para.rstar / r_cyl_t)
    # Tgas_visc = (fac1 * fac2) ** (1 / 4)

    # Tgas = (Tgas_d ** 4 + Tgas_visc ** 4) ** (1 / 4)

    print('DONE! READ in all the previous iterations and setup Tgas')

    #
    # set the edge values to 1 moleculr / cm^3
    # because they are just the inner most grids that don't contribute in the end
    #
    rhogas_orig = rhogas.copy()
    rhogas[rhogas <= const.mu] = const.mu
    
    """Solve the VHSE if needed"""
    if bool_dust_vhse:
        #######################################
        # 1.1). Extended Spherical Coordinate #
        # 1.2). Convert to Cylindrical Cord   #
        #######################################

        rc_cyl, zrc_cyl, qq_cyl, rhogas_cyl_old, Tgas_cyl_old = \
            sph2cyl(rc, thetac, phic, rhogas, Tgas, nr_cyl=mint.nr_cyl_insitu, nz_cyl=mint.ntheta_cyl_insitu)

        rhogas0 = np.zeros(rhogas.shape)
        rhogas0[:, :, 0] = np.loadtxt('rhogas_probset.inp')

        rhogas_sph_setup = rhogas0.copy()

        rhogas_sph_setup[rhogas_sph_setup <= const.mu] = const.mu
        print('Grid in Sphrical Coordinate: ', rhogas_sph_setup.shape)

        rc_cyl, zrc_cyl, qq_cyl, rhogas_cyl_setup, Tgas_cyl_old = \
            sph2cyl(rc, thetac, phic, rhogas_sph_setup, Tgas, nr_cyl=mint.nr_cyl_insitu, nz_cyl=mint.ntheta_cyl_insitu)
            
        ##############################################
        # 2). Calculate the New Density Distribution #
        ##############################################
        # 
        # Compute new density profile
        #
        rr_cyl = qq_cyl[0]
        zrr_cyl = qq_cyl[1]
        zz_cyl = zrr_cyl * rr_cyl
        rhogas_cyl_new = np.zeros([mint.nr_cyl_insitu, mint.ntheta_cyl_insitu, nphi])

        for i in range(0, mint.nr_cyl_insitu):
            # rhogas_cyl_new[i,ntheta_cyl_insitu-1,0] = rhogas_cyl_old[i,ntheta_cyl_insitu-1,0]
            rhogas_cyl_new[i, mint.ntheta_cyl_insitu - 1, 0] = rhogas_cyl_setup[i, mint.ntheta_cyl_insitu - 1, 0]
            for j in range(mint.ntheta_cyl_insitu - 2, -1, -1):
                #
                # Recompute the density structure based on the updated gas temp
                #
                """NOTE: changed rr_cyl**2 - zz_cyl**2 to rr_cyl**2 + zz_cyl**2"""
                rad3 = (rr_cyl[i,j+1,0]**2 + zz_cyl[i,j+1,0]**2.0)**1.5 # cyl r^3
                facjp1 = const.gg*para.mstar*const.mu*zz_cyl[i,j+1,0]/(rad3*const.kk*Tgas_cyl_old[i,j+1]) # j+1
                rad3 = (rr_cyl[i,j,0]**2 + zz_cyl[i,j,0]**2.0)**1.5 # cyl r^3
                facj = const.gg*para.mstar*const.mu*zz_cyl[i,j,0]/(rad3*const.kk*Tgas_cyl_old[i,j]) # j
                fac = 0.5*(facj + facjp1) # average of [ z/h^2 over cells]
                dz = abs(zz_cyl[i,j+1,0] - zz_cyl[i,j,0])
                rhogas_cyl_new[i,j,0] = rhogas_cyl_new[i,j+1,0]*\
                                        (Tgas_cyl_old[i,j+1,0]/Tgas_cyl_old[i,j,0])*\
                                        np.exp(-fac*dz)

        # print('the rhogas in cyl coordinate')
        # print(rhogas_cyl_new)
        
        rhogas_cyl_new[np.isnan(rhogas_cyl_new)] = 0
        rhogas_cyl_new[rhogas_cyl_new <= const.mu] = const.mu
        
        # do renormalization in cylindrical coordinate to capture the real surface density
        tmpnew = np.trapz(-rhogas_cyl_setup, zz_cyl, axis=1) / np.trapz(-rhogas_cyl_new, zz_cyl, axis=1)

        rhogas_cyl_new_norm = np.zeros([mint.nr_cyl_insitu, mint.ntheta_cyl_insitu, nphi])
        for j in range(mint.nr_cyl_insitu):
            # rho_scal = rho_prime * sigma_0 / sigma_prime
            rhogas_cyl_new_norm[j,:,0] = rhogas_cyl_new[j,:,0]*tmpnew[j,0]

        # rhodust_cyl_new_norm = rhogas_cyl_new_norm.copy() / ratio_g2d # From gas to dust

        ########################################
        # 3.1). Convert back to Spherical Cord #
        ########################################
        # rhogas_sphr_new, Tg_sphr_new = cyl2sph(rc, thetac, phic, rr_cyl, zrr_cyl, rhogas_cyl_new, Tgas_cyl_old)
        # ngas_sphr_new = rhogas_sphr_new/mu

        #
        # the gas density
        #
        rhogas_sphr_new_norm, Tg_sphr_new = cyl2sph(rc, thetac, phic, rr_cyl, zrr_cyl, rhogas_cyl_new_norm, Tgas_cyl_old)
        
        # clean out the nan and inf
        # rhogas_sphr_new_norm[np.isnan(rhogas_sphr_new_norm)] = 0.0
        # rhogas_sphr_new_norm[np.isinf(rhogas_sphr_new_norm)] = 0.0
        # rhogas_sphr_new_norm[rhogas_sphr_new_norm <= const.mu] = const.mu
        
        mass_gas = (rhogas_sphr_new_norm * mint.vol).sum(0).sum(0).sum(0)
        print('gas mass after solving VHSE = %.3e'%(mass_gas * 2 / const.ms))
        
        #
        # the dust density
        #
        # if bool_dust_settling_in_iteration:
        #     """apply dust settling"""
        #     print('start solving dust settling based on the new rho gas')
        #     """now the dust is only affected by the stokes and no longer coupled with gas"""
        #     # **NOTE** (developing)
        #     # record the masses for each grain sizes
        #     rhodust = dustdens.rhodust.copy()
        #     mdust_orig_foreachsize = np.zeros(nrspec)
        #     for i_dust in range(nrspec):
        #         mdust_orig_foreachsize[i_dust] = (rhodust[:, :, :, i_dust] * mint.vol).sum(0).sum(0).sum(0)
            
        #     # convert dust to cyl coord
        #     _, _, _, rhodust_cyl_old = sph2cyl_dust(rc, thetac, phic, rhodust, nr_cyl=mint.nr_cyl_insitu, nz_cyl=mint.ntheta_cyl_insitu)
            
        #     """new function (with comments from Uma)"""
        #     # apply settling
        #     visc_alpha = para.visc_alpha
        #     rhodust_cyl_settled_norm = set_dust_settling_cyl(mint, qq_cyl, visc_alpha, rr_cyl, zrr_cyl, rhogas_cyl_old, Tgas_cyl_old, rhodust_cyl_old=rhodust_cyl_old, verbose=verbose)
            
        #     """old function
        #     # apply settling
        #     visc_alpha = para.visc_alpha
        #     rhodust_cyl_settled_norm = set_dust_settling_cyl_vert_isothermal(mint, qq_cyl, visc_alpha, nrspec, rr_cyl, zrr_cyl, rhogas_cyl_old, Tgas_cyl_old, rhodust_cyl_input=rhodust_cyl_old, verbose=verbose)
        #     """
            
        #     """2nd part: convert back to sph coordinate and reset the rhogas"""
            
        #     # convert back to sph coord
        #     rhodust_sphr_settled_norm = cyl2sph_dust(rc, thetac, phic, rr_cyl, zrr_cyl, rhodust_cyl_settled_norm)
            
        #     #
        #     # clean data on the edge
        #     #
        #     rhodust_sphr_settled_norm[np.isinf(rhodust_sphr_settled_norm)] = 0.0
        #     rhodust_sphr_settled_norm[np.isnan(rhodust_sphr_settled_norm)] = 0.0
        #     rhodust_sphr_settled_norm[rhodust_sphr_settled_norm <= 1e-199] = 1e-199
            
        #     # check mass conservation
        #     mass_t = (rhodust_sphr_settled_norm.sum(3) * mint.vol).sum(0).sum(0).sum(0)
        #     if verbose:
        #         print('dust mass after dust settling = %.3e' % (2 * mass_t / const.ms))
            
        #     #
        #     # scale again to get the correct mdiskd (if it is off)
        #     #
        #     # rhodust_sphr_settled_norm_intphi = rhodust_sphr_settled_norm_intphi * (para.mdiskd / 2 / mass_t)
        #     for i_dust in range(nrspec):
                
        #         # do renormalization in cylindrical coordinate to capture the real surface density
        #         tmpnew = np.trapz(-rhodust[:, :, :, i_dust], zz, axis=1) / np.trapz(-rhodust_sphr_settled_norm[:, :, :, i_dust], zz, axis=1)
                
        #         # Apply renormalization factor
        #         rhodust_sphr_settled_norm[:, :, :, i_dust] *= tmpnew[:, np.newaxis]
                
        #         # mdust_orig_foreachsize
        #         mdust_t_thissize = (rhodust_sphr_settled_norm[:, :, :, i_dust] * mint.vol).sum(0).sum(0).sum(0)
        #         rhodust_sphr_settled_norm[:, :, :, i_dust] = rhodust_sphr_settled_norm[:, :, :, i_dust] * (mdust_orig_foreachsize[i_dust]/mdust_t_thissize)
            
        #     #
        #     # update the gas-to-dust mass ratio
        #     #
        #     rhodust_all_sphr_settled_norm = np.nansum(rhodust_sphr_settled_norm, axis=3)
        #     rhodust_all_sphr_settled_norm[np.isnan(rhodust_all_sphr_settled_norm)] = 0.0
        #     rhodust_all_sphr_settled_norm[rhodust_all_sphr_settled_norm <= 1e-199] = 1e-199
        #     ratio_g2d_grid = rhogas_sphr_new_norm / rhodust_all_sphr_settled_norm
        #     ratio_g2d_grid[ratio_g2d_grid == np.inf] = np.nan
        #     # ratio_g2d_grid[np.isnan(ratio_g2d_grid)] = 0.0
            
        #     mint.ratio_g2d = ratio_g2d_grid
            
        #     if verbose:
        #         print('placing the new ratio_g2d_grid in vhse loop and after dust settling')
        #         print('the new ratio_g2d_grid at top inner disk is')
        #         print(ratio_g2d_grid[:10, 0, 0])
                
        #         print('the new ratio_g2d_grid at midplane inner disk is')
        #         print(ratio_g2d_grid[:10, -1, 0])
        #         # print(ratio_g2d_grid)
                
        #         print('the new ratio_g2d_grid at the first radius is')
        #         print(ratio_g2d_grid[0, -10:, 0])
                
        #     # save the dust distribution
        #     rhodust_sphr_new_norm = rhodust_sphr_settled_norm
        #     rhodust_a_t = rhodust_sphr_new_norm
        #     rhodust_all_sphr_new_norm = np.nansum(rhodust_sphr_new_norm, axis=3)
            
        #     mass_corr = (rhodust_all_sphr_new_norm * mint.vol).sum(0).sum(0).sum(0)
        #     print('dust mass after vhse and settling (normalized) = %.3e' % (mass_corr * 2 / const.ms))
                       
        #     # # NOTE: earlier version below
        #     # rhodust_sphr_new_settled_norm = cyl2sph_dust(rc, thetac, phic, rr_cyl, zrr_cyl, rhodust_cyl_settled_normed)
            
        #     # rhodust_sphr_new_settled_norm[np.isinf(rhodust_sphr_new_settled_norm)] = 0.0
        #     # rhodust_sphr_new_settled_norm[np.isnan(rhodust_sphr_new_settled_norm)] = 0.0
        #     # rhodust_sphr_new_settled_norm_intphi = rhodust_sphr_new_settled_norm.sum(3)
            
        #     # # check mass conservation
        #     # mass = (rhodust_sphr_new_settled_norm_intphi * mint.vol).sum(0).sum(0).sum(0)
        #     # print('mass after dust settling = %.3f' % (mass / const.ms))
            
        #     # rhodust_sphr_new_norm = rhodust_sphr_new_settled_norm_intphi
        
        # else:
        """
        do not apply dust settling inside function
        because now we solve dust settling before the vhse
        so even if we want to include the dust settling solver
        inside each iteration, we include those part in the 
        runmodel function instead of here
        
        but there are two cases:
        (a) we want the gas and dust are well-mixed: 
            then the rhodust is updated together with rhogas here
        (b) we decouple the gas and dust:
            then the rhodust is only set by the viscous
            this is called in the set_dust_settling function
        """
        if bool_dust_settling_in_iteration:
            """
            if we want to include dust settling
            here dust and gas are de-coupled
            the rhodust will not be updated in this function
            NOTE: the settling is not included here! 
                  settling needs to be called in set_dust_settling function
                  and the gas and dust are decoupled in this case
                  the code below is just to make sure this function 
                  does not update rhodust
            """
            print('dust and gas are de-coupled, here rhodust is not updated in VHSE solver')
            
            # this is the rhodust for all sizes summed up
            # we keep the same ratio_g2d if we do not apply dust settling here
            rhodust_a_t = dustdens.rhodust.copy()
            rhodust_all_sphr_new_norm = np.nansum(rhodust_a_t, axis=3) # rhogas_sphr_new_norm.copy() / mint.ratio_g2d 
        
            #
            # clean data on the edge
            #
            # rhodust_all_sphr_new_norm[np.isnan(rhodust_all_sphr_new_norm)] = 0
            # rhodust_all_sphr_new_norm[np.isinf(rhodust_all_sphr_new_norm)] = 0
            
            #
            # scale again to get the correct mdiskd (if it is off)
            #
            mass = (rhodust_all_sphr_new_norm * mint.vol).sum(0).sum(0).sum(0)
            print('dust mass after solving VHSE = %.3e' % (mass * 2 / const.ms))

            rhodust_all_sphr_new_norm = rhodust_all_sphr_new_norm * (para.mdiskd / 2 / mass)

            mass_corr = (rhodust_all_sphr_new_norm * mint.vol).sum(0).sum(0).sum(0)
            print('dust mass after re-normalization = %.3e' % (mass_corr * 2 / const.ms))
            
            #
            # update the gas-to-dust mass ratio
            #
            rhodust_all_sphr_new_norm[np.isnan(rhodust_all_sphr_new_norm)] = 0.0
            rhodust_all_sphr_new_norm[rhodust_all_sphr_new_norm <= const.mu * 0.01] = const.mu * 0.01
            ratio_g2d_grid = rhogas_sphr_new_norm / rhodust_all_sphr_new_norm
            ratio_g2d_grid[ratio_g2d_grid == np.inf] = np.nan
            # ratio_g2d_grid[np.isnan(ratio_g2d_grid)] = 0.0
            
            mint.ratio_g2d = ratio_g2d_grid
            
            if verbose:
                print('placing the new ratio_g2d_grid in vhse loop')
                print('the new ratio_g2d_grid at top inner disk is')
                print(ratio_g2d_grid[:10, 0, 0])
                
                print('the new ratio_g2d_grid at midplane inner disk is')
                print(ratio_g2d_grid[:10, -1, 0])
                # print(ratio_g2d_grid)
                
                print('the new ratio_g2d_grid at the first radius is')
                print(ratio_g2d_grid[0, -10:, 0])
                
        else:
            """
            if we do not include dust settling
            then gas and dust are well-mixed vertically
            """
            print('dust and gas are well-mixed vertically')
            
            # this is the rhodust for all sizes summed up
            rhodust_all_sphr_new_norm = rhogas_sphr_new_norm.copy() / mint.ratio_g2d
            # all dust sizes are well-mixed, so we do not assign specific rhodust for each size
            rhodust_a_t = None
            
            #
            # clean data on the edge
            #
            rhodust_all_sphr_new_norm[np.isnan(rhodust_all_sphr_new_norm)] = 0
            rhodust_all_sphr_new_norm[np.isinf(rhodust_all_sphr_new_norm)] = 0
            #
            # scale again to get the correct mdiskd (if it is off)
            #
            mass = (rhodust_all_sphr_new_norm * mint.vol).sum(0).sum(0).sum(0)
            print('dust mass after solving VHSE = %.3e' % (mass * 2 / const.ms))

            rhodust_all_sphr_new_norm = rhodust_all_sphr_new_norm * (para.mdiskd / 2 / mass)

            mass_corr = (rhodust_all_sphr_new_norm * mint.vol).sum(0).sum(0).sum(0)
            print('dust mass after re-normalization = %.3e' % (mass_corr * 2 / const.ms))

        #
        # then get the correct gas density with the scaled dust disk
        #
        rhogas_sphr_new_norm = rhodust_all_sphr_new_norm.copy() * mint.ratio_g2d
        rhogas_sphr_new_norm[np.isnan(rhogas_sphr_new_norm)] = 0.0
        rhogas_sphr_new_norm[rhogas_sphr_new_norm <= const.mu] = const.mu
        
        mass_gas = (rhogas_sphr_new_norm * mint.vol).sum(0).sum(0).sum(0)
        print('gas mass after solving VHSE and casted from the corrected dust mass = %.3e'%(mass_gas * 2 / const.ms))
        
        rhogas_sphr_new_norm = rhogas_sphr_new_norm  * (mint.mdiskg_setup / 2 / mass_gas)
        
        mass_gas_corr = (rhogas_sphr_new_norm * mint.vol).sum(0).sum(0).sum(0)
        print('gas mass after re-normalization = %3e'%(mass_gas_corr * 2 / const.ms))
        
        ############################################
        # 3.2). Renormalize by the surface density #
        ############################################
        """
        NOTE: the renormalization was done above in the cylindrical coordinates now
        """
        # rhogas_sphr_new_norm = np.zeros([nr,ntheta,nphi])
        # #
        # # Do the renormalization according to the original rhogas to avoid potential mass loss
        # #
        # tmpnew = np.trapz(-rhogas_sph_setup,zz,axis=1)/np.trapz(-rhogas_sphr_new,zz,axis=1)
        # # tmpnew = np.trapz(-rhogas_orig,zz,axis=1)/np.trapz(-rhogas_sphr_new,zz,axis=1)
        # for j in range(nr):
        #     rhogas_sphr_new_norm[j,:,0] = rhogas_sphr_new[j,:,0]*tmpnew[j,0] # rho_scal = rho_prime * sigma_0 / sigma_prime
        # rhodust_all_sphr_new_norm = rhogas_sphr_new_norm.copy() / ratio_g2d # From gas to dust

        #########################################################
        # 4) Save the Data and Put the dust radial distribution #
        #########################################################

        # save the rhonew for dust mass density
        rhodust_all = rhodust_all_sphr_new_norm.copy()
        # rhonew = rhodust_input.copy()
        # rhonew[np.isnan(rhonew)] = 0
        # rhonew[np.isinf(rhonew)] = 0
        # rhodust_input = rhonew.copy()
        
    else:
        """
        in case we do not need to solve the VHSE here
        we just put everything as read in
        """
        rhodust = dustdens.rhodust.copy()
        rhodust_all = np.nansum(rhodust, axis=3)
        rhogas_sphr_new_norm = rhogas.copy()
        rhodust_a_t = rhodust
    
    """
    continuue to place the dust radial distribution
    this is mainly used to apply the dust fragmentation and radial drift
    """
    rhodust_a = place_rhodust_all(mint, rhodust_all, 
                                  rhodust_a_input=rhodust_a_t,
                                  bool_dust_fragmentation=bool_dust_fragmentation_set, 
                                  bool_dust_radial_drifting=bool_dust_radial_drifting_set,
                                  tgas_input=Tgas)
    
    # **NOTE** because the dust radial drift will change the sigmad, but we want to keep the same sigmag (which should not be affected by the dust radial drifting). Therefore the gas to dust mass ratio (mint.ratio_g2d) should be updated as well.
        
    # **NOTE** (developing)
    # update the gas-to-dust mass ratio
    ratio_g2d_grid = rhogas_sphr_new_norm / np.nansum(rhodust_a, axis=3)
    ratio_g2d_grid[ratio_g2d_grid == np.inf] = np.nan
    ratio_g2d_grid[np.isnan(ratio_g2d_grid)] = 0.0
    
    mint.ratio_g2d = ratio_g2d_grid
    
    if verbose:
        print('the ratio g2d grid after this vhse iteration')
        print('the new ratio_g2d_grid at top is')
        print(ratio_g2d_grid[:10, 0, 0])
        
        print('the new ratio_g2d_grid at midplane is')
        print(ratio_g2d_grid[:10, -1, 0])
        
        print('the new ratio_g2d_grid at the first radius is')
        print(ratio_g2d_grid[0, -10:, 0])
        
    # if verbose:
    print('the gas-to-dust mass ratio is also updated, and the numbers are saved to mint.ratio_g2d and a local file xxx_ratio_g2d_grid.dat')
    save_3dgrid_to_dat(ratio_g2d_grid, 'ratio_g2d_grid.new')
    
    # clean up the rhodust_a
    # rhodust_a[np.isnan(rhodust_a)] = 0.0
    # rhodust_a[np.isinf(rhodust_a)] = 0.0
    # for ia in range(nrspec):
    #     rhodust_a_ia = rhodust_a[:, :, :, ia]
    #     rhodust_min_ia = const.mu * 0.01 * mint.fracs[ia]
    #     rhodust_a_ia[rhodust_a_ia <= rhodust_min_ia] = rhodust_min_ia
    #     rhodust_a[:, :, :, ia] = rhodust_a_ia
    
    # Clean up rhodust_a by setting a bottom for minimum values
    rhodust_a = clean_rhodust_array(mint, rhodust_a)
    with open('dust_density.new', 'w+') as f:
        f.write('1\n')  # Format number
        f.write('%d\n' % (nr * ntheta * nphi))  # Nr of cells
        f.write('%i\n' % (nrspec))  # Nr of dust species
        for i in range(nrspec):
            # data = mint.fracs[i] * rhonew.ravel(order='F')  # Create a 1-D view, fortran-style indexing
            data = rhodust_a[:,:,:,i].ravel(order='F')         # Create a 1-D view, fortran-style indexing
            data.tofile(f, sep='\n', format="%13.6e")
            f.write('\n')

    # with open('gas_temperature.new', 'w+') as f:
    #     f.write('1\n')  # Format number
    #     f.write('%d\n' % (nr * ntheta * nphi))  # Nr of cells
    #     data = Tgas.ravel(order='F')  # Create a 1-D view, fortran-style indexing
    #     data.tofile(f, sep='\n', format="%13.6e")
    #     f.write('\n')

    with open('gas_temperature.inp', 'w+') as f:
        f.write('1\n')  # Format number
        f.write('%d\n' % (nr * ntheta * nphi))  # Nr of cells
        data = Tgas.ravel(order='F')  # Create a 1-D view, fortran-style indexing
        data.tofile(f, sep='\n', format="%13.6e")
        f.write('\n')

    print('Done with Current Iteration, the new gas_temperature.inp for the dust_density.inp and dust_density.new are written')

    return 1


def fill_nan_in_2d_array_inlog(array, method='1D', method_interpolate='linear', verbose=False):
    """developing"""
    
    array = np.log10(array)
    
    # Create a meshgrid of the row and column indices
    rows, cols = np.indices(array.shape)
    
    if method == '1D':

        # Mask valid (non-NaN) values
        valid_mask = ~np.isnan(array)

        # Interpolate NaN values along the first dimension (rows for each column)
        array_filled = array.copy()
        for col in range(array.shape[1]):  # Loop over each column
            # Get the current column
            column = array[:, col]
            
            # Mask for valid (non-NaN) values in the column
            valid_mask = ~np.isnan(column)
            
            # Perform interpolation only if there are NaN values
            if np.any(~valid_mask):  # Check if there are NaN values
                # Interpolate NaN values in the column
                array_filled[:, col] = np.interp(
                    np.arange(len(column)),         # x-coordinates of all values
                    np.arange(len(column))[valid_mask],  # x-coordinates of valid values
                    column[valid_mask]              # Corresponding valid values
                )
        
    # elif method == '2D':
    else:
        # Flatten the grid and the array
        flat_rows = rows.ravel()
        flat_cols = cols.ravel()
        flat_array = array.ravel()

        # Mask for valid (non-NaN) values
        valid_mask = ~np.isnan(flat_array)

        # Extract valid points and their corresponding values
        valid_points = np.array([flat_rows[valid_mask], flat_cols[valid_mask]]).T
        valid_values = flat_array[valid_mask]

        # Points for NaN values
        nan_points = np.array([flat_rows[~valid_mask], flat_cols[~valid_mask]]).T

        # Perform 2D interpolation using griddata
        filled_values = interpolate.griddata(valid_points, valid_values, nan_points, method=method_interpolate)

        # Fill the NaN positions in the original array
        array_filled = array.copy()
        array_filled[np.isnan(array)] = filled_values
    
    if verbose:
        print("Original array:")
        print(array)
        print("\nArray with NaN values filled:")
        print(array_filled)
    
    return 10**array_filled


def set_dust_settling_cyl(mint, qq_cyl, visc_alpha, rr_cyl, zrr_cyl, rhogas_cyl_old, Tgas_cyl_old, rhodust_cyl_old=[], rhodust_cyl_probset=None, verbose=False):
    """ 
    **NOTE** (Developing)
    Solve the dust settling in the cylindrical coordinate
    This will be a subroutine in the get new dust density since it
    needs to be done in the cylindrical coordinate
    
    Args:
        mint (object): The DiskMINT model object containing the parameters and data.
        qq_cyl: 2d meshgrid consisting with information of radius (r) and height/radius (z/r) in cylindrical coordinates
        visc_alpha: the parameter describing the efficiency of the dust settling
        verbose (bool): whether to print out for debugging
        
    Returns:
        rhodust_cyl_new_settled_norm: the array with the shape of (nr, nz, ndust), describes the rhodust in cylindrical coordinate of how much dust is settled
    """
    
    para = mint.parameters
    
    rr_cyl  = qq_cyl[0]
    zrr_cyl = qq_cyl[1]
    zz_cyl  = zrr_cyl * rr_cyl
    
    # 
    # Apply Settling
    #
    
    # Defining variables I think you already have in your code:
    #....................................................
    # rad : Radial distance defined from r, z cylindrical coordinates
    # omg : Keplerian angular velocity
    # cs  : Isothermal sound speed at this (r,z) from gas temperature tg
    # hp  : Gas pressure scale height at this (r,z)

    # ALSO Need rhogas(r,z), Sigmad(a,r), particle density rhogr(a)
    #....................................................
    rad = np.sqrt(rr_cyl**2 + zz_cyl**2)
    omg = np.sqrt(const.GG * para.mstar / (rad**3))
    cs  = np.sqrt(const.kk * Tgas_cyl_old / const.mu)
    hp  = cs/omg
    
    alpha = visc_alpha
    
    #....................................................
    # This below needs to be in loop to do this for each particle size
    #....................................................
    # define the output rhodust
    rhodust_cyl_settled_norm = np.ones(rhodust_cyl_old.shape) * np.nan
    # define delta_zz
    delta_zz_cyl = np.ones(zz_cyl.shape) * np.nan
    delta_zz_cyl[:, 1:, :] = zz_cyl[:, :-1, :] - zz_cyl[:, 1:, :]
    delta_zz_cyl[:, 0, :]  = delta_zz_cyl[:, 1, :]
    # loop through all sizes
    for i_dust in range(mint.dust_spec_nr):
        adust = mint.a_ave[i_dust] # Dust size in cm 

        if verbose:
            print('solving the dust settling for dust #%i with size %.2e[cm]'%(i_dust, adust))
        
        rhogr = para.rhobulk # the bulk density of the grain
        # Define Stokes for this particle size
        st_z = rhogr * adust/( rhogas_cyl_old * hp ) # function of (r,z)

        # Integrating Eq 18 of Fromang & Nelson 2009
        # d_z[ ln (rho_d(a)/rhogas) ] = - omg^2 tau_s z/D = - Stokes/alpha * (z/hp^2)

        # Define delta_z = z[j+1] - z[j] in au and so on for the grid if you don't have it
        # delta_z = delta_zz_cyl
        
        # fractional mass in bin, check for factor of two for both sides of disk if not in grid
        # this needs to be done at each radii
        # fracd = sigmad / np.sum(rhogas_cyl_old * delta_zz_cyl)
        sigmad_idust = np.nansum(rhodust_cyl_old[:, :, :, i_dust] * delta_zz_cyl, axis=1)
        fracd    = sigmad_idust / np.nansum(rhogas_cyl_old * delta_zz_cyl, axis=1)
        # Define the quantity ln rho_d(a)/rhogas which is ln of d2g ratio for this particle size
        # ln_d2ga = np.log(fracd) * np.ones(np.size(z)) # Default value set to well mixed ratio
    
        # ln_d2ga = np.ones(zz_cyl.shape)
        # for i_r in range(zz_cyl.shape[0]):
        #     ln_d2ga[i_r, :, 0] = np.log(fracd[i_r, 0]) * np.ones(zz_cyl.shape[1])
        # the above can be simplified as below
        ln_d2ga = np.log(fracd[:, 0])[:, np.newaxis, np.newaxis] * np.ones_like(zz_cyl)
        
        # # Solve the PDE at different height j (solve from mid-plane)
        # for j_z in np.arange(zz_cyl.shape[1]-2, -1, -1):
        #     dum_x = zz_cyl[:, j_z+1, 0]/hp[:, -1, 0]
        #     for i_r in range(zz_cyl.shape[0]):
        #         if dum_x[i_r] <= 0.33:
        #             ln_d2ga[i_r, j_z, 0]=ln_d2ga[i_r, j_z+1, 0]
        #         else:
        #             ln_d2ga[i_r, j_z, 0]=ln_d2ga[i_r, j_z+1, 0] - delta_zz_cyl[i_r, j_z+1, 0] * (st_z[i_r, j_z+1, 0]/alpha) * zz_cyl[i_r, j_z+1, 0]/hp[i_r, j_z+1, 0]**2

        # Solve the PDE at different heights j (solve from mid-plane)
        for j_z in range(zz_cyl.shape[1] - 2, -1, -1):
            dum_x = zz_cyl[:, j_z + 1, 0] / hp[:, -1, 0]
            mask = dum_x <= 0.33  # Boolean mask for the condition
                                  # we only solve the settling to a certain scale height
                                  # because the vertical turbulent will still bring up 
                                  # the dust that is very close to the midplane. 
                                  # this prevents the dust to collapse into a infinitely small layer
            ln_d2ga[:, j_z, 0] = ln_d2ga[:, j_z + 1, 0]  # Default assignment
            ln_d2ga[~mask, j_z, 0] -= (
                delta_zz_cyl[~mask, j_z + 1, 0]
                * (st_z[~mask, j_z + 1, 0] / alpha)
                * zz_cyl[~mask, j_z + 1, 0]
                / hp[~mask, j_z + 1, 0] ** 2
            )
        
        rhodust_cyl_settled_idust =  np.exp(ln_d2ga) * rhogas_cyl_old # convert ratio to density
        
        # clean up the nan and inf small values
        # rhodust_cyl_settled_idust[np.isnan(rhodust_cyl_settled_idust)] = 0.0
        # rhodust_cyl_settled_idust[np.isinf(rhodust_cyl_settled_idust)] = 0.0
        rhodust_cyl_settled_idust[~np.isfinite(rhodust_cyl_settled_idust)] = 0.0
        rhodust_min_idust = const.mu * 0.01 * mint.fracs[i_dust]
        rhodust_cyl_settled_idust[rhodust_cyl_settled_idust <= rhodust_min_idust] = rhodust_min_idust
        
        # Renormalize to conserve mass in this bin
        # for i_r in range(zz_cyl.shape[0]):
        #     rhodust_cyl_settled_norm[i_r, :, 0, i_dust] = rhodust_cyl_settled[i_r, :, 0] * (sigmad_idust/np.nansum(rhodust_cyl_settled * delta_zz_cyl, axis=1)) # check for factor of two here again!
        if rhodust_cyl_probset is None:
            sigmad_idust_norm = sigmad_idust
        else:
            sigmad_idust_norm = np.nansum(rhodust_cyl_probset[:, :, :, i_dust] * delta_zz_cyl, axis=1)
            
        norm_factor = sigmad_idust_norm / np.nansum(rhodust_cyl_settled_idust * delta_zz_cyl, axis=1)
        rhodust_cyl_settled_norm[:, :, 0, i_dust] = rhodust_cyl_settled_idust[:, :, 0] * norm_factor[:, np.newaxis, 0]
    
    return rhodust_cyl_settled_norm


def set_dust_settling_cyl_vert_isothermal(mint, qq_cyl, visc_alpha, nrspec, rr_cyl, zrr_cyl, rhogas_cyl_new_norm, Tgas_cyl_old, rhodust_cyl_input=[], rhodust_cyl_probset=None, verbose=False):
    """ 
    **NOTE** (Developing)
    Solve the dust settling in the cylindrical coordinate
    This will be a subroutine in the get new dust density since it
    needs to be done in the cylindrical coordinate
    
    This method assume that the disk is vertically isothermal
    
    Args:
        mint (object): The DiskMINT model object containing the parameters and data.
        qq_cyl: 2d meshgrid consisting with information of radius (r) and height/radius (z/r) in cylindrical coordinates
        visc_alpha: the parameter describing the efficiency of the dust settling
        verbose (bool): whether to print out for debugging
        
    Returns:
        rhodust_cyl_new_settled_norm: the array with the shape of (nr, nz, ndust), describes the rhodust in cylindrical coordinate of how much dust is settled
    """
    
    para = mint.parameters
    
    rr_cyl  = qq_cyl[0]
    zrr_cyl = qq_cyl[1]
    zz_cyl  = zrr_cyl * rr_cyl

    #
    # 1. calculate the pre factors for the dust grains.
    #
    # find out the mid-plane scale height (this is the minimum settling z)
    # Tgas_cyl_midplane_old = np.zeros(Tgas_cyl_old.shape)
    # Tgas_cyl_midplane_old[:,:,0] = np.repeat(Tgas_cyl_old[:,-1,:],int(Tgas_cyl_old.shape[1]),axis=1)
    hp = np.sqrt(const.kk / (const.mu*const.GG*para.mstar)) \
    * np.sqrt(Tgas_cyl_old[:,:,:] * (rr_cyl**2+zz_cyl**2)**(3/2))
    # print('hp shape ', hp.shape)
    
    hp_midplane = np.zeros(hp.shape)
    hp_midplane[:,:,0] = np.repeat(hp[:, -1, :], int(hp.shape[1]), axis = 1)
    
    the_shape = Tgas_cyl_old.shape
    if rhodust_cyl_input == []:
        rhodust_cyl = np.ones([the_shape[0], the_shape[1], the_shape[2], nrspec])
    else:
        rhodust_cyl = rhodust_cyl_input
    
    fac_rhogas_int = np.trapz(-rhogas_cyl_new_norm[:, :, 0], zz_cyl[:, :, 0], axis=1)
    
    stokes = np.zeros(rhodust_cyl.shape)
    zmax_all = np.zeros(rhodust_cyl.shape)
    # zmax_all_notwithhp = np.zeros(rhodust_cyl.shape)
    d2g = np.zeros(rhodust_cyl.shape)
    d2g_settled = np.zeros(rhodust_cyl.shape)
    rhodust_cyl_non_settled = np.zeros(rhodust_cyl.shape)
    rhodust_cyl_settled = np.zeros(rhodust_cyl.shape)
    rhodust_cyl_settled_normed = np.zeros(rhodust_cyl.shape)
    
    rhobulk = np.ones(nrspec) * para.rhobulk
    
    for i_dust in range(nrspec):
    # for i_dust in [0, 1]:
        print('calculate stokes number for species # %i' % i_dust)
    
        #
        # 2. calculate the Stokes numbers
        # here the a_ave is in the unit of cm
        #
        fac1 = mint.a_ave[i_dust] * rhobulk[i_dust] # rhodust[:, :, 0, i_dust]
        fac2 = hp * rhogas_cyl_new_norm[:, :, :] # * 2 * np.pi
        stokes_i = fac1/fac2
        stokes[:, :, :, i_dust] = stokes_i
        
        # if verbose:
        #     print('stokes for dust #%i'%(i_dust))
        #     print(stokes[0:20, -1, 0, i_dust])
        #     # print(stokes)
        #     # print(stokes.shape)

        # set the zmax accordingly for each dust grain
        # ref. Estrada+2016
        # zmax_i = hp/ \
        #          np.sqrt(1.0 \
        #          + stokes_i*(1.0 + stokes_i**2)/visc_alpha)
        zmax_i = hp/\
                 np.sqrt(1.0 + stokes_i/visc_alpha)
        
        # if verbose:
        #     print('zmax for dust #%i'%(i_dust))
        #     print(zmax_i[0:20, -1, 0])
        #     # print(zmax_i)

        zmax_i[np.isnan(zmax_i)] = 0.0

        # print('zmax_i shape ', zmax_i.shape)
        for i_coord in range(zmax_i.shape[0]):
            for j_coord in range(zmax_i.shape[1]):
                i_coord, j_coord = int(i_coord), int(j_coord)

                # print(hp[i_coord, j_coord])
                # print(zmax_i[i_coord, j_coord])

                # zmax_all[i_coord, j_coord, 0, i_dust] = \
                #     np.max([hp_midplane[i_coord, j_coord, 0], zmax_i[i_coord, j_coord, 0]])
                """
                **NOTE**
                this is to prevent the large grains collapse to inf thin layer
                """
                zmax_all[i_coord, j_coord, 0, i_dust] = \
                    np.max([hp_midplane[i_coord, j_coord, 0], zmax_i[i_coord, j_coord, 0]])
                
                # zmax_all_notwithhp[i_coord, j_coord, 0, i_dust] =\
                #     zmax_i[i_coord, j_coord, 0]
        
        # if verbose:
        #     print('hpmid for dust #%i'%(i_dust))
        #     print(hp_midplane[0:20, -1, 0])
            
        #     print('max(hpmid, zmax) for dust #%i'%(i_dust))
        #     print(zmax_all[0:20, -1, 0, i_dust])
        
        # set up the d2g before the settling
        """
        **NOTE**
        Here we distribute the d2g from a constant value everytime
        and re-distribute the dust grains to new positions accordingly.
        
        (we ignored the original d2g calculated from the last iteration)
        """
        # d2g[:, :, 0, i_dust] = rhodust[:, :, 0, i_dust] / rhogas_cyl_new_norm[:, :, 0]
        if rhodust_cyl_input == []:
            rhodust_cyl_non_settled[:, :, :, i_dust] = mint.fracs[i_dust] * (rhogas_cyl_new_norm[:, :, :]/mint.ratio_g2d_global)
        else:
            rhodust_cyl_non_settled[:, :, :, i_dust] = rhodust_cyl_input[:, :, :, i_dust]
        
        # d2g[:, :, :, i_dust] = rhodust_cyl_non_settled[:, :, :, i_dust]/rhogas_cyl_new_norm[:, :, :]
        
        # if verbose:
        #     print('rhodust_cyl_non_settled for dust #%i'%(i_dust))
        #     print(rhodust_cyl_non_settled[0:20, -1, 0, i_dust])
            
        #     print('rhogas_cyl_new_norm')
        #     print(rhogas_cyl_new_norm[0:20, -1, 0])
        
        #     print('d2g for dust #%i'%(i_dust))
        #     print(d2g[0:20, -1, 0, i_dust])
        
        #
        # 3. distribute the dust grains for different sizes
        #

    # for i_dust in range(nrspec):
        # calculate new eps (d2g ratio)
        # d2g = rhodust / rhogas
        
        """the previous version"""
        # d2g_settled[:, :, :, i_dust] = d2g[:, :, :, i_dust] \
        #                                * np.exp(-0.5 * (zz_cyl**2 / hp**2) \
        #                                * ((hp/zmax_all[:, :, :, i_dust])**2.0 - 1.0))

        """the new version"""
        fac_rhodust_int = np.trapz(-rhodust_cyl_non_settled[:, :, 0, i_dust], zz_cyl[:, :, 0], axis=1)
        f = fac_rhodust_int/fac_rhogas_int # sigmadust[i_dust] / sigmagas
        # for i_r in range(d2g_settled.shape[0]):
        #     d2g_settled[i_r, :, 0, i_dust] = f[i_r] * (( 1.0 - hp/zmax_all[i_r, :, 0, i_dust])**2.0)
        d2g_settled[:, :, 0, i_dust] = f[:, np.newaxis] * (1.0 - hp[:, :, 0] / zmax_all[:, :, 0, i_dust])**2.0
        
        # calculate new dust mass density = rho_gas * d2g
        rhodust_cyl_settled[:, :, :, i_dust] = \
            d2g_settled[:, :, :, i_dust] * rhogas_cyl_new_norm[:, :, :]
        rhodust_cyl_settled[np.isnan(rhodust_cyl_settled)] = 0.
        
        # if verbose:
        #     print('d2g_settled for dust #%i'%(i_dust))
        #     print(d2g_settled[0:20, -1, 0, i_dust])
            
        #     print('rhodust_cyl_settled for dust #%i'%(i_dust))
        #     print(rhodust_cyl_settled[0:20, -1, 0, i_dust])
        
        #
        # 4. rescale the density distribution
        #
        # Here I deciede to ignore the original d2g calculated from the last iteration
        # so the dust settling is just based on the d2g in each iteration
    
    # for i_dust in range(nrspec):
    
        # do renormalization in cylindrical coordinate to capture the real surface density
        # tmpnew = np.trapz(-rhogas_cyl_setup, zz_cyl, axis=1) / np.trapz(-rhogas_cyl_new, zz_cyl, axis=1)
        
        # Clean up the NaNs
        rhodust_cyl_non_settled_clean = clean_rhodust_array(mint, rhodust_cyl_non_settled_clean)
        rhodust_cyl_settled_clean = clean_rhodust_array(mint, rhodust_cyl_settled_clean)
        
        # do the integral to normalize the data
        if rhodust_cyl_probset is None:
            fac1 = np.trapz(-rhodust_cyl_non_settled_clean[:, :, 0, i_dust], zz_cyl[:, :, 0], axis=1)
        else:
            fac1 = np.trapz(-rhodust_cyl_probset[:, :, 0, i_dust], zz_cyl[:, :, 0], axis=1)
        fac2 = np.trapz(-rhodust_cyl_settled_clean[:, :, 0, i_dust], zz_cyl[:, :, 0], axis=1)
        tmpnew_dust = fac1 / fac2
        tmpnew_dust[np.isinf(tmpnew_dust)] = 1.0
        
        # if verbose:
        #     print('tmpnew_dust for dust #%i'%(i_dust))
        #     print(tmpnew_dust[:20])
        
        # rhogas_cyl_new_norm = np.zeros([mint.nr_cyl_insitu, mint.ntheta_cyl_insitu, nphi])
        for i_r in range(mint.nr_cyl_insitu):
            # rho_scal = rho_prime * sigma_0 / sigma_prime
            # rhogas_cyl_new_norm[i_r,:,0] = rhogas_cyl_new[i_r,:,0]*tmpnew[i_r,0]
            rhodust_cyl_settled_normed[i_r, :, 0, i_dust] = rhodust_cyl_settled[i_r, :, 0, i_dust]*tmpnew_dust[i_r]
            
        for i_dust in range(nrspec):
            rhodust_cyl_settled_normed[:, :, 0, i_dust] = fill_nan_in_2d_array_inlog(rhodust_cyl_settled_normed[:, :, 0, i_dust])
    
    return rhodust_cyl_settled_normed


def read_in_rhodust(file_name = 'rhodust_probset.inp'):
    # Open the file and read the header
    print('read in %s'%(file_name))
    with open(file_name, 'r') as f:
        # Read header lines
        format_number = int(f.readline().strip())  # Format number (e.g., 1)
        nr = int(f.readline().strip())             # Number of r cells
        ntheta = int(f.readline().strip())         # Number of theta cells
        nphi = int(f.readline().strip())           # Number of phi cells
        nrspec = int(f.readline().strip())         # Number of dust species

        # Read the flattened data from the file
        flattened_data = np.loadtxt(f, dtype=np.float64)

    # Reshape the flattened data back into a 4D array.
    # We use order='F' to match the Fortran-order flattening
    rhodust = flattened_data.reshape((nr, ntheta, nphi, nrspec), order='F')

    # Optional: print the header info and shape of the resulting array
    print("Format number:", format_number)
    print("Dimensions: nr =", nr, ", ntheta =", ntheta, ", nphi =", nphi, ", nrspec =", nrspec)
    print("Reshaped array shape:", rhodust.shape)
    
    return rhodust


def set_dust_settling(mint, verbose=False):
    """
    **NOTE** (Developing)
    Aims to have a function that could settle the dust
    from reading the radmc3dPy datafile and then output
    as a stand-alone package
    
    setting the dust settling. Introducing the setlling to the dust species.

    Args:
        mint (object): The DiskMINT model object containing the parameters and data.
        verbose (bool): set true to print out texts for debugging
        
    Returns:
        None.
    """
    
    para = mint.parameters

    # print the current working dir
    if verbose:
        print('the current working dir for set dust settling')
        print(os.getcwd())

    #
    # Read in the rhodust_probset
    # this is used to do the re-normalization to make sure the sigma dust is conserved
    #
    rhodust_probset = read_in_rhodust()
    
    #
    # Read dust_density.inp
    #
    # dustdens = read_dustdens()
    # dustdens = radmc3dPy.analyze.readData(dtemp=True, ddens=True, gtemp=True, binary=False)
    dustdens = modelgrid.readData(dtemp=True, ddens=True, gtemp=True) # , binary=False)
    # nrspec = int(dustdens.nrspec)
    nrspec = int(dustdens.rhodust.shape[-1])
    # get the rhodust
    rhodust = dustdens.rhodust.copy()
    # Use this new defination for the rhogas
    rhogas = mint.ratio_g2d * np.sum(dustdens.rhodust, axis=3)
    rhogas[np.isnan(rhogas)] = 1e-199
    rhogas[rhogas <= const.mu] = const.mu
    
    # rhogas = dustdens.rhodust[:,:,:,0]/fracs[0] * ratio_g2d
    # ngas = rhogas / const.mu

    # grid info
    nr = dustdens.grid.nx
    ntheta = dustdens.grid.ny
    nphi = dustdens.grid.nz
    rc = dustdens.grid.x
    thetac = dustdens.grid.y
    phic = dustdens.grid.z
    qq = np.meshgrid(rc, thetac, phic, indexing='ij')
    rr = qq[0]  # Spherical R
    tt = qq[1]  # Angle with z axis ( = pi/2 for the z=0 axis)
    zr = np.pi / 2.e0 - qq[1]
    zz = np.sin(zr) * rr

    # Tgas and Tdust
    Tgas = dustdens.gastemp.copy()
    Tgas = Tgas[:, :, :, 0]
    Tdust = dustdens.dusttemp.copy()
    if verbose:
        print('shape for Tgas: ', Tgas.shape)
        print('shape for Tdust: ', Tdust.shape)

    # record the masses for each grain sizes
    mdust_orig_foreachsize = np.zeros(nrspec)
    for i_dust in range(nrspec):
        mdust_orig_foreachsize[i_dust] = (rhodust[:, :, :, i_dust] * mint.vol).sum(0).sum(0).sum(0)
    
    if verbose:
        mass_orig = (rhodust.sum(3) * mint.vol).sum(0).sum(0).sum(0)
        print('input dust mass before dust settling = %.3e' % (2 * mass_orig / const.ms))
        
        rhogas_t = rhodust.sum(3) * mint.ratio_g2d
        rhogas_t[np.isnan(rhogas_t)] = 1e-199
        rhogas_t[rhogas_t <= const.mu] = const.mu
        mgas_orig = (rhogas_t * mint.vol).sum(0).sum(0).sum(0) * 2
        print('input gas mass before dust settling = %.3e'%(mgas_orig/const.ms))
    
    #
    # 0. convert from sph to cyl coord
    #
    rhogas[rhogas <= const.mu] = const.mu
    _, _, _, rhogas_cyl_old, Tgas_cyl_old = \
        sph2cyl(rc, thetac, phic, rhogas, Tgas, nr_cyl=mint.nr_cyl_insitu, nz_cyl=mint.ntheta_cyl_insitu)
    
    _, _, qq_cyl, rhodust_cyl_old = \
        sph2cyl_dust(rc, thetac, phic, rhodust, nr_cyl=mint.nr_cyl_insitu, nz_cyl=mint.ntheta_cyl_insitu)
    
    _, _, _, rhodust_cyl_probset = \
        sph2cyl_dust(rc, thetac, phic, rhodust_probset, nr_cyl=mint.nr_cyl_insitu, nz_cyl=mint.ntheta_cyl_insitu)
    
    # rhodust_cyl_old[np.isnan(rhodust_cyl_old)] = 1e-199
    # rhodust_cyl_old[rhodust_cyl_old == 0] = 1e-199

    rhodust_cyl_old = clean_rhodust_array(mint, rhodust_cyl_old)
    rhodust_cyl_probset = clean_rhodust_array(mint, rhodust_cyl_probset)
    
    rr_cyl  = qq_cyl[0]
    zrr_cyl = qq_cyl[1]
    # zz_cyl  = zrr_cyl * rr_cyl
     
    """new function (with comments from Uma)"""
    visc_alpha = para.visc_alpha 
    # rhogas_cyl_new_norm = rhogas_cyl_old.copy() 
    rhodust_cyl_settled_norm = set_dust_settling_cyl(mint, qq_cyl, visc_alpha, rr_cyl, zrr_cyl, rhogas_cyl_old, Tgas_cyl_old, rhodust_cyl_old=rhodust_cyl_old, rhodust_cyl_probset=rhodust_cyl_probset, verbose=verbose)
    
    """the method assuming vertically isothermal"""
    # visc_alpha = para.visc_alpha 
    # rhodust_cyl_settled_norm = set_dust_settling_cyl_vert_isothermal(mint, qq_cyl, visc_alpha, nrspec, rr_cyl, zrr_cyl, rhogas_cyl_old, Tgas_cyl_old, rhodust_cyl_input=rhodust_cyl_old, verbose=verbose)

    rhodust_sphr_settled_norm = cyl2sph_dust(rc, thetac, phic, rr_cyl, zrr_cyl, rhodust_cyl_settled_norm)
    
    # rhodust_sphr_settled_norm[np.isinf(rhodust_sphr_settled_norm)] = 0.0
    # rhodust_sphr_settled_norm[np.isnan(rhodust_sphr_settled_norm)] = 0.0
    rhodust_sphr_settled_norm = clean_rhodust_array(mint, rhodust_sphr_settled_norm)
    
    # renormalization in sph again
    # the mass in rhodust_sph simply from the cyl is not conserved from the problem setup
    for i_dust in range(nrspec):
        rhodust_probset_idust_t = rhodust_probset[:, :, :, i_dust]
        fac1 = np.array([np.trapz(rhodust_probset_idust_t[i_r, :, 0], -zz[i_r, :, 0]) for i_r in range(nr)])
        rhodust_settled_idust_t = rhodust_sphr_settled_norm[:, :, :, i_dust]
        fac2 = np.array([np.trapz(rhodust_settled_idust_t[i_r, :, 0], -zz[i_r, :, 0]) for i_r in range(nr)])
        norm_factor = fac1/fac2
        rhodust_sphr_settled_norm[:, :, 0, i_dust] = rhodust_sphr_settled_norm[:, :, 0, i_dust] * norm_factor[:, np.newaxis]
        
    # check mass conservation
    mass_t = (rhodust_sphr_settled_norm.sum(3) * mint.vol).sum(0).sum(0).sum(0)
    if verbose:
        print('dust mass after dust settling = %.3e' % (2 * mass_t / const.ms))

    #
    # scale again to get the correct mdiskd (if it is off)
    #
    # rhodust_sphr_settled_norm_intphi = rhodust_sphr_settled_norm.sum(2)
    # rhodust_sphr_settled_norm_intphi = rhodust_sphr_settled_norm_intphi * (para.mdiskd / 2 / mass_t)
    for i_dust in range(nrspec):
        # mdust_orig_foreachsize
        mdust_t_thissize = (rhodust_sphr_settled_norm[:, :, :, i_dust] * mint.vol).sum(0).sum(0).sum(0)
        rhodust_sphr_settled_norm[:, :, :, i_dust] = rhodust_sphr_settled_norm[:, :, :, i_dust] * (mdust_orig_foreachsize[i_dust]/mdust_t_thissize)

    # **NOTE** (developing)
    # update the gas-to-dust mass ratio
    ratio_g2d_grid = rhogas / np.nansum(rhodust_sphr_settled_norm, axis=3) # rhodust_a.sum(3)
    ratio_g2d_grid[ratio_g2d_grid == np.inf] = np.nan
    # ratio_g2d_grid[np.isnan(ratio_g2d_grid)] = 0.0
    
    mint.ratio_g2d = ratio_g2d_grid
    
    print('the new ratio_g2d_grid at top inner disk is')
    print(ratio_g2d_grid[:10, 0, 0])
    
    print('the new ratio_g2d_grid at midplane inner disk is')
    print(ratio_g2d_grid[:10, -1, 0])
    # print(ratio_g2d_grid)
    
    #
    # output
    #
    # if verbose:
    print('the gas-to-dust mass ratio is also updated, and the numbers are saved to mint.ratio_g2d and a local file xxx_ratio_g2d_grid.dat')
    save_3dgrid_to_dat(ratio_g2d_grid, 'ratio_g2d_grid.new')
    
    if verbose:
        mass_corr = (rhodust_sphr_settled_norm.sum(3) * mint.vol).sum(0).sum(0).sum(0)
        print('dust mass after dust settling and re-normalization = %.3e' % (2 * mass_corr / const.ms))

        rhogas_t = np.nansum(rhodust_sphr_settled_norm, axis=3) * mint.ratio_g2d
        rhogas_t[np.isnan(rhogas_t)] = 1e-199
        rhogas_t[rhogas_t <= const.mu] = const.mu
        mgas_out = (rhogas_t * mint.vol).sum(0).sum(0).sum(0) * 2
        print('output gas mass after dust settling and re-normalization = %.3e'%(mgas_out/const.ms))
        
        print('the current working dir for adding dust rim')
        print(os.getcwd())
    
    # save the rhonew for dust mass density
    # rhonew = rhodust_sphr_settled_norm_intphi.copy()
    # rhonew = rhodust_sphr_settled_norm
    # rhonew[np.isnan(rhonew)] = 0
    # rhonew[np.isinf(rhonew)] = 0

    # rhodust_input = rhonew.copy()
    # rhodust_a = place_rhodust_all(mint, rhodust_input)
    # rhodust_a = place_rhodust_all(mint, rhodust_input, bool_dust_fragmentation=mint.bool_dust_fragmentation, 
    # bool_dust_radial_drifting=mint.bool_dust_radial_drifting,
    # tgas_input=Tgas)
    rhodust_a = np.ones(rhodust.shape) * np.nan
    rhodust_a[:, :, :, :] = rhodust_sphr_settled_norm
    with open('dust_density_settled.new', 'w+') as f:
        f.write('1\n')  # Format number
        f.write('%d\n' % (nr * ntheta * nphi))  # Nr of cells
        f.write('%i\n' % (nrspec))  # Nr of dust species
        for i in range(nrspec):
            # data = mint.fracs[i] * rhonew.ravel(order='F')  # Create a 1-D view, fortran-style indexing
            data = rhodust_a[:,:,:,i].ravel(order='F')         # Create a 1-D view, fortran-style indexing
            data.tofile(f, sep='\n', format="%13.6e")
            f.write('\n')

    print('Done with dust settling, the new dust_density_settled.new are written')
    
    return 1


def add_dust_rim(mint, verbose=False):
    """
    **NOTE** (Developing)
    Aims to have a function that could add the dust rim
    from reading the radmc3dPy datafile and then output
    as a stand-alone package

    Args:
        mint (object): The DiskMINT model object containing the parameters and data.
        verbose (bool): set true to print out texts for debugging
        
    Returns:
        None.
    """
    
    para = mint.parameters

    # print the current working dir
    if verbose:
        print('the current working dir for set dust settling')
        print(os.getcwd())
    
    #
    # Read dust_density.inp
    #
    # dustdens = read_dustdens()
    # dustdens = radmc3dPy.analyze.readData(dtemp=True, ddens=True, gtemp=True, binary=False)
    dustdens = modelgrid.readData(dtemp=True, ddens=True, gtemp=True) # , binary=False)
    # nrspec = int(dustdens.nrspec)
    nrspec = int(dustdens.rhodust.shape[-1])
    # get the rhodust
    rhodust = dustdens.rhodust.copy()
    # Use this new defination for the rhogas
    rhogas = mint.ratio_g2d * np.sum(dustdens.rhodust, axis=3)
    rhogas[np.isnan(rhogas)] = 1e-199
    rhogas[rhogas <= const.mu] = const.mu
    
    # rhogas = dustdens.rhodust[:,:,:,0]/fracs[0] * ratio_g2d
    # ngas = rhogas / const.mu

    # grid info
    nr = dustdens.grid.nx
    ntheta = dustdens.grid.ny
    nphi = dustdens.grid.nz
    rc = dustdens.grid.x
    thetac = dustdens.grid.y
    phic = dustdens.grid.z
    qq = np.meshgrid(rc, thetac, phic, indexing='ij')
    rr = qq[0]  # Spherical R
    tt = qq[1]  # Angle with z axis ( = pi/2 for the z=0 axis)
    zr = np.pi / 2.e0 - qq[1]
    zz = np.sin(zr) * rr

    # Tgas and Tdust
    Tgas = dustdens.gastemp.copy()
    Tgas = Tgas[:, :, :, 0]
    Tdust = dustdens.dusttemp.copy()
    if verbose:
        print('shape for Tgas: ', Tgas.shape)
        print('shape for Tdust: ', Tdust.shape)

    """
    # record the masses for each grain sizes
    mdust_orig_foreachsize = np.zeros(nrspec)
    for i_dust in range(nrspec):
        mdust_orig_foreachsize[i_dust] = (rhodust[:, :, :, i_dust] * mint.vol).sum(0).sum(0).sum(0)
    """
    
    # find out the rwall midplane and rwall surface
    # also need to define this factor in the class
    # r_wall_mid  = rc[0] # no need to define mid because it is the rin.
    r_wall_surf = para.r_wall_surf # rc[5]

    # choose the select range based on r
    select_range = rc <= r_wall_surf
    
    """
    # we simply do this in the spherical coordinates, 
    # because in the inner disk the z are very close to theta
    """
    
    # get the mass before adding rim
    rhodust_all = rhodust.sum(3).copy()
    # mass_orig_inside = (rhodust_all[select_range, :, 0] * mint.vol[select_range, :, 0]).sum(0).sum(0) * 2.0

    #
    # read the set up rhodust
    #
    rhodust_all0 = np.zeros(rhodust_all.shape)
    rhodust_all0[:, :, 0] = np.loadtxt('rhodust_all_probset.inp')
    
    if verbose:
        mdust_probset = (rhodust_all0 * mint.vol).sum(0).sum(0).sum(0) * 2
        print('input dust mass as the problem setup = %.3e'%(mdust_probset / const.ms))
        
        mdust_orig = (rhodust_all * mint.vol).sum(0).sum(0).sum(0) * 2
        print('input dust mass before adding a rim = %.3e'%(mdust_orig / const.ms))
        
        rhogas_t = rhodust_all * mint.ratio_g2d
        rhogas_t[np.isnan(rhogas_t)] = 1e-199
        rhogas_t[rhogas_t <= const.mu] = const.mu
        mgas_orig = (rhogas_t * mint.vol).sum(0).sum(0).sum(0) * 2
        print('input gas mass before adding a rim = %.3e'%(mgas_orig/const.ms))
    
    # calculate the hp
    hp = np.sqrt(const.kk / (const.mu*const.GG*para.mstar)) \
        * np.sqrt(Tgas[:,:,:] * (rr**3))

    # get the new zwall based on the wall factor xi
    # need to define this factor in the class
    factor_wall_height = para.factor_wall_height # 4.0
    zmax_wall = hp * factor_wall_height

    # redistributie the masses in the wall range
    rhodust_addwall = rhodust.copy()

    # according to the previous hp and the new zwall
    for i_dust in range(int(rhodust.shape[-1])):
        rhodust_addwall[select_range, :, 0, i_dust] = rhodust[select_range, :, 0, i_dust] \
            * np.exp(-0.5 * ((zz[select_range, :, 0])**2 / (hp[select_range, :, 0])**2) \
            * ((hp[select_range, :, 0]/zmax_wall[select_range, :, 0])**2.0 - 1.0))

    rhodust_addwall[np.isnan(rhodust_addwall)] = 0.0
    rhodust_addwall[np.isinf(rhodust_addwall)] = 0.0
    
    """
    # simply distributing using Gaussian distribution
    rhodust_addwall_all_nonorm = 1/(np.sqrt(2*np.pi) * zmax_wall) * np.exp(-((zz)**2 / ( 2 * zmax_wall**2)))

    for i_dust in range(int(rhodust.shape[-1])):
        rhodust_addwall[select_range, :, :, i_dust] = rhodust_addwall_all_nonorm[select_range, :, :] * mint.fracs[i_dust]
    """
    
    # do normalization based on the surface density distribution in side
    select_r_index = np.arange(len(select_range[select_range==True]))

    # renormalize the output
    for i_dust in range(int(rhodust.shape[-1])):
        # rho_t = rhodust[:, :, :, i_dust]
        rho_t = rhodust_all0 * mint.fracs[i_dust]
        integral_rho_origin = np.array([np.trapz(rho_t[i_r, :, 0], -zz[i_r, :, 0]) for i_r in select_r_index])
        
        rhodust_addwall_t = rhodust_addwall[:, :, :, i_dust]
        integral_rho_output = np.array([np.trapz(rhodust_addwall_t[i_r, :, 0], -zz[i_r, :, 0]) for i_r in select_r_index])
        
        for i_r in select_r_index:
            rhodust_addwall[i_r, :, :, i_dust] = rhodust_addwall[i_r, :, :, i_dust] * integral_rho_origin[i_r]/integral_rho_output[i_r]
    
    # **NOTE** (developing)
    # update the gas-to-dust mass ratio
    ratio_g2d_grid = rhogas / np.nansum(rhodust_addwall, axis=3) # rhodust_a.sum(3)
    ratio_g2d_grid[ratio_g2d_grid == np.inf] = np.nan
    # ratio_g2d_grid[np.isnan(ratio_g2d_grid)] = 0.0
    
    mint.ratio_g2d = ratio_g2d_grid
    
    print('the new ratio_g2d_grid at top inner disk is')
    print(ratio_g2d_grid[:10, 0, 0])
    
    print('the new ratio_g2d_grid at midplane inner disk is')
    print(ratio_g2d_grid[:10, -1, 0])
    # print(ratio_g2d_grid)
    
    #
    # output
    #
    # if verbose:
    print('the gas-to-dust mass ratio is also updated, and the numbers are saved to mint.ratio_g2d and a local file xxx_ratio_g2d_grid.dat')
    save_3dgrid_to_dat(ratio_g2d_grid, 'ratio_g2d_grid.new')
    
    if verbose:
        # test the total mass to see whether the mass is conservative to what it should be set upackage_position
        # mdust_out = (rhodust_addwall_nonorm * mint.vol).sum(0).sum(0).sum(0) * 2
        mdust_out_normed = (np.nansum(rhodust_addwall, axis=3) * mint.vol).sum(0).sum(0).sum(0) * 2
        print('dust mass before adding rim = %.2e [ms]'%(mdust_orig/const.ms))
        # print('dust mass after adding rim = %.2e [ms]'%(mdust_out/const.ms))
        print('dust mass after adding rim = %.2e [ms]'%(mdust_out_normed/const.ms))
        
        rhogas_t = np.nansum(rhodust_addwall, axis=3) * mint.ratio_g2d
        rhogas_t[np.isnan(rhogas_t)] = 1e-199
        rhogas_t[rhogas_t <= const.mu] = const.mu
        mgas_out = (rhogas_t * mint.vol).sum(0).sum(0).sum(0) * 2
        print('output gas mass after adding a rim = %.3e'%(mgas_out/const.ms))
        
        print('the current working dir for adding dust rim')
        print(os.getcwd())
    
    rhodust_a = rhodust_addwall.copy()
    rhodust_a = clean_rhodust_array(mint, rhodust_a)
    with open('dust_add_rim.new', 'w+') as f:
        f.write('1\n')  # Format number
        f.write('%d\n' % (nr * ntheta * nphi))  # Nr of cells
        f.write('%i\n' % (nrspec))  # Nr of dust species
        for i in range(nrspec):
            # data = mint.fracs[i] * rhonew.ravel(order='F')  # Create a 1-D view, fortran-style indexing
            data = rhodust_a[:,:,:,i].ravel(order='F')         # Create a 1-D view, fortran-style indexing
            data.tofile(f, sep='\n', format="%13.6e")
            f.write('\n')

    # if verbose:
    print('Done with adding dust rim, the new dust_add_rim.new are written')
    
    return 1
    

def chemistry_setup(mint, save_name, bool_same_rc_as_radmc3d=False):
    """
    setting up the initial grids and density for chemical network
    for now, we assume dust-gas are well coupled, so that
    the gas is simply get from dust no matter dust settled or not.
    
    Args:
        mint (object): The DiskMINT model object containing the parameters and data.
        save_name (str): The model name to be used for saving the output file.
        
        bool_same_rc_as_radmc3d (bool, optional): Whether to use the same rc grid as from the radmc3d files

    Returns:
        None.
        
        The files needed to run the chemical network is saved to current dir (where the code is running)
    """

    para = mint.parameters
    #  save_name = para.chemical_save_name
    
    nr_cyl = para.nr_cyl_LIME
    nz_cyl = para.nz_cyl_LIME

    """
    #############################
    # Setup the Coordinate Info #
    #############################
    """

    # nr_cyl = 100
    # nz_cyl = 200

    # save_name = chemical_save_name

    # print('Chemistry File To be Saved at: ', save_name)
    # print('Grid for LIME and Chemistry: (nr, nz) = ', nr_cyl, nz_cyl)

    ##################################
    #                                #
    # Read Data as the Initial Input #
    #                                #
    ##################################

    # dustdens = radmc3dPy.analyze.readData(dtemp=True, ddens=True, gtemp=False, binary=False)
    dustdens = modelgrid.readData(dtemp=True, ddens=True, gtemp=False) # , binary=False)
    nrspec = int(dustdens.rhodust.shape[-1])

    rhodust = dustdens.rhodust.copy()
    rhodust_all = np.nansum(rhodust, axis=3)
    
    ratio_g2d_grid_file = 'ratio_g2d_grid.dat'
    if os.path.exists(ratio_g2d_grid_file):
        print('read in the reference ratio_g2d_grid from %s'%(ratio_g2d_grid_file))
        ratio_g2d_grid_readin = load_3dgrid_from_dat(ratio_g2d_grid_file)
        ratio_g2d_grid = ratio_g2d_grid_readin
        ratio_g2d_grid[ratio_g2d_grid == np.inf] = np.nan

        mint.ratio_g2d = ratio_g2d_grid
    
    rhogas = mint.ratio_g2d * rhodust_all
    rhogas[np.isnan(rhogas)] = 0.0
    rhogas[rhogas <= const.mu] = const.mu
    
    ngas = rhogas / const.mu  # number density of molecule; 2.3
    nH = 2.0 * ngas.copy()  # number density of H nuclei; assume all H2
    # nH = rhogas / mu_atom # number density of H nuclei;
    nr = dustdens.grid.nx
    ntheta = dustdens.grid.ny
    nphi = dustdens.grid.nz
    rc = dustdens.grid.x
    thetac = dustdens.grid.y
    phic = dustdens.grid.z
    qq = np.meshgrid(rc, thetac, phic, indexing='ij')
    rr = qq[0]  # Spherical R 
    tt = qq[1]  # Angle with z axis ( = pi/2 for the z=0 axis)
    zr = np.pi / 2.e0 - qq[1]
    zz = np.sin(zr) * rr

    Tdust = dustdens.dusttemp.copy()
    Tdust[Tdust < 2.7] = 2.7
    
    Tgas = solve_gastemp(mint, dustdens, Tdust, 
                         bool_temp_decouple=mint.bool_temp_decouple)
    
    #
    # Compute the gas temperature
    #
    # Tgas = np.zeros([nr, ntheta, nphi])
    # #
    # # Compute the ndsd
    # #
    # # ndsd = nd * ng * a_ave**2
    # # ndsd = fracs_nb * a_ave**2 # number fractions for fracs
    # # ndsd_t = mint.ndsd 
    # # Compute ndust in a vectorized way
    # ndust = rhodust / (para.rhobulk * (4/3) * np.pi * (np.array(mint.a_ave) / 2) ** 3)

    # # Compute fracs_nb_eachsize only for the third dimension index 0
    # fracs_nb_eachsize = np.zeros_like(ndust)
    # ndsd_eachsize = np.zeros_like(ndust)

    # # Compute fraction and ndsd only for the third dimension index 0
    # sum_ndust = np.nansum(ndust[:, :, 0, :], axis=-1, keepdims=True)  # Sum over last axis
    # fracs_nb_eachsize[:, :, 0, :] = ndust[:, :, 0, :] / sum_ndust  # Broadcasting applies
    # ndsd_eachsize[:, :, 0, :] = fracs_nb_eachsize[:, :, 0, :] * np.array(mint.a_ave)**2  # Broadcasting applies
    
    # #
    # # solve the Tgas at different location
    # #
    # for i in range(0, nr):
    #     for j in range(0, ntheta):
    #         # ndsd = nd * ng * a_ave**2
    #         # ndsd = fracs_nb * a_ave**2 # number fractions for fracs
    #         # ndsd_t = mint.ndsd

    #         Tdmax_here = np.nanmax(Tdust[i, j, 0, :])
    #         ndsd_t = ndsd_eachsize[i, j, 0, :]
            
    #         Tgas[i, j, 0] = scipy.optimize.bisect(func_Tgas, 2.7, Tdmax_here, args=(nrspec, Tdust[i, j, 0, :], ndsd_t))

    Tgas_orig = Tgas.copy()

    """
    NOTE: Oct-18, 2022
    Whether add the Tgas_visc or not will not change the chemical network nor the final results
    """

    # 
    # NOTE: add the Tvisc heating here.
    #
    # Compute the Viscous Heating Tgas
    #

    #     Tgas_d = Tgas.copy()

    #     r_cyl_t = rr * np.cos(zr)
    #     fac1 = 3. * gg * para.mstar * para.mdotacc / (8. * np.pi * sigma_SB * r_cyl_t**3.)
    #     fac2 = 1. - np.sqrt(para.rstar/r_cyl_t)
    #     Tgas_visc = (fac1 * fac2)**(1./4.)

    #     Tgas = (Tgas_d**4. + Tgas_visc**4.)**(1./4.)

    """
    # NOTE: feed Tgas_d as Tdust to chemical code
    # NOTE: don't add Tgasvisc for chemical network
    # add . to all int
    """

    #####################################
    # Convert to Cylindrical Coordinate #
    #####################################

    #
    # grid info has been defined in the begginning
    #
    # rhogas_orig = rhogas.copy()
    rhogas_sph = rhogas.copy()
    rhogas_sph[rhogas_sph == 0] = 1e-199
    rhogas_sph[rhogas_sph <= const.mu] = const.mu
    
    # put most grid in the inner disk
    # but with denser grid in the outer part
    # first sample in log grid
    if bool_same_rc_as_radmc3d:
        rc_cyl_set = rc.copy()
        nr_cyl = len(rc_cyl_set)
        mint.parameters.nr_cyl_LIME = nr_cyl
        
    else:
        rc_cyl_set = None
    
    #
    # Convert both dust and gas into cylindrical coordinate
    #
    rc_cyl, zrc_cyl, qq_cyl, rhogas_cyl_old, Tgas_cyl_old = \
        sph2cyl(rc, thetac, phic, rhogas_sph, Tgas, nr_cyl=nr_cyl, nz_cyl=nz_cyl, rc_cyl_set=rc_cyl_set)
    
    _, _, _, rhodust_cyl_old = \
        sph2cyl_dust(rc, thetac, phic, rhodust, nr_cyl=nr_cyl, nz_cyl=nz_cyl, rc_cyl_set=rc_cyl_set)

    rhogas_cyl_old[np.isnan(rhogas_cyl_old)] = 1e-199
    rhogas_cyl_old[rhogas_cyl_old <= const.mu] = const.mu

    ngas_cyl_old = rhogas_cyl_old / const.mu
    
    rhodust_cyl_old_all = np.nansum(rhodust_cyl_old, axis=3)
        
    #
    # setup the ratio g2d for LIME grid
    #
    ratio_g2d_LIMEgrid = rhogas_cyl_old / (rhodust_cyl_old_all + 1e-199)
    ratio_g2d_LIMEgrid[rhogas_cyl_old <= 1e-198] = np.nanmax(ratio_g2d_LIMEgrid)
    
    # set up a grid of g2d in the shape of the model grid
    # 
    """
    **TODO** now this g2d grid for LIME (cylindrical) 
    doesn't work very well when the inner grid is small
    since the differences among the cylindrical and spherical
    coordinates. Will have a better function in the future.
    Now it is properly working for constant g2d, but 
    be cautious when you want to set a g2d as a function of r (cylindrical-r) or R (spherical-R).
    """
    # if hasattr(para, 'ratio_g2d') and para.ratio_g2d is not None:
    #     #  print('case 1')
    #     ratio_g2d = para.ratio_g2d
    # else:
    #     #  print('case 2')
    #     ratio_g2d = para.ratio_g2d_global
    
    # now the g2d is set as a 2d grid
    # ratio_g2d = mint.ratio_g2d

    #  print(para.dummy, 'this is the dummy')

    # ratio_g2d_set = copy.deepcopy(ratio_g2d)[:, :, 0]
    # ratio_g2d_LIMEgrid = place_ratio_g2d(ratio_g2d_set, shape=[nr_cyl, nz_cyl, 1], rc_cyl=rc_cyl, coord='cyl')

    # rhodust_cyl_old = np.zeros(list(rhogas_cyl_old.shape) + [nrspec])
    # for i_dust in range(nrspec):
    #     rhodust_cyl_old[:, :, :, i_dust] = rhogas_cyl_old / ratio_g2d_LIMEgrid * mint.fracs[i_dust]

    print('Grid in Cylindrical Coordinate: ', rhogas_cyl_old.shape, rc_cyl.shape, zrc_cyl.shape, rhodust_cyl_old.shape)

    #############################################
    # Compute tau in the cylindrical coordinate #
    #############################################
    #
    # Note here taux in spherical coordinate, tauz and dtauz in cylindrical coordinate
    #
    taux, tauz_up_cyl, tauz_dn_cyl, dtauz_cyl = calc_tau_new(dustdens, rc_cyl, zrc_cyl, rhodust_cyl_old, nrspec)

    #
    # If the taux in cylindrical coordinate is wanted, run the following block
    #
    taux_in = np.zeros([taux.shape[0], taux.shape[1], 1])
    taux_in[:, :, 0] = taux.copy()
    taux_in[taux_in == 0] = 1e-300

    ###########################################
    # also setup the Tdust as the Tdust(aave) #
    # where aave = np.sqrt(amax * amin)       #
    ###########################################

    a_ave_whole = np.sqrt(para.amin_all * para.amax_all)

    i_a_ave_whole = np.where(np.abs(mint.a_ave - a_ave_whole) == np.min(np.abs(mint.a_ave - a_ave_whole)))[0]
    n_iainrange = len(i_a_ave_whole)

    Tdust_aave_in_all = np.zeros(list(Tdust.shape[0:3]) + [n_iainrange])

    i_a_t = -1
    for i_a_ave_whole in np.where(np.abs(mint.a_ave - a_ave_whole) == np.min(np.abs(mint.a_ave - a_ave_whole)))[0]:
        i_a_t = i_a_t + 1
        Tdust_aave_in_all[:, :, :, i_a_t] = Tdust[:, :, :, i_a_ave_whole]

    """
    NOTE: now I am simply take an average of the Tdust
    that is close enough to the averaged a.
    """
    assert Tdust_aave_in_all.shape[3] == n_iainrange, 'not correct number of dust species are saved into the Tdust aave [%i] != [%i]' % \
        (Tdust_aave_in_all.shape[3], n_iainrange)
    Tdust_aave_in = Tdust_aave_in_all.sum(3) / n_iainrange

    rc_cyl, zrc_cyl, qq_cyl, taux_cyl_old, Tdust_aave_cyl_old = \
        sph2cyl(rc, thetac, phic, taux_in, Tdust_aave_in, nr_cyl=nr_cyl, nz_cyl=nz_cyl, rc_cyl_set=rc_cyl_set)

    ###############################################
    # Clean the output file for running chemistry #
    ###############################################

    # Replace both infinite and NaN values with 1e-299
    tauz_up_cyl  = np.where(np.isfinite(tauz_up_cyl), tauz_up_cyl, 1e-299)
    tauz_dn_cyl  = np.where(np.isfinite(tauz_dn_cyl), tauz_dn_cyl, 1e-299)
    dtauz_cyl    = np.where(np.isfinite(dtauz_cyl), dtauz_cyl, 1e-299)
    ngas_cyl_old = np.where(np.isfinite(ngas_cyl_old), ngas_cyl_old, 1e-299)
    
    # tauz_up_cyl[tauz_up_cyl == np.inf] = 1e-300
    # tauz_dn_cyl[tauz_dn_cyl == np.inf] = 1e-300
    # dtauz_cyl[dtauz_cyl == np.inf] = 1e-300

    # ngas_cyl_old[ngas_cyl_old == np.inf] = 1e-300
    # ngas_cyl_old[ngas_cyl_old == np.nan] = 1e-300

    #################################
    # Convert ngas (molecule) to nH #
    #################################

    # each H2 molecule has 2 H nuclei
    # and assume gas are all H2 molecule 
    # and H nuclei is mainly in H2 molecule
    # so that the total H nuclei number is simply H2 * 2.
    nH_cyl = ngas_cyl_old.copy() * 2.
    
    ##############################
    # calculate the dust moments #
    ##############################
    """ (new feature) 01/08/2025
    also save the dust moment here
    these will be used in the chemical network
    to get the correct grain surface chemistry
    """ 
    # convert the dust to cyl coordinate
    _, _, _, rhodust_cyl_old = \
        sph2cyl_dust(rc, thetac, phic, rhodust, nr_cyl=nr_cyl, nz_cyl=nz_cyl, rc_cyl_set=rc_cyl_set)
    
    # rhodust_cyl_old[np.isnan(rhodust_cyl_old)] = 1e-199
    # rhodust_cyl_old[rhodust_cyl_old == 0] = 1e-199
    
    rhodust_cyl_old = clean_rhodust_array(mint, rhodust_cyl_old)
    
    # rr_cyl  = qq_cyl[0]
    # zrr_cyl = qq_cyl[1]
    # zz_cyl  = zrr_cyl * rr_cyl
    
    # set up the empty arrays
    ndust = np.zeros(rhodust_cyl_old.shape) # number density of dust grains
    ndusta = np.zeros(rhodust_cyl_old.shape) # ndust * a
    ndusta2 = np.zeros(rhodust_cyl_old.shape) # ndust * a^2
    
    rhogr = para.rhobulk # the bulk density of the grain
    
    for i_dust in range(nrspec):
        adust = mint.a_ave[i_dust] # Dust size in cm 
        rhodust_cyl_a = rhodust_cyl_old[:, :, :, i_dust] # the mass density of the grain at size a
        
        massgr = rhogr * (4/3) * np.pi * adust**3 # the mass of the grain (note, the size a is the radii of the grain)
        ndust_cyl_a = rhodust_cyl_a / (massgr) # the number density of the grain at size a
        
        ndust[:, :, :, i_dust] = ndust_cyl_a
        ndusta[:, :, :, i_dust] = ndust_cyl_a * adust
        ndusta2[:, :, :, i_dust] = ndust_cyl_a * adust * adust
        
    # then calculate the sum of these dust moments
    ndust_all = np.nansum(ndust, axis=3)
    ndusta_all = np.nansum(ndusta, axis=3)
    ndusta2_all = np.nansum(ndusta2, axis=3)  
    
    #########################################
    # Save and Create the CO init Grid file #
    #########################################

    #
    # In Cyl Coordinate
    #
    row_numb = len(rc_cyl) * len(zrc_cyl)
    # col_numb = 11

    # X_H2   = 0.5             # Fraction of H2 molecule
    # ratio_ng2CO = 1.4e-4     # ratio from gas density to CO
    # conv_factor_12co = 1     # ratio from CO to 12CO
    # conv_factor_13co = 1/50  # ratio from CO to 13CO
    # conv_factor_c18o = 1/500 # ratio from CO to C18O
    # X_12CO = ratio_ng2CO
    # X_13CO = ratio_ng2CO * conv_factor_13co
    # X_C18O = ratio_ng2CO * conv_factor_c18o

    print('Chemistry File To be Saved at: ', save_name)
    print('Grid for LIME and Chemistry: (nr, nz) = ', nr_cyl, nz_cyl)
    
    save_file = 'COinitgrid-' + save_name + '-cyl.dat'
    print('The Saving Name: ', save_file)

    with open(save_file, 'w+') as f:
        f.write('ROW NUMBER=' + str(int(row_numb)) + ' (cgs units)\n')
        f.write('Nr=' + str(int(len(rc_cyl))) + '\n')
        f.write('Nz=' + str(int(len(zrc_cyl))) + '\n')
        f.write('r[cm],z/r,z[cm],n_H[cm**-3],Tgas[K],tauv_star,tauv_zup,tauv_zdn,dtauv_z,Tdust[K],ndust[cm**-3],ndusta[cm**-2],ndusta2[cm**-1],g2d_mass_ratio\n')
        for i_r in range(len(rc_cyl)):
            # for i_zr in range(len(zrc_cyl)):
            for i_zr in np.arange(len(zrc_cyl) - 1, -1, -1):
                # all in cgs units
                rincm_t = rc_cyl[i_r]
                zrr_t   = zrc_cyl[i_zr]
                zincm_t = zrc_cyl[i_zr] * rc_cyl[i_r]
                n_H_t   = nH_cyl[i_r, i_zr, 0]
                T_gas_t = Tgas_cyl_old[i_r, i_zr, 0]
                """
                NOTE: Tdust is only input to LIME and not to chemical network
                if you are not using LIME, then this column is not needed.
                """
                T_dust_t  = Tdust_aave_cyl_old[i_r, i_zr, 0]  # only input to LIME
                taux_t    = taux_cyl_old[i_r, i_zr, 0]
                tauz_up_t = tauz_up_cyl[i_r, i_zr]
                tauz_dn_t = tauz_dn_cyl[i_r, i_zr]
                dtauz_t   = dtauz_cyl[i_r, i_zr]
                """also save the dust moments here"""
                ndust_t   = ndust_all[i_r, i_zr, 0]
                ndusta_t  = ndusta_all[i_r, i_zr, 0]
                ndusta2_t = ndusta2_all[i_r, i_zr, 0]
                """now also save the gas-to-dust mass ratio"""
                ratiogtd_t = ratio_g2d_LIMEgrid[i_r, i_zr, 0]

                this_line = \
                    ' {:.3E},'.format(rincm_t) + \
                    ' {:.3E},'.format(zrr_t) + \
                    ' {:.3E},'.format(zincm_t) + \
                    ' {:.3E},'.format(n_H_t) + \
                    ' {:.3E},'.format(T_gas_t) + \
                    ' {:.3E},'.format(taux_t) + \
                    ' {:.3E},'.format(tauz_up_t) + \
                    ' {:.3E},'.format(tauz_dn_t) + \
                    ' {:.3E},'.format(dtauz_t) + \
                    ' {:.3E},'.format(T_dust_t) + \
                    ' {:.3E},'.format(ndust_t) + \
                    ' {:.3E},'.format(ndusta_t) + \
                    ' {:.3E},'.format(ndusta2_t) + \
                    ' {:.3E}'.format(ratiogtd_t) + \
                    '\n'
                
                # NOTE: last col does not contain ','

                # ' {:.3E},'.format(T_dust_t)  +\
                # ' {:.3E},'.format(X_H2)    +\
                # ' {:.3E},'.format(X_12CO)  +\
                # ' {:.3E},'.format(X_13CO)  +\
                # ' {:.3E}'.format(X_C18O)  +\

                f.write(this_line)
    
    """
    #
    # Sphr Coordinate is not used now
    #
    #
    # In Sphr Coordinate
    #
    alphac = np.pi/2 - thetac
    row_numb = len(rc) * len(alphac)

    # X_H2   = 0.5             # Fraction of H2 molecule
    # ratio_ng2CO = 1.4e-4     # ratio from gas density to CO
    # conv_factor_12co = 1     # ratio from CO to 12CO
    # conv_factor_13co = 1/50  # ratio from CO to 13CO
    # conv_factor_c18o = 1/500 # ratio from CO to C18O
    # X_12CO = ratio_ng2CO
    # X_13CO = ratio_ng2CO * conv_factor_13co
    # X_C18O = ratio_ng2CO * conv_factor_c18o

    with open('COinitgrid-' + save_name + '-sph.csv', 'w+') as f:
    #     f.write('ROW NUMBER:   ' + str(int(row_numb)) + ' (cgs units)\n')
        # f.write('R[cm],pi/2-theta,n_H[cm**-3],Tgas[K],Av_star,T_dust[K],X_H2,X_12CO,X_13CO,X_C18O\n')
        f.write('R[cm],pi/2-theta,n_H[cm**-3],Tgas[K],tauv_star\n')
        for i_r in range(len(rc)):
            # for i_zr in range(len(zrc_cyl)):
            for i_zr in np.arange(len(alphac)-1, -1, -1):
                # all in cgs units
                rincm_t   = rc[i_r]
                alpha_t   = alphac[i_zr]
                n_H_t     = ngas[i_r, i_zr, 0]
                T_gas_t   = Tgas[i_r, i_zr, 0]
                # T_dust_t  = Tdust[i_r, i_zr, 0, -1]
                taux_t    = taux[i_r, i_zr]

                this_line = \
                ' {:.3E},'.format(rincm_t)   +\
                ' {:.3E},'.format(alpha_t) +\
                ' {:.3E},'.format(n_H_t)     +\
                ' {:.3E},'.format(T_gas_t)   +\
                ' {:.3E}'.format(taux_t)    +\
                '\n'

                # ' {:.3E},'.format(T_dust_t)  +\
                # ' {:.3E},'.format(X_H2)    +\
                # ' {:.3E},'.format(X_12CO)  +\
                # ' {:.3E},'.format(X_13CO)  +\
                # ' {:.3E}'.format(X_C18O)  +\

                f.write(this_line)
    """
    return 1


def copy_files_with_extension(source_dir, destination_dir, extension, verbose=False):
    """
    The copy_files_with_extension function copies files with a specified extension from a source directory to a destination directory. It takes the following arguments:
    
    Args:
        source_dir (str): The path to the source directory where the files are located.
        destination_dir (str): The path to the destination directory where the files will be copied.
        extension (str): The extension of the files to be copied (e.g., ".txt", ".jpg").
        verbose (bool, optional): If set to True, the function will print the source and destination file paths for each copied file. Default is False.
    
    Returns:
        None.
        
        It copies files from source_dir to destination_dir
        
    Note:
        The function uses `shutil` to copy the files.
    """
    # Create the destination directory if it doesn't exist
    if not os.path.exists(destination_dir):
        os.makedirs(destination_dir)

    # Iterate over all files in the source directory
    for filename in os.listdir(source_dir):
        # Check if the file ends with the specified extension
        if filename.endswith(extension):
            source_file = os.path.join(source_dir, filename)
            destination_file = os.path.join(destination_dir, filename)
            shutil.copyfile(source_file, destination_file)
            if verbose:
                print(f"Copied: {source_file} -> {destination_file}")


def savemodel(mint, save_dir, model_dir='.'):
    """
    The savemodel function saves a model by copying specific files with given extensions from a source directory to a destination directory.
    
    Args:
        model (object): The DiskMINT model object containing the data and parameters to be saved.
        save_dir (str): The path to the destination directory where the model files will be saved.
        model_dir (str, optional): The path to the source directory where the model files are located. Default is '.' (current directory).
    """

    if not os.path.isdir(save_dir):
        os.makedirs(save_dir)

    if mint.bool_VHSE:
        extension_list = ['.inp', '.dat', '.out', '.py', '.new']
    else:  # no .new files if not using VHSE
        extension_list = ['.inp', '.dat', '.out', '.py']

    #  for ic in range(len(filein_list)):
    #  filein_name = modeldir + filein_list[ic]
    #  os.system("cp " + filein_name + " " + savedir)

    source_dir = model_dir
    destination_dir = save_dir
    print('saving the model from \n%s\nto %s' % (source_dir, destination_dir))
    # Iterate over the files and copy them to the destination directory
    for extension in extension_list:
        copy_files_with_extension(source_dir, destination_dir, extension)
        #  source_file = os.path.join(source_dir, file)
        #  destination_file = os.path.join(destination_dir)
        #  if os.path.isfile(source_file):
        #  print(source_file)
        #  print(destination_file)
        #  shutil.copy2(source_file, destination_file)

    # NOTE (Developing)
    # also save the set up parameter file
    para = mint.parameters
    para.write_parameters_to_csv(para.chemical_save_name+'_parameters', directory=save_dir, extension='.csv')
    
    if mint.bool_chemistry:
        # also save the chemistry part if the chemistry is runned
         
        src = os.path.join(model_dir, 'chemistry')
        dst = os.path.join(save_dir, 'chemistry')

        # Make sure the src (model_dir/chemistry) exists
        os.makedirs(src, exist_ok=True)  # <--- This line creates it if missing

        # Make sure the parent of dst exists (i.e., save_dir)
        os.makedirs(os.path.dirname(dst), exist_ok=True)

        # Remove existing save_dir/chemistry if it exists
        if os.path.exists(dst):
            shutil.rmtree(dst)

        # Copy folder
        shutil.copytree(src, dst)

        print(f"Copied {src} to {dst}")
    
    return 1


def write_chem_input(mint, chem_code_dir, G0Hab_set=0.0):
    """
    write the input files for the DiskMINT Fortran chemical network module
    
    Args:
        mint (object): The DiskMINT model object containing the data and parameters to be saved.
        chem_code_dir (str): where the chemistry code is put
        G0Hab_set (0.0 or other float > 0.0): the G0 at stellar surface in Hab units
    """

    para = mint.parameters

    #  chemical_save_dir = para.chemical_save_dir
    save_name = para.chemical_save_name
    print('the model named as \n%s' % (save_name))
    read_initgrid_file = 'COinitgrid-' + save_name + '-cyl.dat'

    # read in the grid info
    with open(read_initgrid_file, 'r') as f:
        lines = f.readlines()
        nr = int(lines[1][3:-1])
        nz = int(lines[2][3:-1])

    ###############################################################
    # Compute few quantities needed by the chemical model
    ###############################################################
    #
    # 1. Compute 0th, 1st and 2nd moment of the dust distribution
    #
    def dust_fun(nn):
        return (para.amax_all ** (nn - para.pla_dustsize) - para.amin_all ** (nn - para.pla_dustsize)) / \
               (nn - para.pla_dustsize)

    eta = 1.4 * const.mp / para.rhobulk / para.ratio_g2d_global
    nd = 3.0 * eta * dust_fun(1) / 4.0 / np.pi / dust_fun(4)
    nda = 3.0 * eta * dust_fun(2) / 4.0 / np.pi / dust_fun(4)
    nda2 = 3.0 * eta * dust_fun(3) / 4.0 / np.pi / dust_fun(4)
    
    #
    # 2. Compute G0 (UV radiation field) from spectrum (in Habing units)
    #
    #  f        = np.loadtxt('../../'+fmodel_filename)
    #  f        = np.loadtxt(mint.file_dir+para.fmodel_filename)
    #  wl       = f[:,0]
    #  Fnu      = f[:,1]*(3.08572e18/para.rstar)**2
    # (developing)
    # the one I was using 
    if G0Hab_set == 0.0:
        wl  = mint.lam
        Fnu = mint.Fnu * (3.08572e18 / para.rstar) ** 2
        nu  = const.clum / wl / 1.0e-4
        wll = const.hh * const.clum / 13.6 / const.eV / 1.0e-4
        wlu = const.hh * const.clum / 6.0 / const.eV / 1.0e-4
        Fnu[wl < wll] = 0.0
        Fnu[wl > wlu] = 0.0
        G0Hab = np.trapz(Fnu[::-1], nu[::-1]) / 1.6e-3
    else:
        G0Hab = G0Hab_set
        
    # 
    ###############################################################
    # Read structure computed with RADMC-3D and write it in a 
    # format compatible with the chemical model
    ###############################################################
    #
    f = np.loadtxt(read_initgrid_file, skiprows=4, delimiter=',')
    rdisk = np.reshape(f[:, 0], (nr, nz))
    zrdisk = np.reshape(f[:, 1], (nr, nz))
    zdisk = np.reshape(f[:, 2], (nr, nz))
    ndisk = np.reshape(f[:, 3], (nr, nz))  # nH
    Tgdisk = np.reshape(f[:, 4], (nr, nz))
    Avdisk = np.reshape(1.086 * f[:, 5], (nr, nz))
    Avup = np.reshape(1.086 * f[:, 6], (nr, nz))
    Avdn = np.reshape(1.086 * f[:, 7], (nr, nz))
    Tddisk = Tgdisk.copy()  # the area weighted Tdust = Tgas.

    """(developing) 01/08/2025
    in the new version, we need to input the <nd>, <nda>, and <nda2>
    for different locations (as function of r, z)
    they are calculated here and then be read in in the chemical network
    the file to be save/read is called dust.moment
    """ 
    nd_2d   = np.reshape(f[:, 10], (nr, nz))
    nda_2d  = np.reshape(f[:, 11], (nr, nz))
    nda2_2d = np.reshape(f[:, 12], (nr, nz))

    #
    # write out the initial required files for chemistry
    #
    #  directory = 'chemistry/data/struc/'
    directory = os.path.join(mint.file_dir, 'chemistry', 'data', 'struc')
    # Create the directory and any missing parent directories if they don't exist
    os.makedirs(directory, exist_ok=True)

    nrtmp = 0
    nztmp = np.zeros(nr, dtype=int)
    for i in range(0, nr):
        for j in range(nz - 1, -1, -1):
            if (ndisk[i, j] > 1.0):
                nztmp[i] += 1
        if (nztmp[i] > 0): nrtmp += 1
    # Write parameters files
    with open(os.path.join(directory, 'diskparameters'), 'w+') as f:
        f.write('%s \n' % para.rstar)
        f.write('%s \n' % G0Hab)
        f.write('%s \n' % nd)
        f.write('%s \n' % nda)
        f.write('%s \n' % nda2)
    
    # Write out the dust moment 
    with open(os.path.join(directory, 'dustmoment'), 'w+') as f:
        f.write('%s \n' % nrtmp)
        for i in range(0, nr):
            if (nztmp[i] > 0):
                f.write('%s \n' % nztmp[i])
                for j in range(0, nz):
                    if (ndisk[i, j] > 1.0):
                        f.write('%13.6e %13.6e %13.6e\n' % (
                        nd_2d[i, j], nda_2d[i, j], nda2_2d[i, j]))
        
    # Write the disk structure
    with open(os.path.join(directory, 'diskstructure'), 'w+') as f:
        f.write('%s \n' % nrtmp)
        for i in range(0, nr):
            if (nztmp[i] > 0):
                f.write('%s \n' % nztmp[i])
                for j in range(0, nz):
                    if (ndisk[i, j] > 1.0):
                        f.write('%13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e \n' % (
                        np.log10(rdisk[i, j] / const.au),
                        np.log10(zdisk[i, j] / const.au), 
                        np.log10(ndisk[i, j]), np.log10(Tgdisk[i, j]),
                        np.log10(Tddisk[i, j]),
                        Avdisk[i, j], Avup[i, j], Avdn[i, j]))
                        
    # Write the UV field
    with open(os.path.join(directory, 'uvfield'), 'w+') as f:
        for i in range(0, nr):
            for j in range(0, nz):
                if (ndisk[i, j] > 1.0):
                    f.write('%13.6e %13.6e %13.6e \n' % (Avdisk[i, j], Avup[i, j], Avdn[i, j]))

    """
    NOTE 
    This is used in the earlier version of the write_chem_input
    where we execute the chemistry code inside the diskmint package
    now we would execute the chemistry code inside each data file
    so we will need to copy other input files back instead.
    
    #
    # copy the strcuture to the chemical code directory
    #
    source_dir = directory
    destination_dir = os.path.join(chem_code_dir, 'data', 'struc')

    # Get a list of all files in the source directory
    #  files = os.listdir(source_dir)
    files = ['diskparameters', 'dustmoment', 'diskstructure', 'uvfield']

    # Iterate over the files and copy them to the destination directory
    for file in files:
        source_file = os.path.join(source_dir, file)
        destination_file = os.path.join(destination_dir, file)
        if os.path.isfile(source_file):
            shutil.copy2(source_file, destination_file)
    """
    
    #
    # copy the chem input from the chemical code dir
    #
    source_dir = os.path.join(chem_code_dir, 'data', 'chem')
    destination_dir = os.path.join(mint.file_dir, 'chemistry', 'data', 'chem')

    print("Checking paths:")
    print("SOURCE:     ", source_dir)
    print("DESTINATION:", destination_dir)

    # Make sure the source exists
    if not os.path.exists(source_dir):
        raise FileNotFoundError(f"Source folder does not exist: {source_dir}")

    # Create destination's parent directory
    os.makedirs(os.path.dirname(destination_dir), exist_ok=True)

    # If destination exists, delete it (clean copy)
    if os.path.exists(destination_dir):
        print("Removing existing destination folder...")
        shutil.rmtree(destination_dir)

    # Copy folder
    print("Copying folder...")
    shutil.copytree(source_dir, destination_dir)
    print(f"Copied chemistry input folder from:\n  {source_dir}\nto:\n  {destination_dir}")

    # Confirm contents
    print("Contents of copied folder:")
    for root, dirs, files in os.walk(destination_dir):
        for name in files:
            print(os.path.join(root, name))
    
    print(f"Copied chemistry input folder from:\n  {source_dir}\nto:\n  {destination_dir}")
    
    return 1


def write_chem_output(mint, chem_code_exe_dir, chem_network_name='reducedRGH22', bool_save_allab=False):
    """
    write the chemical output (Fortran style) in the Python style back to the model directory
    
    Args:
        mint (object): The DiskMINT model object containing the data and parameters to be saved.
        chem_code_exe_dir (str): where the chemistry code is executed
        chem_network_name (str): the name of the chemical network (just for recording)
        bool_save_allab (bool): if true, then save all of the abundance files from the chemistry output to the main data folder
    """

    para = mint.parameters

    write_dir = mint.file_dir
    save_name = para.chemical_save_name

    print('saving the model of %s to\n%s' % (save_name, write_dir))

    # read in the grid info
    file_COinitgrid = os.path.join(write_dir, 'COinitgrid-' + save_name + '-cyl.dat')
    with open(file_COinitgrid, 'r') as f:
        lines = f.readlines()
        nr = int(lines[1][3:-1])
        nz = int(lines[2][3:-1])

    print('load initgrid ' + file_COinitgrid)

    # setup the grid for output
    f = np.loadtxt(file_COinitgrid, skiprows=4, delimiter=',')
    rdisk = f[:, 0]
    zdisk = f[:, 2]
    ndisk = f[:, 3]
    abh2 = np.zeros(nr * nz)
    abh2[:] = -99
    abc18o = np.zeros(nr * nz)
    abc18o[:] = -99

    print(ndisk.shape)
    print(abh2.shape)

    print('Reading/saving the files from %s'%(chem_code_exe_dir))
    
    outab_h2 = os.path.join(chem_code_exe_dir, 'out', 'ab', 'H2.txt')
    #  f = np.loadtxt('out/ab/H2.txt',skiprows=1)
    f = np.loadtxt(outab_h2, skiprows=1)
    abh2tmp = f[:, 2]
    outab_c18o = os.path.join(chem_code_exe_dir, 'out', 'ab', 'CO@.txt')
    #  f = np.loadtxt('out/ab/CO@.txt',skiprows=1)
    f = np.loadtxt(outab_c18o, skiprows=1)
    abc18otmp = f[:, 2]

    print(abh2tmp.shape)
    print(abh2[ndisk > 1.0].shape)

    abh2[ndisk > 1.0] = abh2tmp
    abc18o[ndisk > 1.0] = abc18otmp

    file_COendgrid = os.path.join(write_dir, 'COendgrid-' + chem_network_name + '-' + save_name + '-cyl.chem')
    with open(file_COendgrid, 'w+') as f:
        k = 0
        for i in range(nr):
            for j in range(nz):
                f.write('%13.5e %13.5e %13.5e %13.5e\n' % (rdisk[k], zdisk[k], abh2[k], abc18o[k]))
                k += 1

    print('Done, saving file to %s' % (file_COendgrid))
     
    if bool_save_allab:
        # save all the ab files from the chemistry output to the main data/chemistry/ab
        #
        # copy the strcuture to the chemical code directory
        #
        source_dir = os.path.join(chem_code_exe_dir, 'out', 'ab')
        destination_dir = os.path.join(write_dir, 'chemistry', 'out', 'ab')
        
        # make (sure exist) of the dir requested is there; if not, make them
        os.makedirs(destination_dir, exist_ok=True)
        print(f"Directory ready: {destination_dir}")
        
        # Get a list of all files in the source directory
        files = os.listdir(source_dir)

        # Iterate over the files and copy them to the destination directory
        for file in files:
            source_file = os.path.join(source_dir, file)
            destination_file = os.path.join(destination_dir, file)
            if os.path.isfile(source_file):
                shutil.copy2(source_file, destination_file)

    return 1


def runchemistry(mint, chem_code_dir, G0Hab_set=0.0, bool_save_allab=False, verbose=True):
    """
    Run the chemical network of the DiskMINT model.

    Args:
        mint (object): The DiskMINT model object containing the data and parameters to be saved.
        chem_code_dir (str): where the chemistry code is put
        G0Hab_set (float): The G0 at the stellar surface in Habing units (optional).
        bool_save_allab (bool): Whether to save all abundance outputs to the chemistry data directory.
        verbose (bool): If true, printout some texts for debuggin and monitoring
    """

    # Absolute paths and working directories
    target_dir = radmc_data_dir = mint.file_dir
    os.chdir(target_dir)
    
    # Print initial directory info
    current_dir = os.getcwd() # os.path.dirname(os.path.abspath(__file__))

    if verbose:
        print('changed to radmc3d data directory:\n' + current_dir)
    
    # data_dir = mint.file_dir # os.path.join(target_dir, 'data')
    # chem_code_dir = os.path.join(package_position, 'chemistry', 'reducedRGH22')
    """NOTE: now all the chemistry will be running inside the local folder"""
    chem_code_exe_dir = os.path.join(mint.file_dir, 'chemistry')
    
    fortran_executable_disk_main = os.path.abspath(os.path.join(chem_code_dir, '..', 'bin', 'disk_main'))
    fortran_executable_disk_extract = os.path.abspath(os.path.join(chem_code_dir, '..', 'bin', 'disk_extract'))

    # Write input files for the chemical network
    write_chem_input(mint, chem_code_dir, G0Hab_set=G0Hab_set)

    # Print initial directory info again to make sure
    current_dir = os.getcwd()
    if verbose:
        print('Original script directory:\n' + current_dir)

    # Change to chemistry working directory
    if verbose:
        print('Changing working directory to the chemical network dir:')
        print('This dir should be X/chemistry/')
        print(chem_code_exe_dir)
    os.chdir(chem_code_exe_dir)

    current_dir = os.getcwd()
    if verbose:
        print('Current working directory (after change):\n' + current_dir)

    # Timing the run
    if verbose:
        time0 = time.time()
        print('--------------------')
        print(f'Solving Chemistry')
        print(f'START at {time0:.2f} [s]')
        print('Wall time:', time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())) 
        print(f'Elapsed = {time0 / 3600.0:.2f} HOURS')
        print('--------------------')

    # Run Fortran chemistry code
    subprocess.run([fortran_executable_disk_main])
    subprocess.run([fortran_executable_disk_extract])

    # Output timing info
    if verbose:
        time_elapsed = time.time() - time0
        print('--------------------')
        print(f'END after {time_elapsed:.2f} [s]')
        print('Wall time:', time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())) 
        print(f'Cost = {time_elapsed / 3600.0:.2f} HOURS')
        print('--------------------')

    # Write chemical network output
    write_chem_output(mint, chem_code_exe_dir, chem_network_name=os.path.basename(chem_code_dir), bool_save_allab=bool_save_allab)

    # Return to original working directory
    os.chdir(radmc_data_dir)
    if verbose:
        print('Returned to original working directory:\n' + os.getcwd())

        print('Done! Chemical network finished. Use LIME or RADMC-3D for radiative transfer.')

    return 1


def runchemistry_prev(mint, chem_code_dir, G0Hab_set=0.0, bool_save_allab=False):
    """
    run the chemical network of DiskMINT model
    This is the previous version that the code is executed in the diskmint package location
    
    Args:
        mint (object): The DiskMINT model object containing the data and parameters to be saved.
        chem_code_dir (str): where the chemistry code is put
        G0Hab_set (None or float): the G0 at stellar surface in Hab units
        bool_save_allab (bool): if true, save all the ab output to the data/chem dir
    """

    # Get the absolute path to the directory
    current_dir = os.path.dirname(os.path.abspath(__file__))
    print('current working directory:\n' + current_dir)
    code_start_dir = current_dir

    #
    # Write the input files to Fortran format 
    # feeding the chemical network
    #
    write_chem_input(mint, chem_code_dir, G0Hab_set=G0Hab_set)

    # Change the current working directory to the path containing the radmc3d model files
    print('changing working dir to the chemical network dir: ')
    print(chem_code_dir)
    target_dir = chem_code_dir
    os.chdir(target_dir)

    current_dir = os.path.dirname(os.path.abspath(__file__))
    print('current working directory (after changed):\n' + current_dir)

    #
    # run the chemical network 
    #
    # Specify the path to the Fortran executable
    fortran_executable_disk_main = os.path.join('..', 'bin', 'disk_main')
    fortran_executable_disk_extract = os.path.join('..', 'bin', 'disk_extract')

    # Run the Fortran executable
    subprocess.run([fortran_executable_disk_main])
    subprocess.run([fortran_executable_disk_extract])

    #  os.system('../bin/disk_main')
    #  os.system('../bin/disk_extract')

    # 
    # write the output
    # NOTE: Warnning - in this version here is a bug that the chem_code_name is not specified
    # and will always output 'reducedRGH22'
    write_chem_output(mint, chem_code_dir, bool_save_allab=bool_save_allab)

    #
    # go back to the model folder
    #
    os.chdir(code_start_dir)

    current_dir = os.path.dirname(os.path.abspath(__file__))
    print('current working directory (after changed):\n' + current_dir)

    print('Done! chemical network finished. Please use LIME or RADMC3D or other code to do the following line radiative transfer')

    return 1


# def runmodel(mint, test_alliteration=False):
#     """
#     Run the DiskMINT model
    
#     Args:
#         mint (object): The DiskMINT model object containing the data and parameters to be saved.
#         test_alliteration (bool): If True, force to run the all the iterations set up for VHSE to see whether it is truly converged at the last few iterations.
#     """

#     proceed_code = True
    
#     para = mint.parameters
    
#     if not os.path.isdir(mint.file_dir):
#         os.mkdir(mint.file_dir)

#     # Get the absolute path to the directory
#     current_dir = os.path.dirname(os.path.abspath(__file__))
#     print('current working directory:\n' + current_dir)
#     code_start_dir = current_dir

#     # Change the current working directory to the path containing the radmc3d model files
#     target_dir = mint.file_dir
#     os.chdir(target_dir)

#     current_dir = os.path.dirname(os.path.abspath(__file__))
#     print('current working directory (after changed):\n' + current_dir)

#     if mint.bool_MakeDustKappa:
#         print('the DiskMINT module does not support making the opacity files here')
#         print('Please make your own Dust Kappa and feed into the Model')
#         print('Here the setup amin and amax are used to make the neccessary information')
#         #  os.system("rm dustkappa_*.inp")
#         #  os.system("python makedustkappa.py")
#         mint.setup_dust_info(para)
#         #  elif len(para.ndsd)
#     else:
#         print('reading pre set a_ave.inp, fracs.inp, fracs_nb.inp, and ndsd.inp')
#         #  then the dust info needs to in read from .inp files
#         mint.readin_dust_info(file_dir=mint.file_dir)
#         print('pre set a_ave.inp, fracs.inp, fracs_nb.inp, and ndsd.inp are read')

#     #
#     # Clean old files
#     #
#     # os.system("rm dust_density_*.dat")
#     # os.system("rm dust_temperature_*.dat")
#     # os.system("rm gas_temperature_*.new")
#     file_list = ['dust_density.inp', 'dust_temperature.dat', 'gas_temperature.inp', 'gas_temperature.new', 'dust_density_*.dat', 'dust_temperature_*.dat', 'gas_temperature_*.new', 'gas_temperature_*.inp', 'dust_density_*.inp', 'dust_density_*.new', 'dust_density_settled.new', 'ratio_g2d_grid_*.dat', 'ratio_g2d_grid_*.new','ratio_g2d_grid.*', 'rhodust_all_probset.inp', 'rhogas_probset.inp', 'rhodust_probset.inp', 'dust_add_rim.new', 'taux.dat', 'tauz_dn.dat', 'tauz_up.dat', 'spectrum_*.out', 'spectrum.out']
#     for file_name_i in file_list:
#         os.system("rm %s"%(file_name_i))

#     #
#     # Call problem_setup
#     #
#     command_mctherm = 'radmc3d mctherm'
    
#     # start over or read in the previous files
#     if mint.bool_startover:
#         #
#         # The Module Below relies on reading in the new dust opacity, 
#         # so import it after the dust opacity setup and old files removed.
#         #
#         #  import wrapper_scripts.wrapper_disk_density as dd

#         problem_setup(mint)
#         print("Done with problem_setup setup")

#         nr_redundant_itr = 3
#         i_redundant_iter = 1
        
#         # comp dust temperature
#         os.system(command_mctherm)
        
#         fname = 'dust_temperature.dat'
#         key_out = 0
#         while not key_out:
#             if os.path.isfile(fname):
#                 print('successfully obtained the dust_temperature.dat')
#                 key_out = 1
            
#             elif i_redundant_iter <= nr_redundant_itr:
#                 i_redundant_iter += 1
#                 print("run the initial radmc3d mctherm for #%i time"%(i_redundant_iter))
                
#                 # os.system(command_mctherm)
#                 get_new_dust_temperature(mint, command=command_mctherm, bool_clean_dusttemp_midplane=mint.bool_clean_dusttemp_midplane)
                
#             else:
#                 print('Something wrong with the initial iteration, check the code before proceeding')
#                 key_out = 1
#                 proceed_code = False
    
#     # read in the previous files        
#     else:
#         print('Using previous setups as the initial Gaussian distribution')

#         problem_setup(mint)
#         print("Done with problem_setup setup")
#         print('Copy the files needed from previous runs')

#         # file_dir_readin = '../result_radmcfiles/' + prev_setup_name + '/'
#         file_dir_readin = mint.data_dir_readin + mint.prev_setup_name + '/'

#         print(file_dir_readin)

#         if mint.read_in_type == 'Gaussian':
#             filein_list = ['stars.inp', 'dustopac.inp', 'fracs.inp', 'aave.inp', 'fracs_numb.inp', 'amr_grid.inp',
#                            'wavelength_micron.inp', 'dust_density.inp', 'dust_temperature.dat', 'spectrum_init_GS.out']

#         elif mint.read_in_type == 'VHSE':
#             filein_list = ['stars.inp', 'dustopac.inp', 'fracs.inp', 'aave.inp', 'fracs_numb.inp', 'amr_grid.inp',
#                            'wavelength_micron.inp', 'dust_density_mdust_VHSE_1.dat',
#                            'dust_temperature_mdust_VHSE_1.dat', 'spectrum_init_GS.out']

#         elif mint.read_in_type == 'VHSEend':
#             filein_list = ['stars.inp', 'dustopac.inp', 'fracs.inp', 'aave.inp', 'fracs_numb.inp', 'amr_grid.inp',
#                            'wavelength_micron.inp', 'dust_density.inp', 'dust_temperature.dat', 'spectrum_init_GS.out']

#         file_list = ['stars.inp', 'dustopac.inp', 'fracs.inp', 'aave.inp', 'fracs_numb.inp', 'amr_grid.inp',
#                      'wavelength_micron.inp', 'dust_density.inp', 'dust_temperature.dat', 'spectrum.out']

#         # reset the rhogas_probset.inp
#         #         rhogas0_prev        = np.loadtxt(file_dir_readin+'rhogas_probset.inp')
#         #         ratio_g2d_prev_grid = place_ratio_g2d(mint.ratio_g2d_prev)
#         #         rhogas0             = rhogas0_prev/mint.ratio_g2d_prev*mint.ratio_g2d

#         #         np.savetxt('rhogas_probset.inp', rhogas0)
#         #         mint.rhogas0 = rhogas0

#         for ic, file_name_i in enumerate(file_list):
#             filename = file_dir_readin + filein_list[ic]
#             os.system("cp " + filename + " " + file_name_i)

#         #
#         # The Module Below relies on reading in the new dust opacity, 
#         # so import it after the dust opacity setup and old files removed.
#         #
#         #  import wrapper_scripts.wrapper_disk_density as dd

#         os.system('rm dustkappa_*.inp')
#         os.system('cp ' + file_dir_readin + 'dustkappa_*.inp ./')

#     command_sed = "radmc3d sed incl %.1f phi %.1f" % (para.incl, para.phi)  # setthreads 10"
#     # Compute the sed for the Gaussian Model
#     if mint.bool_SED and mint.bool_startover and proceed_code:
#         if not mint.bool_VHSE:
#             os.system(command_sed)
#             filename = "spectrum_init_GS.out"
#             os.system("cp spectrum.out " + filename)

#     # Save the initial Gaussian disk before VHSE calculation    
#     if mint.bool_chemistry and proceed_code:
#         chemistry_setup(mint, save_name='GSinit_' + para.chemical_save_name, bool_same_rc_as_radmc3d=mint.bool_same_rc_as_radmc3d)

#     # set up the test values to see when to stop the vhse iterations
#     difference_test_max_prev_itr = np.inf
#     difference_test_mean_prev_itr = np.inf
    
#     # start the main loop for solving VHSE
#     if mint.bool_VHSE and proceed_code:
#         # results above needs to be read in
#         # when import the following package
#         # so import here

#         nloops = mint.n_vhse_loop  # set for testing now
#         nr_redundant_itr = int(nloops / 1.0) # set the maximum roduntent iteration
#         key_out = 0  # see whether to stop iterations

#         # now we start from i_iter 0 for the initial Gaussian distribution
#         i_iter = -1 
#         i_redundant_iter = 0
#         # for i_iter in range(1, nloops + 1):
#         # while i_iter < nloops and key_out == 0:
#         while key_out == 0:
            
#             i_iter = int(i_iter + 1)
#             print('start iteration #%i'%(i_iter))
            
#             """
#             previous version,
#             solve VHSE and dust settling both based on 
#             the previous density and thermal structure
#             """
            
#             # #
#             # # solve vhse
#             # #
#             # if i_iter >= 1:
#             #     bool_dust_settling_in_iteration_t = mint.bool_dust_settling
#             # else:
#             #     # we do not apply dust settling in the first iteration
#             #     bool_dust_settling_in_iteration_t = False

#             # #
#             # # compute the dust density structure
#             # #
#             # """NOTE (developing) now we integrate the dust settling code inside the vhse loop"""
#             # get_new_dust_density(mint, bool_dust_settling_in_iteration=bool_dust_settling_in_iteration_t)
            
#             """
#             01/08/2025
#             now in each iteration after i_iter 1,
#                 first, we include dust settling,
#                 then, we run radmc3d mctherm to get new Tdust (and Tgas),
#                 fiannly, we solve vhse on the settled dust and new Tgas
#             then go back to first step above
#             """
#             if i_iter >= 1:
                
#                 if mint.bool_dust_settling:
#                     #
#                     # solve the dust settling
#                     #
#                     set_dust_settling(mint, verbose=True)
#                     print('dust and gas are de-coupled, rhodust is updated according to the settling solver at itr #%i'%(i_iter))
                    
#                     # move the new files as the input
#                     os.system('mv dust_density_settled.new dust_density.inp')
#                     os.system("mv ratio_g2d_grid.new ratio_g2d_grid.dat")
                
#                 # then we need to do the thermal calculation again to get the new dust Temperature
#                 # os.system("radmc3d mctherm")
#                 get_new_dust_temperature(mint, command=command_mctherm, bool_clean_dusttemp_midplane=mint.bool_clean_dusttemp_midplane)
            
#             # else: continue
#             # we do not apply dust settling in the first iteration
            
#             """
#             # read in the dust density 1 (previous density from last iteration)
#             # this is for comparison and see whether the density converges
#             """
#             # dustdens1 = radmc3dPy.analyze.readData(dtemp=False, ddens=True, gtemp=True)
#             dustdens1 = modelgrid.readData(dtemp=False, ddens=True, gtemp=True)
#             rhodust1  = dustdens1.rhodust.copy()
#             rhogas1   = np.nansum(rhodust1, axis=3) * mint.ratio_g2d
            
#             if i_iter >= 1:
#                 # in the first iteration before solving VHSE
#                 # Tgas is not available.
#                 Tgas1_0   = dustdens1.gastemp.copy()
#                 Tgas1_0   = Tgas1_0[:, :, :, 0]
            
#             #
#             # solve VHSE
#             # in the case when dust settling is turned on (dust-gas de-coupled)
#             # this function is getting the new gas density by solving VHSE
#             # rather than really get the dust density
#             get_new_dust_density(mint, 
#                                  bool_dust_settling_in_iteration=mint.bool_dust_settling) # ,
#                                 #  bool_dust_fragmentation_set=mint.bool_dust_fragmentation,
#                                 #  bool_dust_radial_drifting_set=mint.bool_dust_radial_drifting)
            
#             #
#             # Save old files and update dust density file
#             #
#             # the dust density for the previous iteration
#             filename0 = "dust_density_mdust_VHSE_" + str(i_iter) + ".dat"
#             os.system("cp dust_density.inp " + filename0)
#             # the dust temperature for the previous iteration
#             filename1 = "dust_temperature_mdust_VHSE_" + str(i_iter) + ".dat"
#             os.system("cp dust_temperature.dat " + filename1)
#             # the gas temperature for the previous iteration
#             filename2 = "gas_temperature_mdust_VHSE_" + str(i_iter) + ".inp"
#             os.system("cp gas_temperature.inp " + filename2)
            
#             print('the input files for this iteration (= the files for previous iteration) are saved as\n%s\n%s\n%s'%(filename0, filename1, filename2))
            
#             # the ratio g2d grid
#             if os.path.exists('ratio_g2d_grid.dat'):
#                 filename3 = "ratio_g2d_grid_VHSE_" + str(i_iter) + ".dat"
#                 os.system("cp ratio_g2d_grid.dat " + filename3)
                
#                 print('%s is also saved'%(filename3))
                
#             if mint.bool_clean_dusttemp_midplane:
#                 # the original dust temperature for the previous iteration
#                 filename4 = "dust_temperature_originalmidplane_mdust_VHSE_" + str(i_iter) + ".dat"
#                 os.system("cp dust_temperature_originalmidplane.dat " + filename4)

#                 print('%s is also saved'%(filename4))

#             # 
#             # # the NEW ratio g2d
#             # filename4 = "ratio_g2d_grid_VHSE_" + str(i_iter) + ".new"
#             # os.system("cp ratio_g2d_grid.new " + filename4)
            
#             # put the newly calculated dust density and g2d grid to the inp files.
#             os.system("mv dust_density.new dust_density.inp")
#             os.system("mv ratio_g2d_grid.new ratio_g2d_grid.dat")
            
#             #
#             # calculate dust temperature
#             #
#             # os.system("radmc3d mctherm")
            
#             # """put in the dust settling during the VHSE iterations"""
#             # if mint.bool_dust_settling & proceed_code:
#             #     #
#             #     # solve the dust settling
#             #     #
#             #     set_dust_settling(mint, verbose=True)
                
#             #     # the dust density for the previous iteration
#             #     filename0 = "dust_density_mdust_VHSE_" + str(i_iter+1) + "_before_dustset.dat"
#             #     os.system("cp dust_density.inp " + filename0)
#             #     # the dust temperature for the previous iteration
#             #     filename1 = "dust_temperature_mdust_VHSE_" + str(i_iter+1) + "_before_dustset.dat"
#             #     os.system("cp dust_temperature.dat " + filename1)
#             #     # the gas temperature for the previous iteration
#             #     filename2 = "gas_temperature_mdust_VHSE_" + str(i_iter+1) + "_before_dustset.inp"
#             #     os.system("cp gas_temperature.inp " + filename2)
#             #     # the ratio g2d grid
#             #     if os.path.exists('ratio_g2d_grid.dat'):
#             #         filename3 = "ratio_g2d_grid_VHSE_" + str(i_iter+1) + "_before_dustset.dat"
#             #         os.system("cp ratio_g2d_grid.dat " + filename3)

#             #     print('the files for this iteration are saved as\n%s\n%s\n%s'%(filename0, filename1, filename2))
#             #     if os.path.exists('ratio_g2d_grid.dat'):
#             #             print('%s is also saved'%(filename3))
                
#             #     # move the new files as the input
#             #     os.system('mv dust_density_settled.new dust_density.inp')
#             #     os.system("mv ratio_g2d_grid.new ratio_g2d_grid.dat")
                
#             #     # then we need to do the thermal calculation again to get the new dust Temperature
#             #     os.system("radmc3d mctherm")
#             #     # then if we need to get the new gas temperature based on this new dust thermal structure
#             #     # dustdens = radmc3dPy.analyze.readData(dtemp=True, ddens=True, binary=False)
#             #     # Tdust = dustdens.dusttemp.copy()
#             #     # Tdust[Tdust < 2.7] = 2.7
#             #     # solve_gastemp(mint, dustdens, Tdust, bool_output_gastemp=True)
#             #     # os.system('mv gas_temperature_direct_output.inp gas_temperature.inp')
                
#             # """add the dust inner rim"""
#             # if mint.bool_dust_inner_rim:
#             #     """(devloping)"""
#             #     add_dust_rim(mint, verbose=True)
#             #     os.system('cp dust_add_rim.new dust_density.inp')
#             #     os.system("mv ratio_g2d_grid.new ratio_g2d_grid.dat")
            
#             """
#             # read in the dust density 2 (new density this iteration)
#             """
#             # dustdens2 = radmc3dPy.analyze.readData(dtemp=False, ddens=True, gtemp=True)
#             dustdens2 = modelgrid.readData(dtemp=False, ddens=True, gtemp=True)
#             rhodust2  = dustdens2.rhodust.copy()
#             rhogas2   = np.nansum(rhodust2, axis=3) * mint.ratio_g2d
            
#             if i_iter >= 1:
#                 # in the first iteration before solving VHSE
#                 # Tgas is not available.
#                 Tgas2_0   = dustdens2.gastemp.copy()
#                 Tgas2_0   = Tgas2_0[:, :, :, 0]
            
#             # rc      = dustdens1.grid.x.copy()
#             # thetac  = dustdens1.grid.y.copy()
#             # rhodust_all_1 = dustdens1.rhodust.sum(3)
#             # rhodust_all_2 = dustdens2.rhodust.sum(3)
#             # difference_test = np.abs((rhodust_all_2 - rhodust_all_1) / rhodust_all_1)
#             # difference_test[rhodust_all_2 <= 1e-20] = 0
            
#             if i_iter >= 1:
#                 # in the first iteration before solving VHSE
#                 # Tgas is not available.
#                 difference_test = np.abs((Tgas2_0 - Tgas1_0)/Tgas1_0)
#                 difference_test[rhogas2 <= 1e-20] = 0
#                 difference_test[np.isnan(difference_test)] = 0
#                 difference_test_max_t = np.nanmax(difference_test)
#                 difference_test_mean_t = np.nanmean(difference_test)
#                 print('-------------')
#                 print('test the Tgas differences of rho between the two iterations')
#                 print('delta_T/T with mean %.3e and max %.3e' % (difference_test_mean_t, difference_test_max_t))

#             difference_test = np.abs((rhogas2 - rhogas1)/rhogas1)
#             difference_test[rhogas2 <= 1e-20] = 0
#             difference_test[np.isnan(difference_test)] = 0
#             difference_test_max_t = np.nanmax(difference_test)
#             difference_test_mean_t = np.nanmean(difference_test)
#             print('-------------')
#             print('test the rhogas differences of rho between the two iterations')
#             print('delta_rho/rho with mean %.3e and max %.3e' % (difference_test_mean_t, difference_test_max_t))
#             print('-------------')
            
#             if not test_alliteration:
#                 if np.nanmean(difference_test) < 1e-20:
#                     print('the new density seem to be identical to the previous iteration')
                    
#                     if i_redundant_iter < nr_redundant_itr:
#                         print('something wrong in this iteration, repeat this iteration #%i'%(i_iter))
                        
#                         i_iter = int(i_iter - 1)
#                         i_redundant_iter = int(i_redundant_iter + 1)
                        
#                         print('this will be the the #%i redundant iteration'%(i_redundant_iter))
                    
#                     else:
#                         print('a lot iterations have issues, meet the stop criteria, please check the code before proceeding')
#                         key_out = 1
#                         # proceed_code = False
                
#                 else:
                    
#                     different_criteria = 0.05 # 0.1 and 0.05 tried in the past
#                     if np.nanmax(difference_test) < different_criteria:
#                         print('the delta_rho/rho is less than criteria of %.2f'%(different_criteria))
#                         print('meet the criteria at itr %i' % (i_iter))
#                         proceed_code = True
#                         key_out = 1 
                        
#                     elif (difference_test_mean_t >= difference_test_mean_prev_itr) & (difference_test_max_t >= difference_test_max_prev_itr):
#                         print('the differences starts to increase from this iteration')
#                         print('the delta_rho/rho in previous iterations are with\nmean %.3e and max %.3e' % (difference_test_mean_prev_itr, difference_test_max_prev_itr))
#                         print('the current set up might not be able to decrease this difference and converge into the set up criteria')
#                         print('please consider increase nr and nphot')
#                         print('code still proceed; be careful with the results.')
#                         proceed_code = True
#                         key_out = 1
                        
#                     else:
#                         difference_test_max_prev_itr  = difference_test_max_t
#                         difference_test_mean_prev_itr = difference_test_mean_t
            
#             print('-------------')
#             print("i_iter %i done" % (i_iter))
#             print('-------------')
            
#             if i_iter >= nloops:
#                 print('Went through (%i) iterations, larger than the requested number (%i); loop out'%(i_iter, nloops))
#                 key_out = 1
            
#             if key_out == 1:
#                 print('this is the last iteration in VHSE')
#                 print('-------------')
                
#                 # # record the last ratio_g2d_grid as dat instead of new
#                 # filename = "ratio_g2d_grid.dat"
#                 # os.system("cp ratio_g2d_grid.new " + filename)
#                 if not proceed_code:
#                     if i_iter > int(nloops / 2.0):
#                         proceed_code = True
#                         print('Went though (%i) iterations. At least half of the required iterations (requested nloop = %i) have finished, code proceed, but careful on the final results'%(i_iter, nloops))
                        
#                     else:
#                         proceed_code = False
#                         print('The code has some issues, meet the stop criteria, please check the code before proceeding')
                    
#                 if proceed_code:
#                     #
#                     # calculate dust temperature
#                     # solve the dust temperature in the last iteration
#                     # to make sure it is the correct Tdust for this density
#                     # os.system("radmc3d mctherm")
#                     get_new_dust_temperature(mint, command=command_mctherm, bool_clean_dusttemp_midplane=mint.bool_clean_dusttemp_midplane)
    
#     else: 
#         """For the case without using VHSE but want dust settling"""
#         if mint.bool_dust_settling & proceed_code:
#             #
#             # solve the dust settling
#             #
#             set_dust_settling(mint, verbose=True)
            
#             # the dust density for the previous iteration
#             filename0 = "dust_density_mdust_GS_" + "final_before_dustset.dat"
#             os.system("cp dust_density.inp " + filename0)
#             # the dust temperature for the previous iteration
#             filename1 = "dust_temperature_mdust_GS_" + "final_before_dustset.dat"
#             os.system("cp dust_temperature.dat " + filename1)
#             # the gas temperature for the previous iteration
#             filename2 = "gas_temperature_mdust_GS_" + "final_before_dustset.inp"
#             os.system("cp gas_temperature.inp " + filename2)
#             # the ratio g2d grid
#             if os.path.exists('ratio_g2d_grid.dat'):
#                 filename3 = "ratio_g2d_grid_GS_" + "final_before_dustset.dat"
#                 os.system("cp ratio_g2d_grid.dat " + filename3)

#             print('the files for this iteration are saved as\n%s\n%s\n%s'%(filename0, filename1, filename2))
#             if os.path.exists('ratio_g2d_grid.dat'):
#                 print('%s is also saved'%(filename3))
            
#             # move the new files as the input
#             os.system('mv dust_density_settled.new dust_density.inp')
#             os.system("mv ratio_g2d_grid.new ratio_g2d_grid.dat")
            
#             # then we need to do the thermal calculation again to get the new dust Temperature
#             # os.system("radmc3d mctherm")
#             get_new_dust_temperature(mint, command=command_mctherm, bool_clean_dusttemp_midplane=mint.bool_clean_dusttemp_midplane)
            
#             # then if we need to get the new gas temperature based on this new dust thermal structure
#             # dustdens = radmc3dPy.analyze.readData(dtemp=True, ddens=True, binary=False)
#             # Tdust = dustdens.dusttemp.copy()
#             # Tdust[Tdust < 2.7] = 2.7
#             # solve_gastemp(mint, dustdens, Tdust, bool_output_gastemp=True)
#             # os.system('mv gas_temperature_direct_output.inp gas_temperature.inp')
        
#     bool_frag_andor_drift = mint.bool_dust_fragmentation or mint.bool_dust_radial_drifting
#     if bool_frag_andor_drift & proceed_code:
        
#         #
#         # Save old files and update dust density file
#         #
#         # the dust density for the previous iteration
#         filename0 = "dust_density_mdust_VHSE_" + str(i_iter) + "_before_fragordrift.dat"
#         os.system("cp dust_density.inp " + filename0)
#         # the dust temperature for the previous iteration
#         filename1 = "dust_temperature_mdust_VHSE_" + str(i_iter) + "_before_fragordrift.dat"
#         os.system("cp dust_temperature.dat " + filename1)
#         # the gas temperature for the previous iteration
#         filename2 = "gas_temperature_mdust_VHSE_" + str(i_iter) + "_before_fragordrift.inp"
#         os.system("cp gas_temperature.inp " + filename2)
#         # the ratio g2d grid
#         if os.path.exists('ratio_g2d_grid.dat'):
#             filename3 = "ratio_g2d_grid_VHSE_" + str(i_iter) + "_before_fragordrift.dat"
#             os.system("cp ratio_g2d_grid.dat " + filename3)

#         print('the files for this iteration are saved as\n%s\n%s\n%s'%(filename0, filename1, filename2))
#         if os.path.exists('ratio_g2d_grid.dat'):
#                 print('%s is also saved'%(filename3))
                
#         #
#         # Put the fragmentation/radial drift
#         #
#         get_new_dust_density(mint, 
#                             bool_dust_vhse=False,
#                             bool_dust_settling_in_iteration=False, # both vhse and settling is solved previously
#                             bool_dust_fragmentation_set=mint.bool_dust_fragmentation,
#                             bool_dust_radial_drifting_set=mint.bool_dust_radial_drifting)
        
#         # move the new files as the input
#         os.system("mv dust_density.new dust_density.inp")
#         os.system("mv ratio_g2d_grid.new ratio_g2d_grid.dat")
        
#         # if mint.bool_dust_settling:
#         #     #
#         #     # Save old files and update dust density file
#         #     #
#         #     # the dust density for the previous iteration
#         #     filename0 = "dust_density_mdust_VHSE_" + str(i_iter) + "_wfragordrift_before_settl.dat"
#         #     os.system("cp dust_density.inp " + filename0)
#         #     # the dust temperature for the previous iteration
#         #     filename1 = "dust_temperature_mdust_VHSE_" + str(i_iter) + "_wfragordrift_before_settl.dat"
#         #     os.system("cp dust_temperature.dat " + filename1)
#         #     # the gas temperature for the previous iteration
#         #     filename2 = "gas_temperature_mdust_VHSE_" + str(i_iter) + "_wfragordrift_before_settl.inp"
#         #     os.system("cp gas_temperature.inp " + filename2)
#         #     # the ratio g2d grid
#         #     if os.path.exists('ratio_g2d_grid.dat'):
#         #         filename3 = "ratio_g2d_grid_VHSE_" + str(i_iter) + "_wfragordrift_before_settl.dat"
#         #         os.system("cp ratio_g2d_grid.dat " + filename3)

#         #     print('the files for this iteration are saved as\n%s\n%s\n%s'%(filename0, filename1, filename2))
#         #     if os.path.exists('ratio_g2d_grid.dat'):
#         #             print('%s is also saved'%(filename3))
        
#         #     set_dust_settling(mint, verbose=True)
#         #     # move the new files as the input
#         #     os.system('mv dust_density_settled.new dust_density.inp')
#         #     os.system("mv ratio_g2d_grid.new ratio_g2d_grid.dat")
            
#         # if mint.bool_dust_inner_rim:
#         #     """(devloping)"""
#         #     add_dust_rim(mint, verbose=True)
#         #     os.system('cp dust_add_rim.new dust_density.inp')
#         #     os.system("mv ratio_g2d_grid.new ratio_g2d_grid.dat")
        
#         # then we need to do the thermal calculation again to get the new dust Temperature
#         # os.system("radmc3d mctherm")
#         get_new_dust_temperature(mint, command=command_mctherm, bool_clean_dusttemp_midplane=mint.bool_clean_dusttemp_midplane)
        
#     """add the dust inner rim"""
#     if mint.bool_dust_inner_rim & proceed_code:
#         """(devloping)"""
#         add_dust_rim(mint, verbose=True)
        
#         # the dust density for the previous iteration
#         filename0 = "dust_density_mdust_VHSE_" + "final_before_addrim.dat"
#         os.system("cp dust_density.inp " + filename0)
#         # the dust temperature for the previous iteration
#         filename1 = "dust_temperature_mdust_VHSE_" + "final_before_addrim.dat"
#         os.system("cp dust_temperature.dat " + filename1)
#         # the gas temperature for the previous iteration
#         filename2 = "gas_temperature_mdust_VHSE_" + "final_before_addrim.inp"
#         os.system("cp gas_temperature.inp " + filename2)
#         # the ratio g2d grid
#         if os.path.exists('ratio_g2d_grid.dat'):
#             filename3 = "ratio_g2d_grid_VHSE_" + "final_before_addrim.dat"
#             os.system("cp ratio_g2d_grid.dat " + filename3)
        
#         # move the new files as the input
#         os.system('cp dust_add_rim.new dust_density.inp')
#         os.system("mv ratio_g2d_grid.new ratio_g2d_grid.dat")
        
#         # get the new thermal structure
#         # os.system("radmc3d mctherm")
#         get_new_dust_temperature(mint, command=command_mctherm, bool_clean_dusttemp_midplane=mint.bool_clean_dusttemp_midplane)
    
#     if mint.bool_SED & mint.bool_VHSE & proceed_code:
#         #
#         # Compute the sed
#         #
#         os.system(command_sed)
#         #
#         filename = "spectrum_VHSE_final.out"
#         os.system("cp spectrum.out " + filename)
    
#     elif mint.bool_SED & mint.bool_dust_settling & proceed_code:
#         # Compute the sed
#         os.system(command_sed)
#         #
#         filename = "spectrum_GS_settled.out"
#         os.system("cp spectrum.out " + filename)

#     """
#     #
#     # Call Dust Sublimation SED calculation
#     # 03/30/2022, we don't need to capture the perfect dust sublimation curve
#     # this dust i_iteration_dustsublimation just needs to run once
#     # to find the sublimation radius
#     #
#     # os.system('python wrapper_iteration_dustsublimation.py')

#     # **TODO**
#     # now the upper scripts was put in a stand out file
#     # for future upgrade, put the dust sublimation computation into the disk_density scripts
#     """

#     #
#     # After the above computation, save for chemical network
#     #
#     if mint.bool_chemistry & proceed_code:
#         # os.system("python chemical_initial_setup.py")
#         chemistry_setup(mint, save_name=para.chemical_save_name, bool_same_rc_as_radmc3d=mint.bool_same_rc_as_radmc3d)
#         # os.system("./run_chem.sh")

#         # for IM Lup
#         # G0Hab_set = 3.439e+10 # this was what was used
#         # G0Hab_set = 1.0e11 # this is the one that could give G0Hab_1AU = 1e7 (actually G0Hab_set at 2.5Rsolar is 7.4e10)
        
#         runchemistry(mint, chem_code_dir=mint.chem_code_dir, G0Hab_set=para.G0Hab_set)

#         """
#         After the chemical network, the .chem file will be created
#         and then use the write_LIME_grid.py file to create the LIME grids
#         for running LIME and get the synthetic line flux from LIME
#         """

#     if mint.bool_savemodel & proceed_code: 
#         # save the model to the data directory
#         savemodel(mint, save_dir=os.path.join(para.chemical_save_dir, para.chemical_save_name), model_dir=mint.file_dir)
        
#         print('also saving the parameter files to')
#         print('%s'%(os.path.join(para.chemical_save_dir, para.chemical_save_name)))
#         # NOTE (Developing)
#         mint.parameters.write_parameters_to_csv(para.chemical_save_name+'_parameters', directory=os.path.join(para.chemical_save_dir, para.chemical_save_name), extension='.csv')
        
#     print('----------------------')
#     print('finished the model for')
#     print('%s'%(para.chemical_save_name))
#     print('----------------------')
    
#     return proceed_code

