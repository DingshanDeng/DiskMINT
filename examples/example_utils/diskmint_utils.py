#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author: Dingshan Deng @ University of Arizona
# contact: dingshandeng@arizona.edu
# created: 02/28/2024
# updated: 02/23/2026
# The useful tools for using DiskMINT to model disks

import datetime, os, sys, copy

import numpy as np
import pandas as pd

from astropy.io import fits
import astropy.units as u
import astropy.constants as C
from astropy.convolution import convolve_fft, Gaussian2DKernel

import tqdm

# for interpolate
from scipy import interpolate
from scipy.optimize import curve_fit

import matplotlib
from matplotlib import colorbar,patches,colors,cm
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.axes_grid1 import AxesGrid
import cmasher as cmr

import bettermoments as bm

style = [
#    'seaborn-ticks',
    {
        'figure.dpi': 300,
        'figure.figsize': (3.5, 3.5 / 1.618),
        'font.size': 9,  # 12
        'image.cmap': 'inferno',
        'figure.dpi': 200,
        'font.family': 'serif',
        'font.serif': ['Times', 'Times New Roman'] + plt.rcParams['font.serif'],
        'xtick.top': True,
        'xtick.direction': 'in',
        'ytick.right': True,
        'ytick.direction': 'in',
        'mathtext.fontset': 'cm'
        }]

plt.style.use(style)

#
# specific packages
#

# from disklab import opacity
# import dsharp_opac as opacity

# import extinction

# import radmc3dPy

# from spectral_cube import SpectralCube
# from gofish import imagecube

# sys.path.append('/home/dingshandeng/github/')
# from gofish_1p5 import imagecube

sys.path.append(os.getcwd())
import platform
if platform.node() == 'garnet':
    github_dir = '/home/dingshandeng/github/DiskModeling/' # garnet
elif platform.node() == 'graphite':
    github_dir = "/Users/dingshandeng/github/DiskModeling/" # macos
elif platform.node() == 'fld':
    github_dir = "/tank/data/research/software/DiskMINT/"
else:
    github_dir = "/home/u24/dingshandeng/github/DiskModeling/" # hpc
    
package_position = os.path.join(github_dir, "src")
sys.path.append(package_position)
import diskmint.model as model
import diskmint.disk_density as dd
import diskmint.constants as const


def get_total_flux_from_radial_profile(intensity, radius):    
    # intensity is in the unit of mJy/arcsec^2
    # and radius is in the unit of arcsec
    # F_tot = int Intensity * 2*pi * radius * d_radius
    func_y = intensity * 2*np.pi * radius
    return np.trapz(func_y, radius)


def cal_fluxinteg(x2, y2, xrange):
    # the x2 y2 are the arrays for the velocity and integrated flux per veloctiy, xrange in m/s
    # the results are in the unit of Jy*m/s
    mask = (x2 > xrange[0]) & (x2 < xrange[1])
    y2_inrange = y2[mask]
    x2_inrange = x2[mask]
    dx2        = (x2[1:] - x2[:-1])[0]
    flux_integ = np.sum(y2_inrange * dx2)
    return flux_integ


def flux2luminosity(flux_intg, dpc, restfreq):
    conv_factor = restfreq * u.GHz / C.c
    return (4 * np.pi * flux_intg * conv_factor * u.mJy * u.km/u.s * (dpc * u.pc)**2).to(u.L_sun).value


def asc2df(file_dir, file_type = 'Photometry'):
    # read SVO's ascII files into DataFrames

    f_test = open(file_dir, 'r')
    line = f_test.readline()

    Filter_all = []
    Wavelength_all = []
    Flux_all = []
    
    if file_type == 'Photometry':
        while line:
#             print(line)
#             if line[0] == '#':
#                 print(line)

            if line[0] != '#':
#                 print(line[0:30].replace(' ',''))
                Filter_all.append(line[0:30].replace(' ','')) # Filter names in Synthetic Photometry

#                 print(line[31:56].replace(' ',''))
                Wavelength_all.append(line[31:56].replace(' ',''))

#                 print(line[57:79].replace(' ',''))
                Flux_all.append(line[57:79].replace(' ',''))

#                 print('--------------------------')
#             print(line.split(' '))
            line = f_test.readline()
        
        df_thisGrid = pd.DataFrame({'Filter': Filter_all, \
                            'Wavelength': np.array(Wavelength_all), \
                            'Flux': np.array(Flux_all)}) 
                            # Wavelength in Angstrom
                            # Flux in erg/cm2/s/A
            
    elif file_type == 'Spectroscopy':
        while line:
#             print(line)
#             line = f_test.readline()
#             if line[0] == '#':
#                 print(line)

            if line[0] != '#':

#                 print(line[0:16].replace(' ',''))
                Wavelength_all.append(line[0:16].replace(' ',''))

#                 print(line[17:32].replace(' ',''))
                Flux_all.append(line[17:32].replace(' ',''))

#                 print('--------------------------')
#             print(line.split(' '))
            line = f_test.readline()
        
        df_thisGrid = pd.DataFrame({'Wavelength': np.array(Wavelength_all), \
                            'Flux': np.array(Flux_all)}) 
                            # Wavelength in Angstrom
                            # Flux in erg/cm2/s/A
    
    # both type need the following moving
    df_thisGrid[['Wavelength', 'Flux']] = df_thisGrid[['Wavelength', 'Flux']].apply(pd.to_numeric)
    
    return df_thisGrid

def find_C_mod2obs(wave_new, fluxes_obs, wave_mod, fluxes_mod):
    # all of the input are the original value (not the over Gaia one)
    # fluxes_obs = C*fluxes_Phot where C = (d/D)**2
    # and d is the radius of the photosphere and D is the distance between observer and the star.
    # the C is the mean value calculated from the optical bands
    
    df_2fluxes = pd.DataFrame({'wave':np.array(wave_new), 
                                'OBS': np.array(fluxes_obs)}) #, 
                                # 'model': np.array(fluxes_new)})

    # +++++++++ Interplot +++++++++

    f = interpolate.interp1d(wave_mod,fluxes_mod,kind="slinear")
    # Clip the wave_new range
    df_2fluxes = df_2fluxes[(df_2fluxes['wave'] > np.min(wave_mod) ) & (df_2fluxes['wave'] < np.max(wave_mod) )]
    # wave_new = wave_new[(wave_new > np.min(wave_mod))]
    # Set the interploted fluxes
    fluxes_new=f(df_2fluxes['wave'])

    # +++++++++ find the best C ++++++++
    df_2fluxes.loc[:,'model'] = np.array(fluxes_new)

    df_2fluxes_filt = df_2fluxes[(df_2fluxes['OBS'] > 0) &
                            (df_2fluxes['OBS'] < np.inf) &
                            (df_2fluxes['wave'] < 3)]
    #             print(df_2fluxes)
    #             print(df_2fluxes_filt)

    C_mod2obs = np.mean(df_2fluxes_filt['OBS']/df_2fluxes_filt['model'])
    
#     print(df_2fluxes_filt)
    
#     plt.figure(dpi=150)
#     # plt.scatter(df_2fluxes_filt['wave'], df_2fluxes_filt['OBS']/df_2fluxes_filt['model'])
#     plt.scatter(df_2fluxes_filt['OBS'], C_mod2obs*df_2fluxes_filt['model'])
#     plt.show()


#     plt.figure(dpi=150)
#     ax0 = plt.axes()

#     plt.scatter(df_2fluxes['wave'], df_2fluxes['OBS'])

#     plt.plot(wave_mod, C_mod2obs*fluxes_mod)

#     ax0.set_xscale('log')
#     ax0.set_yscale('log')
#     plt.show()
    
    return C_mod2obs

# C_mod2obs = find_C_mod2obs(wave_new, np.array(fluxes_obs)*fluxes_GAIA_all.value, wave_mod.value, fluxes_mod*np.sum(fluxes_mod_atGAIA))

def massdust_i(amax, amin, rhobulk = 2.5, pla = 3.5):
#     rhobulk = 2.5 # assume in the unit of [g cm^-3]
#     pla     = 3.x # fitted from long wavelengths SED
    
    md = 4*np.pi/3*rhobulk*(amax**(4-pla) - amin**(4-pla))/(4-pla)
    
    return md

def numbdust_i(amax, amin, rhobulk = 2.5, pla = 3.5):
#     rhobulk = 2.5 # assume in the unit of [g cm^-3]
#     pla     = 3.x # fitted from long wavelengths SED
    
    nd = (amax**(1-pla) - amin**(1-pla))/(1-pla)
    
    return nd

def ndsd_i(amax, amin, rhobulk = 2.5, pla = 3.5):
    # rhobulk = 2.5 # assume in the unit of [g cm^-3]
    # pla     = 3.x # fitted from long wavelengths SED
    # ndsd = nd * a**2 = a**(-pla) * a**2 = a**(2-pla)
    # ndsd_i = int ndsd = a**(3-pla)/(3-pla)
    
    ndsd = (amax**(3-pla) - amin**(3-pla))/(3-pla) # I don't need nd, we just need ndsd here.
    
    return ndsd

def a_average(amax, amin, pla = 3.5):
    # pla = 3.5 for MRN
    # pla = 3.0 fitted from long wavelengths SED
    
    a_ave = (amax**(2-pla) - amin**(2-pla))/(2-pla) / ((amax**(1-pla) - amin**(1-pla))/(1-pla))
    
    return a_ave
            
def wrapper_optool_opac(data_dir, dust_opac_name, amin=0.01/1e4, amax=1000./1e4, pla=3.5, nr_dust=1, lmin=0.01/1e4, lmax=100000/1e4, nl=500, dust_species_1='pyr-mg70', mass_frac_1=0.87, dust_spcies_2='c', mass_frac_2=0.13, porosity=0.25):
    """
    wrapper of making dust opacity using optool; The default values are the DIANA standard opacity.
    
    Caution the units for length to optool should be cm,
    but here we input as um so that the optool command inside the function will multiply the factor of 1e4 to all length variables.
    """
    
    ai = np.logspace(np.log10(amin), np.log10(amax), nr_dust + 1)
    
    for i in range(len(ai)-1):
        
        amin_i = ai[i]
        amax_i = ai[i+1]
        
        optool_command = 'optool %s %.2f %s %.2f -p %.2f -a %.2f %.2f %.2f -l %.2f %.2f %i -mie -radmc'%(dust_species_1, mass_frac_1, dust_spcies_2, mass_frac_2, porosity, amin_i*1e4, amax_i*1e4, pla, lmin*1e4, lmax*1e4, nl)
        os.system(optool_command)

        # data_dir = './data/' # where to save
        # dust_opac_name = 'optool_test' # what is the name

        os.system('mv dustkappa.inp %s'%(os.path.join(data_dir, 'dustkappa_'+dust_opac_name+'_%i.inp'%(int(i)))))
    
    amin_i = ai[0]
    amax_i = ai[-1]
    optool_command = 'optool %s %.2f %s %.2f -p %.2f -a %.2f %.2f %.2f -l %.2f %.2f %i -mie -radmc'%(dust_species_1, mass_frac_1, dust_spcies_2, mass_frac_2, porosity, amin_i*1e4, amax_i*1e4, pla, lmin*1e4, lmax*1e4, nl)
    os.system(optool_command)
    os.system('mv dustkappa.inp %s'%(os.path.join(data_dir, 'dustkappa_'+dust_opac_name+'_combined.inp')))

    return 1


def wrapper_optool_opac_test(data_dir, dust_opac_name, amin=0.01/1e4, amax=1000./1e4, pla=3.5, nr_dust=1, lmin=0.01/1e4, lmax=100000/1e4, nl=500, dust_species_1='pyr-mg70', mass_frac_1=0.87, dust_spcies_2='c', mass_frac_2=0.13, porosity=0.25):
    """
    **NOTE** Developing version
    
    here we test what if we set the amin to amax to a more seperate distribution instead of continuous.
    
    wrapper of making dust opacity using optool; The default values are the DIANA standard opacity.
    
    Caution the units for length to optool should be cm,
    but here we input as um so that the optool command inside the function will multiply the factor of 1e4 to all length variables.
    """
    
    ai = np.logspace(np.log10(amin), np.log10(amax), nr_dust + 1)
    
    for i in range(len(ai)-1):
         
        amin_i = ai[i]
        amax_i = ai[i+1]
        
        optool_command = 'optool %s %.2f %s %.2f -p %.2f -a %.2f %.2f %.2f -l %.2f %.2f %i -mie -radmc'%(dust_species_1, mass_frac_1, dust_spcies_2, mass_frac_2, porosity, amin_i*1e4, amax_i*1e4, pla, lmin*1e4, lmax*1e4, nl)
        os.system(optool_command)

        # data_dir = './data/' # where to save
        # dust_opac_name = 'optool_test' # what is the name

        os.system('mv dustkappa.inp %s'%(os.path.join(data_dir, 'dustkappa_'+dust_opac_name+'_%i.inp'%(int(i)))))
    
    amin_i = ai[0]
    amax_i = ai[-1]
    optool_command = 'optool %s %.2f %s %.2f -p %.2f -a %.2f %.2f %.2f -l %.2f %.2f %i -mie -radmc'%(dust_species_1, mass_frac_1, dust_spcies_2, mass_frac_2, porosity, amin_i*1e4, amax_i*1e4, pla, lmin*1e4, lmax*1e4, nl)
    os.system(optool_command)
    os.system('mv dustkappa.inp %s'%(os.path.join(data_dir, 'dustkappa_'+dust_opac_name+'_combined.inp')))

    return 1


def massdust_i(amax, amin, rhobulk = 2.5, pla = 3.5):
    # rhobulk = 2.5 # assume in the unit of [g cm^-3]
    # pla     = 3.x # fitted from long wavelengths SED

    md = 4*np.pi/3*rhobulk*(amax**(4-pla) - amin**(4-pla))/(4-pla)

    return md


def numbdust_i(amax, amin, rhobulk = 2.5, pla = 3.5):
    # rhobulk = 2.5 # assume in the unit of [g cm^-3]
    # pla     = 3.x # fitted from long wavelengths SED

    nd   = (amax**(1-pla) - amin**(1-pla))/(1-pla)

    return nd


def ndsd_i(amax, amin, rhobulk = 2.5, pla = 3.5):
    # rhobulk = 2.5 # assume in the unit of [g cm^-3]
    # pla     = 3.x # fitted from long wavelengths SED
    # ndsd = nd * a**2 = a**(-pla) * a**2 = a**(2-pla)
    # ndsd_i = int ndsd = a**(3-pla)/(3-pla)

    ndsd = (amax**(3-pla) - amin**(3-pla))/(3-pla) # I don't need nd, we just need ndsd here.

    return ndsd


def a_average(amax, amin, pla = 3.5):
    # pla = 3.5 for MRN
    # pla = 3.x fitted from long wavelengths SED

    a_ave = (amax**(2-pla) - amin**(2-pla))/(2-pla) / ((amax**(1-pla) - amin**(1-pla))/(1-pla))

    return a_ave


def set_dust_parameters_to_nrdust_1n2(para, amin_set_list, amax_set_list, nr_dust_set_list, pla_set, rhobulk_set = 3.224):
    """
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
    """
    #
    #
    # number of dust species
    # to setup our VHSE models
    # [NOTE] here we support two different dust species
    # which can be further divided into different subgroups
    # with different sizes.
    para.nr_dust_1         = nr_dust_set_list[0]
    para.nr_dust_2         = nr_dust_set_list[1]
    # automatically set the total nr of dust
    para.dust_spec_nr      = para.nr_dust_1 + para.nr_dust_2

    #
    # use following parameters to make dust kappa files
    #
    # here a is the radius of dust grain in the unit of cm.
    #
    # minimum size of the dust
    # use 100 AA as 1.0e-6 so that dust is not as small as molecule
    para.amin_1             = amin_set_list[0] # 1.0e-6
    para.amin_2             = amin_set_list[1] # 1.0e-6
    #
    # maximum size of the dust
    para.amax_1             = amax_set_list[0]
    para.amax_2             = amax_set_list[1]
    #
    # automatically set the overall smallest
    # and largest sizes among the two
    para.amin_all           = np.nanmin([para.amin_1, para.amin_2])
    para.amax_all           = np.nanmax([para.amax_1, para.amax_2])
    #
    # power law index of the dust size distribution
    para.pla_dustsize       = pla_set

    #
    # bulk density
    # the bulk density should be larger instead of this value
    # but it seems like it is going to be cancelled so double check the code
    # the rhobulk is used in chemical network
    # prev used 2.5 (too small)
    #
    para.rhobulk            = rhobulk_set   # assume in the unit of [g cm^-3]
    
    return para
    
        
def set_dust_parameters_to_nrdust1(para, amin_set, amax_set, nr_dust_set, pla_set, rhobulk_set = 3.224):
    """
    Set the dust parameters info to dust species 1, this is used when you only want to set up 1 dust specy in the model.
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
    """
    #
    #
    # number of dust species
    # to setup our VHSE models
    # [NOTE] here we support two different dust species
    # which can be further divided into different subgroups
    # with different sizes.
    para.nr_dust_1         = nr_dust_set
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
    para.amin_1             = amin_set # 1.0e-6
    para.amin_2             = np.nan # 1.0e-6
    #
    # maximum size of the dust
    para.amax_1             = amax_set
    para.amax_2             = np.nan
    #
    # automatically set the overall smallest
    # and largest sizes among the two
    para.amin_all           = np.nanmin([para.amin_1, para.amin_2])
    para.amax_all           = np.nanmax([para.amax_1, para.amax_2])
    #
    # power law index of the dust size distribution
    para.pla_dustsize       = pla_set

    #
    # bulk density
    # the bulk density should be larger instead of this value
    # but it seems like it is going to be cancelled so double check the code
    # the rhobulk is used in chemical network
    # prev used 2.5 (too small)
    #
    para.rhobulk            = rhobulk_set   # assume in the unit of [g cm^-3]
    
    return para


# def read_radmc3d_SED(specname, distance, index = 'nufnu'):
#     #
#     # index can be either nufnu or fnu
#     #
    
#     s     = radmc3dPy.analyze.readSpectrum(specname)
#     lam   = s[:,0]
#     nu    = 1e4*radmc3dPy.natconst.cc/lam
#     fnu   = s[:,1] # in the cgs units of [erg cm**-2 s**-1 Hz**-1]
#     # nufnu_1pc = nu*fnu # the model seen from 1pc
#     # scale according to its distance
#     fnu_here = (fnu * (1/distance)**2)
    
#     if index == 'nufnu':
#         outputy = nu * fnu_here
#     elif index == 'fnu':
#         outputy = fnu_here
#     else:
#         print('please choose a index either nufnu or fnu')
#         print('here assumed to out put fnu_here')
#         outputy = fnu_here
        
#     return lam, outputy


def read_model(chem_save_name, data_dir_init, data_dir_chemical, have_CO2):
    #
    # Read in the CO init Grid
    #
    cord_t = 'cyl' 
    file_readin = os.path.join(data_dir_init, 'COinitgrid-' + chem_save_name + '-' + cord_t + '.dat')
    print('reading in the init file %s'%(file_readin))
    df_COinit_1 = pd.read_csv(file_readin, skiprows=3) 

    # Tgas_ar = np.array(df_COinit_1['Tgas[K]'].values, dtype = np.float64)

    #
    # Read in the CO end Grid
    #
    if have_CO2 in [True, False]:
        if have_CO2:
            file_readin = os.path.join(data_dir_chemical, 'COendgrid-withCO2-' + chem_save_name + '-' + cord_t + '.chem')
        if not have_CO2:
            file_readin = os.path.join(data_dir_chemical, 'COendgrid-noCO2-' + chem_save_name + '-' + cord_t + '.chem')
    else:
        file_readin = os.path.join(data_dir_chemical, 'COendgrid-'+str(have_CO2)+'-' + chem_save_name + '-' + cord_t + '.chem')
        
    print("reading in the chem file %s"%(file_readin))
    xc18o_1_chemical = np.loadtxt(file_readin)

    #
    # Combine and save the LIME grid
    #
    df_COinit_comb = df_COinit_1.copy()
    df_COinit_comb['r-chem[cm]'] = xc18o_1_chemical[:, 0]
    df_COinit_comb['z-chem[cm]'] = xc18o_1_chemical[:, 1]
    df_COinit_comb['lgX(H2)']    = xc18o_1_chemical[:, 2]
    df_COinit_comb['lgX(C18O)']  = xc18o_1_chemical[:, 3]
    
    # print(df_COinit_1['n_H[cm**-3]'].values, df_COinit_comb['n_H[cm**-3]'].values)
    
    return df_COinit_comb

def name_modelpara(df_COinit_comb, nr_LIME = 100, nz_LIME = 200):
    #
    # CO init Grid
    #
    r_init_ar  = df_COinit_comb['r[cm]'].values
    zr_init_ar = df_COinit_comb['z/r'].values
    nH_1_init_ar = df_COinit_comb['n_H[cm**-3]'].values
    Tgas_1_init_ar = df_COinit_comb['Tgas[K]'].values

    rr_2D  = np.reshape(r_init_ar, (nr_LIME, nz_LIME))
    zrr_2D = np.reshape(zr_init_ar, (nr_LIME, nz_LIME))

    rr_grid = rr_2D[:,0]
    zrr_grid = zrr_2D[0,:]

    nH_1_init_2D  = np.reshape(nH_1_init_ar, (nr_LIME, nz_LIME))
    Tgas_1_init_2D  = np.reshape(Tgas_1_init_ar, (nr_LIME, nz_LIME))
    
    tauv_star_2D = np.reshape(df_COinit_comb['tauv_star'].values, (nr_LIME, nz_LIME))
    tauv_zup_2D  = np.reshape(df_COinit_comb['tauv_zup'].values, (nr_LIME, nz_LIME))
    #
    # CO end grid
    #
    r_end_ar  = df_COinit_comb['r-chem[cm]'].values
    z_end_ar  = df_COinit_comb['z-chem[cm]'].values
    lgXH2_1_end_ar = df_COinit_comb['lgX(H2)'].values
    lgXC18O_1_end_ar = df_COinit_comb['lgX(C18O)'].values

    rr_end_2D  = np.reshape(r_end_ar, (nr_LIME, nz_LIME))
    zz_end_2D  = np.reshape(z_end_ar, (nr_LIME, nz_LIME))

    rr_end_grid = rr_end_2D[:,0]
    zz_end_grid = zz_end_2D[0,:]

    lgXH2_1_end_2D  = np.reshape(lgXH2_1_end_ar, (nr_LIME, nz_LIME))
    lgXC18O_1_end_2D  = np.reshape(lgXC18O_1_end_ar, (nr_LIME, nz_LIME))
    
    # print(lgXC18O_1_end_2D)
    
    # Combine the two to show its distribution each molecule
    nH2_1_end_2D = nH_1_init_2D * 10**lgXH2_1_end_2D
    nC18O_1_end_2D = nH_1_init_2D * 10**lgXC18O_1_end_2D
    
    return (rr_2D, zrr_2D), (rr_grid, zrr_grid), (tauv_star_2D, tauv_zup_2D), nH_1_init_2D, Tgas_1_init_2D, nH2_1_end_2D, nC18O_1_end_2D

def compute_emittinglayer(df_COinit_comb, dv_set = 0.0, nr = 100, nz = 200, save_file=True):
    
    # setup dv values
    # dv_set = 0.0 # 1.6 # km/s
    
    au = 1.496e13; 
    lsol = 3.9e33; 
    clight = 2.99e10; 
    hplanck = 6.626e-27
    molmass=30.0 # mass of molecule C18O
    
    r    = df_COinit_comb['r[cm]'].values.reshape(nr,nz)/au
    z    = df_COinit_comb['z[cm]'].values.reshape(nr,nz)/au
    n    = df_COinit_comb['n_H[cm**-3]'].values.reshape(nr,nz)
    tg   = df_COinit_comb['Tgas[K]'].values.reshape(nr,nz)
    xco  = 10**df_COinit_comb['lgX(C18O)'].values.reshape(nr,nz)
    
    g = lambda j : 2*j + 1 # statistical weight function

    A21  = 6.011e-07; nu_21 = 219.560e9
    tex2 = 15.81; tex1 = 5.27

    A32  = 2.172e-06; nu_32 = 329.3305e9
    tex3 = 31.61
    
    w21  = clight/nu_21 /1.0e-4 # in microns
    w32  = clight/nu_32 /1.0e-4 # in microns
    
    if save_file:
        os.system('rm lineop')
        wfile = open('lineop','w+')

    # ESCAPE PROBABILITY (Eq. B9 of Tielens+Hollenbach 1985

    def beta(tau_arr):
        n=np.size(tau_arr)
        result = np.zeros(n)
        for i in range(n):
            tau = tau_arr[i] + 1.0e-30
            if tau > 7:
                 result[i] = 1.0/(4.0*tau*np.sqrt(np.log(tau/(np.sqrt(np.pi)))))
            elif tau > 1.0e-10:
                 result[i] =  (1.0-np.exp(-2.34*tau))/(4.68*tau)
            else:
                 result[i] = 0.5 # above expression underflows..

        return result

    # C18O LUMINOSITY 

    c18olum21 = 0.0
    c18olum32 = 0.0

    wholedata = []

    for i in np.arange(0, nr): # nr
        if i==0: 
            dr = 0.0
        else:
            dr = au * (r[i,0] - r[i-1,0])

        dz = au * np.insert(np.diff(z[i,:]),0,z[i,0])
        dvol = 2.0*np.pi * r[i,0]*au * dr * dz # volume element

        nco = xco[i,:] * n[i,:]  # C18O number density
        zpf = 2.76 * tg[i,:] # partition function [approximate!]

        # level populations
        nco1 = nco * g(1)/zpf * np.exp(-tex1/tg[i,:])
        nco2 = nco * g(2)/zpf * np.exp(-tex2/tg[i,:])
        nco3 = nco * g(3)/zpf * np.exp(-tex3/tg[i,:])

        dvdop = dv_set + 0.13 * np.sqrt(tg[i,:]/molmass) # line width assumed thermal

    # FIRST the 2-1 line

        # tau of line in gridcells
        dtau = 2.2472e-19*A21*w21**3.0*(nco1*g(2)/g(1) -nco2) * dz/dvdop
        dtau = dtau + 1.0e-99
        tauup = np.sum(dtau) - np.cumsum(dtau) # tau in the up direction
        taudn = 2* np.sum(dtau) - tauup # tau in the down direction

        # escape in up and down directions
        betaup = beta(tauup); betadn = beta(taudn)

        # luminosity from this gridcell
        dlum21 = (nco2 * dvol) * A21 * hplanck*nu_21 * (betaup+betadn)

    # THEN the 3-2 line

        # tau of line in gridcells
        dtau  = 2.2472e-19*A32*w32**3.0*(nco2*g(3)/g(2) -nco3) * dz/dvdop
        dtau  = dtau + 1.0e-99
        tauup = np.sum(dtau) - np.cumsum(dtau) # tau in the up direction
        taudn = 2* np.sum(dtau) - tauup # tau in the down direction

        # escape in up and down directions
        betaup = beta(tauup); betadn = beta(taudn)

        # luminosity from this gridcell
        dlum32 = (nco3 * dvol) * A32 * hplanck*nu_32 * (betaup+betadn)

        writinglines = np.c_[r[i,:],z[i,:],dlum21,dlum32,tauup,taudn]

        if len(wholedata) == 0:
            wholedata = writinglines
        else:
            wholedata = np.vstack([wholedata, writinglines])

        writinglines[np.isnan(writinglines)] = 0

        if len(writinglines) != nz:
            print('this length does not equal to nz (%i)'%(nz))
            print('length:', len(writinglines))

        if save_file:
            np.savetxt(wfile, writinglines, fmt='%1.3e')

        c18olum21 = c18olum21 + np.sum(dlum21)
        c18olum32 = c18olum32 + np.sum(dlum32)

    print('C18O 2-1 %.3e Lsolar'%(2*c18olum21/lsol)) # Two sides of disk
    print('C18O 3-2 %.3e Lsolar'%(2*c18olum32/lsol)) # Two sides of disk
    
    return wholedata


def read_emittinglayer(wholedata=None, data_dir_emit = './', nr = 100, nz = 200, nr_skiprows = 0):
    
    if wholedata is None: 
        xc18o_1_emit = np.loadtxt(data_dir_emit + 'lineop') 
    else:
        xc18o_1_emit = wholedata
    
#     r_au_flat = xc18o_1_emit[nr_skiprows:, 0]
#     z_au_flat = xc18o_1_emit[nr_skiprows:, 1]
#     dlum21_flat = xc18o_1_emit[nr_skiprows:, 2]
#     dlum32_flat = xc18o_1_emit[nr_skiprows:, 3]
#     tauup_flat = xc18o_1_emit[nr_skiprows:, 4]

    r_au_2D = np.reshape(xc18o_1_emit[nr_skiprows:, 0], (nr, nz))
    z_au_2D = np.reshape(xc18o_1_emit[nr_skiprows:, 1], (nr, nz))
    dlum21_2D = np.reshape(xc18o_1_emit[nr_skiprows:, 2], (nr, nz))
    dlum32_2D = np.reshape(xc18o_1_emit[nr_skiprows:, 3], (nr, nz))
    
    tauup_2D = np.reshape(xc18o_1_emit[nr_skiprows:, 4], (nr, nz))
    
    dlum21_2D_adjusted = dlum21_2D.copy()
    dlum21_2D_adjusted[dlum21_2D_adjusted == 0] = 1e-99

    dlum32_2D_adjusted = dlum32_2D.copy()
    dlum32_2D_adjusted[dlum32_2D_adjusted == 0] = 1e-99
    
    return r_au_2D, z_au_2D, dlum21_2D_adjusted, dlum32_2D_adjusted, tauup_2D


def fits_reader(filein, rest_freq = 219.5603541 * u.GHz, vlsr = 4.5, datatype = 'MOD'):
    """
    simply read in a LIME output fits file.
    datatype can be chosen from 'MOD' or 'OBS'
    rest_freq = 219.5603541 * u.GHz # C18O
     vlsr = 4.5 # km/s for RU Lup
    """
    
    ############################################################
    hdul = fits.open(filein)

    # Get the RA, DEC and channel values from the header
    # The axis 1 is the RA, the axis 2 is the DEC and the
    # axis 3 is the velocity channels

    header = hdul[0].header    # Primary HDU header
    naxis1 = header['naxis1']  # Number of values in each axis
    naxis2 = header['naxis2']  # naxis1 and naxis2 should be equal
    naxis3 = header['naxis3']

    cdelt1 = header['cdelt1']  # Difference between two pixel
    cdelt2 = header['cdelt2']  # cdelt1 and cdelt2 should be opposite
    cdelt3 = header['cdelt3']  # cdelt1 negative

    crpix1 = header['crpix1']  # Reference pixel
    crpix2 = header['crpix2']  # All of them should be 0.
    crpix3 = header['crpix3']

    # crval1 = header['crval1']  # Value at the reference pixel
    # crval2 = header['crval2']
    crval1 = 0  # Value at the reference pixel set as 0 to set the delta x and y values
    crval2 = 0
    crval3 = header['crval3']

    # Compute the values of RA and DEC

    ra = np.arange(naxis1, dtype = np.float64)
    dec = np.arange(naxis2, dtype = np.float64)
    veloc = np.arange(naxis3, dtype = np.float64)
    for i in range(naxis1):
        ra[i] = (i + 1 - crpix1) * cdelt1 + crval1
    for i in range(naxis2):
        dec[i] = (i + 1 - crpix2) * cdelt2 + crval2
    for i in range(naxis3):
        veloc[i] = (i + 1 - crpix3) * cdelt3 + crval3

    # The values of RA and DEC are in degrees in the fits file
    # Velocities are in m/s. We convert them in arcsec and km/s

    if datatype == 'MOD':
        ra = ra * 3600 # arcsec 
        dec = dec * 3600 # arcsec
        veloc = veloc / 1e3 # km/s

    elif datatype == 'OBS':
        freq = veloc
        ra = ra * 3600 # arcsec 
        dec = dec * 3600 # arcsec

        veloc = C.c.to(u.km/u.s).value * (freq - rest_freq.to(u.Hz).value) / rest_freq.to(u.Hz).value # km/s
        veloc = -(veloc + vlsr) # fix the total velocity for this source
    # Fix keyword for velocity
    # header['ctype3'] = 'VELO-LSR'

    # Now read the convolved image data
    img = hdul[0].data

    moment = 0

    # _, nvkms, nx, ny = img.shape

    img_t = np.where(np.isfinite(img), img, 0.0)
    img_t = img_t[0, :, :, :]
    
    return img_t, ra, dec, veloc


def plot_density_nH_n_Tgas(fig, axs, rr_grid_model0, zrr_grid_model0, nH_1_init_2D_model0, Tgas_1_init_2D_model0, tauv_star_2D_model0, lim_nH=np.arange(4, 18.5, 0.5), ticksrange_nH = np.arange(4, 20, 2), lim_Tgas = np.arange(0, 3.6, 0.1), ticksrange_Tgas = np.arange(0, 4.0, 0.5), colormap_selected='viridis', titlelocator = [0.05, 0.95]):
    
    ax_t = axs[0]
    # ax_t.tick_params(bottom='on', top='on', right='on', left='on', which='major')
    # ax_t.tick_params(bottom='on', top='on', right='on', left='on', which='minor')

    ax_t.set_xticks([0.5, 1, 10, 100])
    ax_t.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])

    # ax_t.set_title(name_infigure_model0)

    # ax_t.text(titlelocator[0], titlelocator[1], r'$\log_{10} n({\mathrm{H}})\,\mathrm{[cm^{-3}]}$', fontsize = 18, color = 'white')
    ax_t.annotate(r'$\log_{10} n({\mathrm{H}})\,\mathrm{[cm^{-3}]}$', xy=(titlelocator[0], titlelocator[1]), xycoords='axes fraction', ha='left', va='top', fontsize = 18, color = 'white')

    # # ax_t.text(100, 0.2, 'dust scatter light', color = 'brown', fontsize = 10, ha='right')
    # ax_t.annotate('dust scatter light', xy = (100, 0.15), xytext = (200, 0.25), \
    #               color = 'brown', fontsize = 10, ha='right',\
    #              arrowprops=dict(arrowstyle = "->", color = 'brown'))

    c = ax_t.contourf(rr_grid_model0/const.au, zrr_grid_model0, np.log10((nH_1_init_2D_model0).T), lim_nH, cmap=colormap_selected, extend = 'both', zorder = -1)

    cb = fig.colorbar(c, ax=ax_t, ticks=ticksrange_nH)
    # cb.set_label(r'$\log_{10}{n_{\mathrm{gas}}}$', labelpad = 15, y = 0.5, rotation=270., fontsize = 12)

    ctau = ax_t.contour(rr_grid_model0/const.au, zrr_grid_model0, np.log10(tauv_star_2D_model0.T), [0, 1], \
                        colors='#C7C7C7', linestyles='solid', zorder = 10)
    # manual_locations = [(1, 1.1), (5, 0.9), (5, 0.5)] 
    # plt.clabel(ctau, inline=0, inline_spacing=1, colors = 'k', fontsize=9, fmt = r'$\log_{10}(A_V) = $'+'%1.0f', rightside_up = True, manual = manual_locations)

    # ax_t.set_xlabel('r [AU]', fontsize = 12)
    ax_t.set_ylabel('z/r', fontsize = 12)
    
    # Initialize minor ticks
    ax_t.minorticks_on()

    ###
    # 1, 0 Tgas VHSE
    ###

    ax_t = axs[1]

    # ax_t.set_title('')

    # ax_t.text(titlelocator[0], titlelocator[1], r'$\log_{10} T_{\mathrm{gas}}\,\mathrm{[K]}$', fontsize = 18, color = 'white')
    ax_t.annotate(r'$\log_{10} T_{\mathrm{gas}}\,\mathrm{[K]}$', xy=(titlelocator[0], titlelocator[1]), xycoords='axes fraction', ha='left', va='top', fontsize = 18, color = 'white')

    c = ax_t.contourf(rr_grid_model0/const.au, zrr_grid_model0, np.log10((Tgas_1_init_2D_model0).T), lim_Tgas, cmap=colormap_selected, extend = 'both', zorder = -1)

    cb = fig.colorbar(c, ax=ax_t, ticks=ticksrange_Tgas)
    # cb.set_label(r'$\log_{10}{n_{\mathrm{gas}}}$', labelpad = 15, y = 0.5, rotation=270., fontsize = 12)

    # ctauTgas = ax_t.contour(rr_grid_model0/const.au, zrr_grid_model0, Tgas_1_init_2D_model0.T, [20, 30], \
    #                     colors='k', linestyles='--', zorder = 10)

    ctau = ax_t.contour(rr_grid_model0/const.au, zrr_grid_model0, np.log10(tauv_star_2D_model0.T), [0, 1], \
                        colors='#C7C7C7', linestyles='solid', zorder = 10)
    # manual_locations = [(1, 1.1), (5, 0.9), (5, 0.5)] 
    # plt.clabel(ctau, inline=0, inline_spacing=1, colors = 'k', fontsize=9, fmt = r'$\log_{10}(A_V) = $'+'%1.0f', rightside_up = True, manual = manual_locations)

    ax_t.set_xlabel('r [AU]', fontsize = 12)
    ax_t.set_ylabel('z/r', fontsize = 12)

    ax_t.minorticks_on()
    
    return 1


def plot_density_emitting_layer(fig, axs, rr_grid_model0, zrr_grid_model0, nH_1_init_2D_model0, Tgas_1_init_2D_model0, tauv_star_2D_model0, nC18O_1_end_2D_model0,dlum21_2D_model0,dlum32_2D_model0=None,x_scale='linear',temp_scale='linear', name_infigure='',
                                limymin = 0, limymax = 0.7,
                                limxmin = 0.035, limxmax = 500,
                                titlelocator = [0.05, 0.95],
                                lim_nH = np.arange(4, 12.1, 0.1),
                                ticksrange_nH = np.arange(4, 14, 2),
                                lim_Tgas = np.arange(10, 80.1, 0.1),
                                ticksrange_Tgas = np.arange(10, 90, 10),
                                lim_Tgas_log = np.arange(0, 3.1, 0.1),
                                ticksrange_Tgas_log = np.arange(0, 3.5, 0.5),
                                lim_nC18O = np.arange(-1, 4.05, 0.05),
                                ticksrange_nC18O = np.arange(-1, 5, 1),
                                ticksrange_lum21 = np.arange(-4.0, -0.5, 0.5),
                                ticksrange_lum32 = np.arange(-4.0, -0.5, 0.5),
                                limstep = 0.05,
                                textsize = 16,
                                cmap_selected = 'viridis',
                                bool_3_2=False,
                                bool_set_title=False):

    lim_lum21 = np.arange(np.min(ticksrange_lum21), np.max(ticksrange_lum21) + limstep, limstep)
    lim_lum32 = np.arange(np.min(ticksrange_lum32), np.max(ticksrange_lum32) + limstep, limstep)
    
    #
    # Clean up the tauz map
    #

    # tauup_CO_2D_model1_t = tauup_CO_2D_model1.copy()
    # tauup_CO_2D_model1_t[nC18O_1_end_2D_model1 < 1e-2] = np.nan

    # tauup_CO_2D_model0_t = tauup_CO_2D_model0.copy()
    # tauup_CO_2D_model0_t[nC18O_1_end_2D_model0 < 1e-2] = np.nan
    
    ###
    # 0, 0 nH
    ###

    # ax_t = ax[0, 0]
    ax_t = axs[0]
    
    if bool_set_title:
        ax_t.set_title(r'$\log_{10} n({\mathrm{H}})\,\mathrm{[cm^{-3}]}$', fontsize = textsize)

    if name_infigure != '':
        ax_t.text(titlelocator[0], titlelocator[1], name_infigure, fontsize = textsize, color = 'white', verticalalignment = 'top', transform = ax_t.transAxes)

    c = ax_t.contourf(rr_grid_model0/const.au, zrr_grid_model0, np.log10((nH_1_init_2D_model0).T), lim_nH, cmap = cmap_selected, extend = 'both', zorder = -1)

    cb = fig.colorbar(c, ax=ax_t, ticks=ticksrange_nH)
    # cb.set_label(r'$\log_{10}{n_{\mathrm{gas}}}$', labelpad = 15, y = 0.5, rotation=270., fontsize = 12)

    ctau = ax_t.contour(rr_grid_model0/const.au, zrr_grid_model0, np.log10(tauv_star_2D_model0.T), [0, 1], \
                        colors='white', linestyles='solid', zorder = 10)
    # manual_locations = [(50, 0.35), (50, 0.25)] 
    # plt.clabel(ctau, inline=0, inline_spacing=1, colors = 'white', fontsize=9, fmt = r'$\log_{10}(A_V) = $'+'%1.0f', rightside_up = True, manual = manual_locations)

    # ctau = ax_t.contour(rr_grid_model0/const.au, zrr_grid_model0, tauv_star_2D_model0.T, [1, 10], \
    #                     colors='white', linestyles='solid', zorder = 10)
    # manual_locations = [(80, 0.35), (80, 0.25)] 
    # plt.clabel(ctau, inline=0, inline_spacing=0., colors = 'white', fontsize=9, fmt = r'$A_V = $'+'%1.0f', rightside_up = True, manual = manual_locations)

    #
    # The text location for Av_zup lines
    #
    # ax_t.text(70, 0.27, r'$A_V = $'+'%1.0f' % (1), zorder = 11, color='white', fontsize=9)
    # ax_t.text(70, 0.17, '%1.0f' % (10), zorder = 11, color='white', fontsize=9)

    ctau2 = ax_t.contour(rr_grid_model0/const.au, zrr_grid_model0, np.log10(nC18O_1_end_2D_model0.T), [-1], \
                        colors='k', linestyles='--', zorder = 10)

    # ax_t.text(100, 2.5, r'$\log_{10}(A_V) = $'+'%1.0f' % (-2), zorder = 11)
    # ax_t.text(100, 2.0, '%1.0f' % (0), zorder = 11)
    # ax_t.text(100, 1.7, '%1.0f' % (1), zorder = 11)

    # ax_t.set_xlabel('r [AU]', fontsize = 12)
    ax_t.set_ylabel('z/r', fontsize = 12)

    ax_t.set_xlim(limxmin, limxmax)
    ax_t.set_ylim(limymin, limymax)

    # ax_t.set_xscale('log')
    # Initialize minor ticks
    ax_t.minorticks_on()

    ###
    # 0, 1 Tgas
    ###

    ax_t = axs[1]
   
    
    if bool_set_title:
        if temp_scale == 'log':
            ax_t.set_title(r'$\log_{10} T_{\mathrm{gas}}\,\mathrm{[K]}$', fontsize = textsize)
        else:
            ax_t.set_title(r'$T_{\mathrm{gas}}\,\mathrm{[K]}$', fontsize = textsize)

    if temp_scale == 'log':
        c = ax_t.contourf(rr_grid_model0/const.au, zrr_grid_model0, (np.log10(Tgas_1_init_2D_model0)).T, lim_Tgas_log, cmap = cmap_selected, extend = 'both', zorder = -1)
        cb = fig.colorbar(c, ax=ax_t, ticks=ticksrange_Tgas_log)
        
    else:
        c = ax_t.contourf(rr_grid_model0/const.au, zrr_grid_model0, (Tgas_1_init_2D_model0).T, lim_Tgas, cmap = cmap_selected, extend = 'both', zorder = -1)
        cb = fig.colorbar(c, ax=ax_t, ticks=ticksrange_Tgas)
        
    # cb.set_label(r'$\log_{10}{n_{\mathrm{gas}}}$', labelpad = 15, y = 0.5, rotation=270., fontsize = 12)

    # ctau = ax_t.contour(rr_grid_model0/const.au, zrr_grid_model0, np.log10(tauup_CO_2D_model0_t.T), [0, 1], \
    #                     colors='k', linestyles='solid', zorder = 10)
    # # manual_locations = [(1, 1.1), (5, 0.9), (5, 0.5)] 
    # # plt.clabel(ctau, inline=0, inline_spacing=1, colors = 'k', fontsize=9, fmt = r'$\log_{10}(A_V) = $'+'%1.0f', rightside_up = True, manual = manual_locations)

    ctau = ax_t.contour(rr_grid_model0/const.au, zrr_grid_model0, np.log10(tauv_star_2D_model0.T), [0, 1], \
                        colors='white', linestyles='solid', zorder = 10)

    ctau2 = ax_t.contour(rr_grid_model0/const.au, zrr_grid_model0, np.log10(nC18O_1_end_2D_model0.T), [-1], \
                        colors='k', linestyles='--', zorder = 10)

    # ax_t.text(100, 2.5, r'$\log_{10}(A_V) = $'+'%1.0f' % (-2), zorder = 11)
    # ax_t.text(100, 2.0, '%1.0f' % (0), zorder = 11)
    # ax_t.text(100, 1.7, '%1.0f' % (1), zorder = 11)

    # ax_t.set_xlabel('r [AU]', fontsize = 12)
    # ax_t.set_ylabel('z/r', fontsize = 12)

    ax_t.set_xlim(limxmin, limxmax)
    ax_t.set_ylim(limymin, limymax)

    # ax_t.set_xscale('log')
    # Initialize minor ticks
    ax_t.minorticks_on()

    ###
    # 0, 2 nC18O w/ CO2
    ###

    ax_t = axs[2]
    
    if bool_set_title:
        ax_t.set_title(r'$\log_{10} n({\mathrm{C^{18}O}})\,\mathrm{[cm^{-3}]}$', fontsize = textsize)

    # ax_t.text(titlelocator[0], titlelocator[1], r'$\log_{10} n({\mathrm{C^{18}O}})\,\mathrm{[cm^{-3}]}$', fontsize = textsize, color = 'white', verticalalignment = 'top', transform = ax_t.transAxes)

    nC18O_1_end_2D_model0_t = nC18O_1_end_2D_model0.copy()
    nC18O_1_end_2D_model0_t[nC18O_1_end_2D_model0_t == 0] = 1e-199

    c = ax_t.contourf(rr_grid_model0/const.au, zrr_grid_model0, np.log10((nC18O_1_end_2D_model0_t).T), lim_nC18O, cmap = cmap_selected, extend = 'both', zorder = -1)

    cb = fig.colorbar(c, ax=ax_t, ticks=ticksrange_nC18O)
    # cb.set_label(r'$\log_{10}{n_{\mathrm{gas}}}$', labelpad = 15, y = 0.5, rotation=270., fontsize = 12)

    # ctau = ax_t.contour(rr_grid_model0/const.au, zrr_grid_model0, np.log10(tauup_CO_2D_model0_t.T), [0, 1], \
    #                     colors='k', linestyles='solid', zorder = 10)
    # # manual_locations = [(1, 1.1), (5, 0.9), (5, 0.5)] 
    # # plt.clabel(ctau, inline=0, inline_spacing=1, colors = 'k', fontsize=9, fmt = r'$\log_{10}(A_V) = $'+'%1.0f', rightside_up = True, manual = manual_locations)

    # ax_t.text(100, 1.0, r'$\log_{10}(A_V) = $'+'%1.0f' % (-2), zorder = 11)
    # ax_t.text(50, 0.50, '%1.0f' % (0), zorder = 11)
    # ax_t.text(20, 0.25, '%1.0f' % (1), zorder = 11)

    # ax_t.set_xlabel('r [AU]', fontsize = 12)
    # ax_t.set_ylabel('z/r', fontsize = 12)

    ax_t.set_xlim(limxmin, limxmax)
    ax_t.set_ylim(limymin, limymax)

    # ax_t.set_xscale('log')
    # Initialize minor ticks
    ax_t.minorticks_on()


    ###
    # 0, 3 C18O (2-1) w/ CO2
    ###

    ax_t = axs[3]
    
    if bool_set_title:
        ax_t.set_title(r'$\log_{10} {\mathrm{C^{18}O}}\ (2-1)$', fontsize = textsize)

    # ax_t.text(titlelocator[0], titlelocator[1], r'$\log_{10} {\mathrm{C^{18}O}}\ (2-1)$', fontsize = textsize, color = 'white', verticalalignment = 'top', transform = ax_t.transAxes)

    dlum21_2D_model0_t = dlum21_2D_model0.copy()
    dlum21_2D_model0_t[dlum21_2D_model0_t == 0] = 1e-199
    ratio_t = np.log10(dlum21_2D_model0_t/np.sum(dlum21_2D_model0_t))
    ratio_t[ratio_t < np.min(lim_lum21)] = np.min(lim_lum21)

    c = ax_t.contourf(rr_grid_model0/const.au, zrr_grid_model0, ratio_t.T, lim_lum21, cmap = cmap_selected, extend = 'both', zorder = -1)

    cb = fig.colorbar(c, ax=ax_t, ticks=ticksrange_lum21)
    # cb.set_label(r'$\log_{10}{n_{\mathrm{gas}}}$', labelpad = 15, y = 0.5, rotation=270., fontsize = 12)
    
    # ctau = ax_t.contour(rr_grid_model0/const.au, zrr_grid_model0, np.log10(tauup_CO_2D_model0_t.T), [0, 1], \
    # #                     colors='k', linestyles='solid', zorder = 10)
    # # manual_locations = [(1, 1.1), (5, 0.9), (5, 0.5)] 
    # # plt.clabel(ctau, inline=0, inline_spacing=1, colors = 'k', fontsize=9, fmt = r'$\log_{10}(A_V) = $'+'%1.0f', rightside_up = True, manual = manual_locations)

    # ax_t.text(100, 1.0, r'$\log_{10}(A_V) = $'+'%1.0f' % (-2), zorder = 11)
    # ax_t.text(50, 0.50, '%1.0f' % (0), zorder = 11)
    # ax_t.text(20, 0.25, '%1.0f' % (1), zorder = 11)

    # ax_t.set_xlabel('r [AU]', fontsize = 12)
    # ax_t.set_ylabel('z/r', fontsize = 12)

    ax_t.set_xlim(limxmin, limxmax)
    ax_t.set_ylim(limymin, limymax)

    # ax_t.set_xscale('log')
    # Initialize minor ticks
    ax_t.minorticks_on()


    ###
    # 0, 4 C18O (3-2) w/ CO2
    ###
    if bool_3_2:

        ax_t = axs[4]
        
        if bool_set_title:
            ax_t.set_title(r'$\log_{10} {\mathrm{C^{18}O}}\ (3-2)$', fontsize = textsize)

        # ax_t.text(titlelocator[0], titlelocator[1], r'$\log_{10} {\mathrm{C^{18}O}}\ (3-2)$', fontsize = textsize, color = 'white', verticalalignment = 'top', transform = ax_t.transAxes)

        dlum32_2D_model0_t = dlum32_2D_model0.copy()
        dlum32_2D_model0_t[dlum32_2D_model0_t == 0] = 1e-199
        ratio_t = np.log10(dlum32_2D_model0_t/np.sum(dlum32_2D_model0_t))
        ratio_t[ratio_t < np.min(lim_lum32)] = np.min(lim_lum32)

        c = ax_t.contourf(rr_grid_model0/const.au, zrr_grid_model0, ratio_t.T, lim_lum32, cmap = cmap_selected, extend = 'both', zorder = -1)

        cb = fig.colorbar(c, ax=ax_t, ticks=ticksrange_lum32)
        # cb.set_label(r'$\log_{10}{n_{\mathrm{gas}}}$', labelpad = 15, y = 0.5, rotation=270., fontsize = 12)

        # ctau = ax_t.contour(rr_grid_model0/const.au, zrr_grid_model0, np.log10(tauup_CO_2D_model0_t.T), [0, 1], \
        #                     colors='k', linestyles='solid', zorder = 10)
        # # manual_locations = [(1, 1.1), (5, 0.9), (5, 0.5)] 
        # # plt.clabel(ctau, inline=0, inline_spacing=1, colors = 'k', fontsize=9, fmt = r'$\log_{10}(A_V) = $'+'%1.0f', rightside_up = True, manual = manual_locations)

        # ax_t.text(100, 1.0, r'$\log_{10}(A_V) = $'+'%1.0f' % (-2), zorder = 11)
        # ax_t.text(50, 0.50, '%1.0f' % (0), zorder = 11)
        # ax_t.text(20, 0.25, '%1.0f' % (1), zorder = 11)

        # ax_t.set_xlabel('r [AU]', fontsize = 12)
        # ax_t.set_ylabel('z/r', fontsize = 12)

        ax_t.set_xlim(limxmin, limxmax)
        ax_t.set_ylim(limymin, limymax)

        # ax_t.set_xscale('log')
        # Initialize minor ticks
        ax_t.minorticks_on()
        
    return 1


def plot_m0_map(linecube, ax, vmin=None, vmax=None, dist=None, image_size=8):
    ## calculate the rm using the firs/last 2 channels
    rms =  bm.estimate_RMS(data=linecube.data,N=2)
    ## calculate the moment zero map
    # moments = bm.collapse_zeroth(velax=linecube.velax,data=linecube.data*kep_mask,rms=rms)
    moments = bm.collapse_zeroth(velax=linecube.velax,data=linecube.data,rms=rms)
    M0,dM0 = moments
    
    #####################
    ## plotting the line/continuum
    # ax = axes[i]
    cmap = cmr.rainforest
    norm = None
    
    # vmin,vmax = 0,np.nanmax(M0[linecube.rgrid <= max(3,Rmax)])
    
    cb = ax.imshow(M0,origin='lower',vmin=vmin,vmax=vmax,extent=linecube.extent,cmap=cmap,norm=norm)
    
    xlim,ylim = (image_size/2.,-image_size/2.),(-image_size/2.,image_size/2.)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    
    ax.tick_params(which='both',color='w',width=1.2,size=3) ## change ticks to white to make them visible
    ax.set_xlabel(r'$\Delta$ R.A. [arcsec]',fontsize=10)
    ax.set_ylabel(r'$\Delta$ Decl. [arcsec]',fontsize=10)
    ax.tick_params(axis='both', which='major', labelsize=10)
    
    # Adjust ticks
    ## making custom ticks
    Nt = 2
    half_int = int(image_size/2.)
    half_ticks = np.linspace(0,half_int,Nt+1)
    ticks = [-t for t in reversed(half_ticks[1:])] + [t for t in half_ticks]
    ax.set_xticks([-t for t in ticks])
    ax.set_yticks(ticks)
    
    # add patches and other things
    xmax,xmin = ax.get_xlim()
    ymin,ymax = ax.get_ylim()
    wx,wy = xmax-xmin,ymax-ymin
    
    ## adding the scale bar
    if dist != None:
        scale_bar = 100/dist
        anchor = (xmin+wx/8,ymax-wy/16)
        ax.errorbar([anchor[0]+scale_bar,anchor[0]],[anchor[1]]*2,lw=1.8,color='w',yerr=0.2)
        ax.annotate('100 au',(anchor[0]+scale_bar/2,anchor[1]-0.3),xycoords='data',
                    ha='center',va='top',fontsize=8,color='w')
    ## adding the beam
    Bmaj,Bmin,Bpa = linecube.beam
    print('Beam: ', Bmin, Bmaj, Bpa)

    radius = Bmaj*0.5
    R_frac = radius/wx
    pad_frac = 0.08
    xy = (xmax - (pad_frac+R_frac)*wx,(pad_frac+R_frac)*wy+ymin)
    
    ax.add_patch(patches.Ellipse(xy,width=Bmin,height=Bmaj,angle=180-Bpa,edgecolor='w',facecolor='None',lw=2,hatch='///////'))
    
    return cb


def sph2cyl(rc, thetac, phic, rhogas, Tgas, nr_cyl=500, nz_cyl=500, nr_dense = 500, ntheta_dense = 500, method_polate1 = 'linear', method_polate2 = 'nearest'):
    """
    # in this def, 1st the values will be extended in the spherical coordiante 
    # the extended spherical coordinate will be both extending outwards and sampling
    # from a COARSE to a DENSE coordinate
    # and then they will be conerted into cylindrical coordiante
    """
    
    rhogas_sph_ext, Tgas_sph_ext, qq_ext = extend_insphr(rc, thetac, phic, rhogas, Tgas)
    
    #
    # Make the meshgrid
    # 
    # For original grid
    qq     = np.meshgrid(rc,thetac,phic,indexing='ij')
    rr     = qq[0] # Spherical R 
    tt     = qq[1] # Angle with z axis ( = pi/2 for the z=0 axis)
    zr     = np.pi/2.e0 - qq[1] 
    zz     = np.sin(zr)*rr
    nphi   = len(phic)
    # For extended grid
    rr_ext     = qq_ext[0] # Spherical R 
    tt_ext     = qq_ext[1] # Angle with z axis ( = pi/2 for the z=0 axis)
    zr_ext     = np.pi/2.e0 - qq_ext[1] 
    zz_ext     = np.sin(zr_ext)*rr_ext
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
    # Make Meshgrid
    rc_dense = np.logspace(np.log10(np.min(rr_ext)), np.log10(np.max(rr_ext)), nr_dense)
    thetac_dense = thetac.copy()
    while len(thetac_dense) <= ntheta_dense:
        thetac_dense = np.hstack([thetac_dense, 0.5*(thetac_dense[:-1]+thetac_dense[1:])])
    thetac_dense = np.sort(thetac_dense)
    
    qq_dense = np.meshgrid(rc_dense,thetac_dense,phic,indexing='ij')
    rr_dense = qq_dense[0] # Spherical R 
    tt_dense = qq_dense[1] # Angle with z axis ( = pi/2 for the z=0 axis)
    # zr_dense = np.pi/2.e0 - qq_dense[1]
    
    nr_dense = len(rc_dense)
    ntheta_dense = len(thetac_dense)
    
    print('the DENSE grid: ', nr_dense, ntheta_dense)
    
    #
    # linear interpolate to setup all the cells in spherical coordinate
    #
    Tgas_sph_dense_old = np.zeros([nr_dense,ntheta_dense,nphi])
    grid_dense_Tg = interpolate.griddata((np.log10(rr_flat_ext_old/au), tt_flat_ext_old), np.log10(Tgas_flat_ext_old), (np.log10(rr_dense/au), tt_dense), method=method_polate1)
    Tgas_sph_dense_old = 10**grid_dense_Tg.copy()
    
    rhogas_sph_dense_old = np.zeros([nr_dense,ntheta_dense,nphi])
    grid_sph_rg = interpolate.griddata((np.log10(rr_flat_ext_old/au), tt_flat_ext_old), np.log10(rhogas_flat_ext_old), (np.log10(rr_dense/au), tt_dense), method=method_polate1)
    rhogas_sph_dense_old = 10**grid_sph_rg.copy()

    #
    # Name all the data in extended DENSE grid
    #
    r_flat_cyl_old = rr_dense[:, :, 0].flatten() * np.sin(tt_dense[:, :, 0].flatten())
    z_flat_cyl_old = rr_dense[:, :, 0].flatten() * np.cos(tt_dense[:, :, 0].flatten())
    zr_flat_cyl_old = z_flat_cyl_old/r_flat_cyl_old
    Tgas_flat_sph_dense_old = Tgas_sph_dense_old[:, :].flatten()
    rhogas_flat_sph_dense_old = rhogas_sph_dense_old[:, :].flatten()
    
    #
    # New Grid for Sampling in the scipy griddata regridding
    # The new grid will exactly in r and z/r grids
    # 
    rc_new = np.logspace(np.log10(0.95*np.min(rc)), np.log10(1.05*np.max(rc)), nr_cyl)
    zrc_new = np.linspace(np.min(1/np.tan(tt)), np.max(1/np.tan(tt)), nz_cyl)
    zrc_new = zrc_new[::-1]
    # zrc_new = 1/np.tan(tt)
    qq_new  = np.meshgrid(rc_new, zrc_new, phic, indexing = 'ij')
    rq_new = qq_new[0]; zrq_new = qq_new[1]

    #
    # Regridding to a New Grid
    # The Method is 'cubic' because it is an Oversampling
    # But the 'cubic' method will bring nan value and further destroy the results
    # so now adopt 'linear'
    #
    #     method_polate = 'cubic' # now it is defined in the input
    
    Tgas_cyl_old = np.zeros([nr_cyl,nz_cyl,nphi])
    grid_rz_Tg = interpolate.griddata((np.log10(r_flat_cyl_old/au), zr_flat_cyl_old), np.log10(Tgas_flat_sph_dense_old), (np.log10(rq_new/au), zrq_new), method=method_polate2)
    Tgas_cyl_old = 10**grid_rz_Tg.copy()
    
    rhogas_cyl_old = np.zeros([nr_cyl,nz_cyl,nphi])
    grid_rz_rg = interpolate.griddata((np.log10(r_flat_cyl_old/au), zr_flat_cyl_old), np.log10(rhogas_flat_sph_dense_old), (np.log10(rq_new/au), zrq_new), method=method_polate2)
    rhogas_cyl_old = 10**grid_rz_rg.copy()
    
    return rc_new, zrc_new, qq_new, rhogas_cyl_old, Tgas_cyl_old


def cyl2sph(rc, thetac, phic, rr_cyl, zrr_cyl, rhogas_cyl_new, Tgas_cyl_old, nr_dense = 500, ntheta_dense = 500, method_polate1 = 'nearest', method_polate2 = 'linear'):
    """
    # Convert From Cylindrical to Sphrical Coordiante
    """
    
    nr,ntheta,nphi = len(rc), len(thetac), len(phic)
    
    R_flat_cyl = np.sqrt(rr_cyl.flatten()**2 + (zrr_cyl.flatten()*rr_cyl.flatten())**2)
    tt_flat_cyl = np.pi/2 - np.arctan(zrr_cyl.flatten())

    rhogas_flat_cyl_new = rhogas_cyl_new.flatten()
    Tgas_flat_cyl = Tgas_cyl_old.flatten()

    #
    # DENSE grid for first sampling back
    #
    # Make Meshgrid
    rc_dense = np.logspace(np.log10(np.min(rc)), np.log10(np.max(rc)), nr_dense)
    thetac_dense = thetac.copy()
    while len(thetac_dense) <= ntheta_dense:
        thetac_dense = np.hstack([thetac_dense, 0.5*(thetac_dense[:-1]+thetac_dense[1:])])
    thetac_dense = np.sort(thetac_dense)
    
    qq_dense = np.meshgrid(rc_dense,thetac_dense,phic,indexing='ij')
    rr_dense = qq_dense[0] # Spherical R 
    tt_dense = qq_dense[1] # Angle with z axis ( = pi/2 for the z=0 axis)
    # zr_dense = np.pi/2.e0 - qq_dense[1]
    
    nr_dense = len(rc_dense)
    ntheta_dense = len(thetac_dense)
    
    # sampling to the DENSE spherical coordinate
    Tgas_sph_dense_new = np.zeros([nr_dense,ntheta_dense,nphi])
    grid_dense_Tg = interpolate.griddata((np.log10(R_flat_cyl/au), tt_flat_cyl), np.log10(Tgas_flat_cyl), (np.log10(rr_dense/au), tt_dense), method=method_polate1)
    Tgas_sph_dense_new = 10**grid_dense_Tg.copy()
    
    rhogas_sph_dense_new = np.zeros([nr_dense,ntheta_dense,nphi])
    grid_sph_rg = interpolate.griddata((np.log10(R_flat_cyl/au), tt_flat_cyl), np.log10(rhogas_flat_cyl_new), (np.log10(rr_dense/au), tt_dense), method=method_polate1)
    rhogas_sph_dense_new = 10**grid_sph_rg.copy()
    
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
    qq     = np.meshgrid(rc,thetac,phic,indexing='ij')
    rr     = qq[0] # Spherical R 
    tt     = qq[1] # Angle with z axis ( = pi/2 for the z=0 axis)
    zr     = np.pi/2.e0 - qq[1] 
    zz     = np.sin(zr)*rr

    #
    # Regridding to a New Grid
    # Adopt method of 'linear' because of undersampling
    # 

    Tg_sphr_new = np.zeros([nr,ntheta,nphi])
    grid_rz_Tg = interpolate.griddata((np.log10(rr_flat_dense_new/au), tt_flat_dense_new), np.log10(Tgas_flat_dense_new), (np.log10(rr/au), tt), method=method_polate2)
    Tg_sphr_new = 10**grid_rz_Tg.copy()

    rhogas_sphr_new = np.zeros([nr,ntheta,nphi])
    grid_rz_rg = interpolate.griddata((np.log10(rr_flat_dense_new/au), tt_flat_dense_new), np.log10(rhogas_flat_dense_new), (np.log10(rr/au), tt), method=method_polate2)
    rhogas_sphr_new = 10**grid_rz_rg.copy()
    rhogas_sphr_new[np.isnan(rhogas_sphr_new)] = 0

    return rhogas_sphr_new, Tg_sphr_new

# os.chdir(data_dir)

def read_mint_file(model_name, file_dir_readin, parameter_file_dir):
    """read in the mint model parameter files"""
    data_dir = file_dir_readin # os.path.join(work_dir, 'data')
    para = model.Parameters()

    para.read_parameters_from_csv(directory=parameter_file_dir, filename=model_name+'_parameters', extension='.csv')

    ratio_g2d_reference_file = os.path.join(file_dir_readin, 'ratio_g2d_reference.dat')
    if os.path.exists(ratio_g2d_reference_file):
        para.ratio_g2d = np.loadtxt(ratio_g2d_reference_file)

    mint = model.Mint(para, file_dir=data_dir)
    mint.setup_dust_info(para, bool_savefile=False)
    
    return mint, para


# def read_radmc3d_density(file_dir_readin, mint):
    """read in the radmc3d density profile"""
    current_dir = os.getcwd()
    os.chdir(file_dir_readin)
    print('reading radmc3d format files from')
    print(os.getcwd())

    dustdens = radmc3dPy.analyze.readData(dtemp=True, ddens=True, gtemp=True) # also read in gastemperature

    rhodust = dustdens.rhodust.copy()
    rhodust_all = rhodust.sum(3).copy()

    ratio_g2d_grid_file = 'ratio_g2d_grid.dat'
    # ratio_g2d_grid_file = 'ratio_g2d_grid.new'
    # ratio_g2d_grid_file = 'ratio_g2d_VHSE_grid_1.dat'
    if os.path.exists(ratio_g2d_grid_file):
        print('read in the reference ratio_g2d_grid from %s'%(ratio_g2d_grid_file))
        ratio_g2d_grid_readin = dd.load_3dgrid_from_dat(ratio_g2d_grid_file)
        ratio_g2d_grid = ratio_g2d_grid_readin
        ratio_g2d_grid[ratio_g2d_grid == np.inf] = np.nan
    #     ratio_g2d_grid[np.isnan(ratio_g2d_grid)] = 0.0

        mint.ratio_g2d = ratio_g2d_grid

    rhogas = rhodust_all * mint.ratio_g2d
    rhogas[np.isnan(rhogas)] = 1e-99

    # 
    # Read dust temperature
    #
    Tdust = dustdens.dusttemp.copy()
    #
    # Read gas temperature
    #
    Tgas_0 = dustdens.gastemp.copy()

    #
    # Fill in those 0 blocks
    #
    rhogas[rhogas == 0.] = 1e-99
    Tgas_0[Tgas_0 == 0.] = 1e-99

    # test the total mass to see whether the mass is conservative to what it should be set upackage_position

    mdust_set = (rhodust_all * mint.vol).sum(0).sum(0).sum(0) * 2
    # mgas_set = (rhodust_all * mint.vol * mint.ratio_g2d).sum(0).sum(0).sum(0) * 2
    mgas_set = (rhogas * mint.vol).sum(0).sum(0).sum(0) * 2

    print('the read in radmc3d mdust = %.2e [ms]'%(mdust_set/const.ms))
    print('the read in radmc3d mgas = %.2e [ms]'%(mgas_set/const.ms))
    print('the set up gas-to-dust mass ratio = %.2f'%(mgas_set/mdust_set))

    # plt.plot(ratio_g2d_grid_readin[:, -1, 0])
    
    os.chdir(current_dir)
    print('changing dir back to %s'%(current_dir))
    
    return dustdens, rhodust, rhogas, Tdust, Tgas_0


def get_sigma_interval(rho_t, dustdens):
    nr     = dustdens.grid.nx
    rc     = dustdens.grid.x
    thetac = dustdens.grid.y
    phic   = dustdens.grid.z
    qq     = np.meshgrid(rc,thetac,phic,indexing='ij')
    rr     = qq[0] # Spherical R 
    zr     = np.pi/2.e0 - qq[1]
    zz     = np.sin(zr)*rr

    # ri = dustdens.grid.xi.copy()
    integral_rho_t = np.array([np.trapz(rho_t[i_r, :, 0], -zz[i_r, :, 0]) for i_r in range(nr)])
    # dmass_t = 2 * 2 * np.pi * rc * (ri[1:] - ri[:-1]) * integral_rho_t
    dsigma_t = 2 * integral_rho_t
    
    return dsigma_t


def image_convolve(filein, fwhm = [0.32, 0.32]):
    """
    convolve the image with the 2D Guassian core
    
    Parameters:
        filein (str): Base name of the FITS file (without '.fits' extension).
        fwhm (list): Full width at half maximum in arcsec for the Gaussian kernel [fwhm_x, fwhm_y].

    Returns:
        int: 1 if the process completes successfully.
        
    Output: 
        write out the file with the Base name of the FITS file and extension of '_conv.fits'.
    """
    
    # Open FITS file
    if '.fits' in filein:
        filein = filein.replace('.fits', '')
        
    filename = filein + '.fits'
    if not os.path.exists(filename):
        raise FileNotFoundError(f"{filename} not found.")
    print(f"[INFO] Reading in original image: {filename}")
    hdul = fits.open(filename)

    # Get the RA, DEC and channel values from the header
    # The axis 1 is the RA, the axis 2 is the DEC and the
    # axis 3 is the velocity channels

    header = hdul[0].header    # Primary HDU header
    naxis1 = header['naxis1']  # Number of values in each axis
    naxis2 = header['naxis2']  # naxis1 and naxis2 should be equal
    naxis3 = header['naxis3']

    cdelt1 = header['cdelt1']  # Difference between two pixel
    cdelt2 = header['cdelt2']  # cdelt1 and cdelt2 should be opposite
    cdelt3 = header['cdelt3']  # cdelt1 negative

    crpix1 = header['crpix1']  # Reference pixel
    crpix2 = header['crpix2'] 
    crpix3 = header['crpix3']

    crval1 = header['crval1']  # Value at the reference pixel
    crval2 = header['crval2']
    crval3 = header['crval3']

    # Compute the values of RA and DEC

    ra = (np.arange(naxis1) + 1 - crpix1) * cdelt1 + crval1
    dec = (np.arange(naxis2) + 1 - crpix2) * cdelt2 + crval2
    veloc = (np.arange(naxis3) + 1 - crpix3) * cdelt3 + crval3

    # The values of RA and DEC are in degrees in the fits file
    # Velocities are in m/s. We convert them in arcsec and km/s

    ra = ra * 3600 # arcsec 
    dec = dec * 3600 # arcsec
    veloc = veloc / 1e3 # km/s

    # Fix keyword for velocity
    # header['ctype3'] = 'VELO-LSR'

    # Now read the brightness temperatures
    tb = hdul[0].data

    # Construct Gaussian kernel
    resolution = abs(ra[1] - ra[0])  # arcsec/pixel
    stdev = [fw / (2 * np.sqrt(2 * np.log(2))) for fw in fwhm]  # convert FWHM to stddev
    beam = Gaussian2DKernel(x_stddev=stdev[0] / resolution, y_stddev=stdev[1] / resolution)

    # Convolve the brightness temperature by the beam pattern
    # in each channel

    dim_ra = len(ra)
    dim_dec = len(dec)
    dim_veloc = len(veloc)

    ta = np.zeros(dim_ra * dim_dec * dim_veloc, dtype = np.float64)
    ta = np.reshape(ta, (1, dim_veloc, dim_ra, dim_dec))

    # Convolve each velocity channel
    # convolved_data = np.zeros_like(tb) 
    convolved_data = np.zeros(dim_ra * dim_dec * dim_veloc, dtype = np.float64)
    convolved_data = np.reshape(convolved_data, (1, dim_veloc, dim_ra, dim_dec))
    
    print(f"[INFO] Starting imaging convlution channel by channel")
    # for k in range(naxis3):
    #     convolved_data[0, k, :, :] = convolve_fft(tb[0, k, :, :], beam, allow_huge=True)
    try:
        from tqdm import tqdm
        loop = tqdm(range(naxis3), desc="Convolving")
    except ImportError:
        print(f"[WARNNING] the Python package tqdm is not found. Recommand using tqdm to show the progress bar.")
        loop = range(naxis3)

    if len(tb.shape) == 3:
        for k in loop:
            convolved_data[0, k, :, :] = convolve_fft(tb[k, :, :], beam, allow_huge=True)
    else:
        for k in loop:
            convolved_data[0, k, :, :] = convolve_fft(tb[0, k, :, :], beam, allow_huge=True)
    
    # Update header
    header_out = copy.deepcopy(header)
    header_out['BMAJ'] = fwhm[0] / 3600.  # arcsec -> deg
    header_out['BMIN'] = fwhm[1] / 3600.
    header_out['BPA']  = 0.0
    header_out['BUNIT'] = 'JY/BEAM'

    # Apply Jy/pix -> Jy/beam conversion
    pixel_area = abs(cdelt1 * 3600 * cdelt2 * 3600)  # arcsec^2
    beam_area = fwhm[0] * fwhm[1] * np.pi / (4.0 * np.log(2.0))  # arcsec^2
    convolved_data *= (beam_area / pixel_area)

    # Save output FITS file
    fileout = filein + '_conv.fits'
    if os.path.exists(fileout):
        os.remove(fileout)
    fits.PrimaryHDU(convolved_data, header=header_out).writeto(fileout)
    print(f"[INFO] Convolved file saved to: {fileout}")

    hdul.close()
    return 1


def luminosity_from_flux_density(
    flux_mJy,
    d_pc,
    wavelength=None,
    frequency=None,
    return_in="cgs"  # "cgs" -> erg/s, "SI" -> W, "Lsun" -> L_sun
):
    """
    Convert observed flux density (mJy) at a given band to monochromatic luminosity.

    Parameters
    ----------
    flux_mJy : array-like or float
        Flux density values in mJy (per Hz).
    d_pc : float
        Distance in parsecs.
    wavelength : Quantity or float, optional
        Band central wavelength (e.g., 2.0*u.um). Provide exactly ONE of wavelength or frequency.
    frequency : Quantity or float, optional
        Band central frequency (e.g., 150*u.GHz). Provide exactly ONE of wavelength or frequency.
    return_in : {"cgs", "SI", "Lsun"}
        Unit for νLν output. Lν is always returned in both cgs and SI.

    Returns
    -------
    Lnu_SI : Quantity
        Monochromatic luminosity Lν in W/Hz (same shape as input).
    Lnu_cgs : Quantity
        Monochromatic luminosity Lν in erg/s/Hz.
    nuLnu : Quantity
        νLν in requested units (W, erg/s, or L_sun).
    """
    # Validate inputs
    if (wavelength is None) == (frequency is None):
        raise ValueError("Provide exactly one of 'wavelength' or 'frequency'.")

    flux_mJy = np.asanyarray(flux_mJy)
    S_nu = (flux_mJy * 1e-3) * u.Jy                   # mJy -> Jy
    d = (d_pc * u.pc).to(u.m)                          # parsec -> meters

    if wavelength is not None:
        wl = wavelength*u.m if not isinstance(wavelength, u.Quantity) else wavelength.to(u.m)
        nu = (C.c / wl).to(u.Hz)
    else:
        nu = frequency*u.Hz if not isinstance(frequency, u.Quantity) else frequency.to(u.Hz)

    # Lν = 4π d^2 Sν
    Lnu_SI = (4 * np.pi * d**2 * S_nu.to(u.W / u.m**2 / u.Hz)).to(u.W / u.Hz)

    # νLν (common continuum luminosity proxy)
    nuLnu_W = (nu * Lnu_SI).to(u.W)

    if return_in.lower() == "cgs":
        Lnu = Lnu_SI.to(u.erg / u.s / u.Hz)
        nuLnu = nuLnu_W.to(u.erg / u.s)
    elif return_in.lower() == "si":
        Lnu = Lnu_SI
        nuLnu = nuLnu_W
    elif return_in.lower() == "lsun":
        Lnu = Lnu_SI.to(u.Lsun / u.Hz)
        nuLnu = nuLnu_W.to(u.Lsun)
    else:
        raise ValueError("return_in must be one of: 'cgs', 'SI', 'Lsun'.")

    return Lnu, nuLnu
