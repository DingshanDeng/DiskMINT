# -*- coding: utf-8 -*-
# author: Dingshan Deng @ University of Arizona
# contact: dingshandeng@arizona.edu
# created: 07/03/2023

"""
Wrapper of using LIME to make a group of synthetic images
This is a very simple example of using LIME v1.9.5
and this script is prepared in `bash` (also work in `zsh`) shell
Please prepare your own script

To modify this script for your own purpose,
Please change both the parameters here in the section
`Changing Parameters`
and also the `colimemodel.c` accordingly in the same directory
"""

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.convolution import convolve_fft, Gaussian2DKernel

import os, copy
# tqdm package is recommanded to see the progress
# because sometime image convolution takes long
import tqdm

#
# Some natural and astronomy constants in the cgs unit
#
au  = 1.49598e13     # Astronomical Unit       [cm]

def image_convolve(LIME_data_dir, file2_name_ind, lineJ = '2-1', fwhm = [0.32, 0.32]):
    
    file2_mod2 = LIME_data_dir + file2_name_ind + '.fits'

    # fwhm = [0.32, 0.32]

    filein = 'LIME_image_c18o_'+ lineJ # 'LIME_image_c18o_3-2' #
    hdul = fits.open(filein + '.fits')

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

    crval1 = header['crval1']  # Value at the reference pixel
    crval2 = header['crval2']
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

    ra = ra * 3600 # arcsec 
    dec = dec * 3600 # arcsec
    veloc = veloc / 1e3 # km/s

    # Fix keyword for velocity
    # header['ctype3'] = 'VELO-LSR'

    # Now read the brightness temperatures
    tb = hdul[0].data

    resol = ra[0] - ra[1]
    stdev = fwhm / (2 * np.sqrt(2 * np.log(2)))

    beam = Gaussian2DKernel(x_stddev = stdev[0] / resol, y_stddev = stdev[1] / resol)

    # Convolve the brightness temperature by the beam pattern
    # in each channel

    dim_ra = len(ra)
    dim_dec = len(dec)
    dim_veloc = len(veloc)

    ta = np.zeros(dim_ra * dim_dec * dim_veloc, dtype = np.float64)
    ta = np.reshape(ta, (1, dim_veloc, dim_ra, dim_dec))

    # Compute the convolution using the fft
    for k in tqdm.tqdm(range(dim_veloc)):
        ta[0,k,:,:] = convolve_fft(tb[0,k,:,:], beam, allow_huge=True)

    header_o = copy.deepcopy(header)
    header_o['BMAJ'] = fwhm[0] / 3600.
    header_o['BMIN'] = fwhm[1] / 3600.
    header_o['BPA']  = 0.0

    # Convert from Jy/pixel to Jy/beam
    pixel_area = (-cdelt1*3600 * cdelt2*3600)
    beam_area = fwhm[0] * fwhm[1] * np.pi / 4. /np.log(2.0)
    ta_o = ta * beam_area/pixel_area

    header_o['BUNIT'] = 'JY/BEAM'

    hdul_o = fits.PrimaryHDU(ta_o, header = header_o)

    file2_mod2_out = LIME_data_dir + file2_name_ind + '_conv.fits'

    print(file2_mod2_out)
    os.system('rm ' + file2_mod2_out)
    # save the convolved file to data_dir
    hdul_o.writeto(file2_mod2_out)
    # also move the original file to data_dir
    os.system('mv ./' + filein + '.fits ' + LIME_data_dir + file2_name_ind + '.fits')
    
    
def write_LIME(chem_save_name, data_dir_init, data_dir_chemical, chem_code_name='reducedRGH22', nr_LIME = 100, nz_LIME = 200):
    #
    # Read in the CO init Grid
    #
    cord_t = 'cyl' 
    df_COinit_1 = pd.read_csv(data_dir_init + 'COinitgrid-' + chem_save_name + '-' + cord_t + '.dat', skiprows=3)

    Tgas_ar = np.array(df_COinit_1['Tgas[K]'].values, dtype = np.float64)

    #
    # Read in the CO end Grid
    #
    xc18o_1_chemical = np.loadtxt(data_dir_chemical + 'COendgrid-' + chem_code_name + '-' + chem_save_name + '-' + cord_t + '.chem') 

    #
    # Combine and save the LIME grid
    #
    df_COinit_comb = df_COinit_1.copy()
    df_COinit_comb['r-chem[cm]'] = xc18o_1_chemical[:, 0]
    df_COinit_comb['z-chem[cm]'] = xc18o_1_chemical[:, 1]
    df_COinit_comb['X(H2)']      = 10**xc18o_1_chemical[:, 2]
    df_COinit_comb['X(C18O)']    = 10**xc18o_1_chemical[:, 3]
    
    #
    # Read in rr and zz
    #
    r_init_ar  = df_COinit_1['r[cm]'].values
    zr_init_ar = df_COinit_1['z/r'].values

    rr_2D  = np.reshape(r_init_ar, (nr_LIME, nz_LIME))
    zrr_2D = np.reshape(zr_init_ar, (nr_LIME, nz_LIME))

    rr_grid = rr_2D[:,0]
    zrr_grid = zrr_2D[0,:]
    
    rr_cyl = rr_grid.copy()
    zr_cyl = zrr_grid.copy()
    zr_chosen = zr_cyl.copy()
    print(zr_chosen.shape)

    row_numb = nr_LIME * nz_LIME
    col_numb = 11

    X_unk  = 9.990E-01
    X_e    = 0.0
    orthoH = 0.75
    paraH  = 0.25

    with open('limegrid.dat', 'w+') as f:
        f.write('        ' + str(int(row_numb)) + '        11 \n')

        for i_row in range(len(df_COinit_comb)):

            rinau  = df_COinit_comb.loc[i_row, 'r[cm]']/au
            zinau  = df_COinit_comb.loc[i_row, 'z[cm]']/au
            n_H    = df_COinit_comb.loc[i_row, 'n_H[cm**-3]'] * 1e6
            T_gas  = np.float64(df_COinit_comb.loc[i_row, 'Tgas[K]'])
            
            if 'Tdust[K]' in df_COinit_comb.columns:
                T_dust = np.float64(df_COinit_comb.loc[i_row, 'Tdust[K]'])
            else:
                T_dust = T_gas
                
            X_H2   = df_COinit_comb.loc[i_row, 'X(H2)']
            # X_H2   = 0.5 # this is assuming all H2 is preserved
            X_C18O = df_COinit_comb.loc[i_row, 'X(C18O)']

            this_line = \
            '  {:.3E}'.format(rinau) +\
            '  {:.3E}'.format(zinau) +\
            '  {:.3E}'.format(n_H)   +\
            '  {:.3E}'.format(T_gas) +\
            '  {:.3E}'.format(T_dust)+\
            '  {:.3E}'.format(X_H2)  +\
            '  {:.3E}'.format(X_unk) +\
            '  {:.3E}'.format(X_e)   +\
            '  {:.3E}'.format(X_C18O)+\
            '  {:.3E}'.format(orthoH)+\
            '  {:.3E}'.format(paraH) +\
            '\n'

            f.write(this_line)

    r_init_ar  = df_COinit_1['r[cm]'].values
    zr_init_ar = df_COinit_1['z/r'].values

    rr_2D    = np.reshape(r_init_ar, (nr_LIME, nz_LIME))
    zrr_2D   = np.reshape(zr_init_ar, (nr_LIME, nz_LIME))

    rr_grid  = rr_2D[:,0]
    zrr_grid = zrr_2D[0,:]

    rr_cyl = rr_grid.copy()
    zr_cyl = zrr_grid.copy()

    row_numb = len(rr_cyl)

    with open('rgrid.dat', 'w+') as f:
        f.write('        ' + str(int(row_numb)) + '\n')

        for i_r in range(len(rr_cyl)):

            rinau = rr_cyl[i_r] / au
            rlbg  = 0 + i_r * len(zr_cyl)
            rled  = len(zr_cyl) - 1 + i_r * len(zr_cyl)

            this_line = \
            '  {:.3E}'.format(rinau) +\
            '  {:.0f}'.format(rlbg) +\
            '  {:.0f}'.format(rled) +\
            '\n'

            f.write(this_line)


    chem_save_name_save = chem_save_name + '_' + chem_code_name
    file_name = chem_save_name_save

    if not os.path.isdir('./' + file_name):
        os.mkdir('./' + file_name)

    output_dir2 = './' + file_name

    os.system('cp limegrid.dat ' + output_dir2)
    os.system('cp rgrid.dat ' + output_dir2)
    
    print('Saved!', chem_save_name_save)
    
    return

"""
Changing Parameters
Caution! Change the LIME parameters in colimemodel.c accordingly!
"""

#
# Where is your LIME bin files
#
LIMEfile_dir = 'Yourpath/lime-1.9.5/lime'

#
# Which line, for example it can be '2-1', or '3-2' 
#
lineJ = '2-1'

#
# Where to save the LIME created Images
# Images will be .fits file and they could be HUGE
# (can be up to a few GB)
# Please make sure you have enough free space
#
LIME_data_dir = './images/RULup_Image_C18O/'

#
# The model name (the name you used for the model)
#
chem_save_name_list = [
    'example_model_RULup_VHSE_bestfit',
]

#
# The dust opacity file
# Because LIME does not support multiple dust species (multiple sizes)
# So here this need to be the same dust composition used for SED
# But IT NEEDS TO BE COMBINED TO ONE file instead of multiple files for different sizes
#
dustkappa_list = [
    'dustkappa_WD01comb.inp',
]

#
# Global ratio of the gas-to-dust mass
# We note that because LIME is not doing the dust continuum radiative transfer properly
# because LIME assumes the density of dust based on the gas collider density with a constant g2d
# that we set up here, while sometimes the real dust density is slightly different.
# so that the continuum emission might be off, and you might need to change the ratio_g2d here to
# perfectly match the dust continuum emission in the image.
#
ratio_g2d_list = np.array([30.0])

#
# The Post Tag Name for LIME generated images
# You might want to change the LIME set ups and generate different images
# such as nonLTE or LTE ones or whatever you want.
#
save_name_post_fix = 'incl19_LTE'

#
# Name for the Chemical Code
#
chem_code_name = 'reducedRGH22'

#
# dir for reading CO init grid
#
data_dir_init     = '../data/'

#
# dir for reading CO end grid
#
data_dir_chemical = '../data/'

#
# The synthetic image FWHM
#
fwhm_set = [0.32, 0.32] # arcsec

"""
Main Code
"""

i_mod = -1
for chem_save_name in chem_save_name_list:
        
    i_mod = i_mod+1
    ratio_g2d_t = ratio_g2d_list[i_mod]
    print('ratio_g2d of %.2f start running'%(ratio_g2d_t))
    np.savetxt('ratio_g2d.dat', [ratio_g2d_t], fmt='%3.2f')
    
    file2_name_ind = 'LIME_C18O' + lineJ + '_' + chem_save_name  + '_' + chem_code_name + '_' + save_name_post_fix
    
    os.system('rm limegrid.dat')
    os.system('rm rgrid.dat')
    
    print('making new files and will be saved to ', chem_save_name  + '_' + chem_code_name)
    # in the example of using LIME, the grid for LIME is hard coded to 100 x 200
    write_LIME(chem_save_name, data_dir_init, data_dir_chemical, chem_code_name, nr_LIME = 100, nz_LIME = 200)
    
    """
    Caution! Change the LIME parameters in colimemodel.c accordingly!
    """
    os.system('cp '+dustkappa_list[i_mod]+' dustkappa_currentlyusing.inp')
    # os.system(LIMEfile_dir + ' -s colimemodel.c') # -s used to Suppresses Output
    os.system(LIMEfile_dir + ' -n colimemodel.c') # -n for normal output
    image_convolve(LIME_data_dir, file2_name_ind, lineJ = lineJ, fwhm = fwhm_set)
    
