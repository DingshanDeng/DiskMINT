#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author: Dingshan Deng @ University of Arizona
# contact: dingshandeng@arizona.edu
# created: 05/18/2023
# new version updated: 04/08/2025
# DiskMINT v1.5.0

import copy, os
import traceback

try:
    import numpy as np
except ImportError:
    np = None
    print(' Numpy cannot be imported ')
    print(' To use this python module,  you need to install Numpy')
    print(traceback.format_exc())

from .constants import *
from . import disk_density as dd

class Parameters:
    """
    diskmint.parameters is the object for DiskMINT model
    to keep all the parameters that are needed.
    
    Here we have two modes to set up the parameters
    
    1. directly set up the parameters in Python code as changing the attributes
    2. read in a csv file to set up the parameters.
    """

    def __init__(self):
        """Initialize an instance."""
        self.parameters = {}
        self.set_default_parameters()

    def set_parameter(self, name, value, units=None, description=None, valuetype='float64'):
        """
        Set the value of a parameter.

        Args:
            name (str): The name of the parameter.
            value (float): The value of the parameter.
            units (str): The units of the parameter.
            valuetype (str): The value type of the parameter, unless specified as 
                It can be set as 'str', 'list', 'int', 'float'
                default (and any other specified name) to be 'float64'. 
                Currently 'int' is also stroed as 'float64'.
            description (str): The description of the parameter.
        """
        self.parameters[name] = {
            'value': value,
            'units': units,
            'valuetype': valuetype,
            'description': description
        }

    def edit_parameter(self, name, new_value=None, new_units=None, new_description=None, new_valuetype=None):
        """
        Edit the value, units, or description of a parameter.

        Args:
            name (str): The name of the parameter to edit.
            new_value (float, optional): The new value for the parameter.
            new_units (str, optional): The new units for the parameter.
            new_description (str, optional): The new description for the parameter.
        """
        if name in self.parameters:
            if new_value is not None:
                self.parameters[name]['value'] = new_value
            if new_units is not None:
                self.parameters[name]['units'] = new_units
            if new_description is not None:
                self.parameters[name]['description'] = new_description
            if new_valuetype is not None:
                self.parameters[name]['valuetype'] = new_valuetype
        else:
            if new_value is not None and new_valuetype is None:
                self.set_parameter(name, new_value, new_units, new_description)
            elif new_value is not None and new_valuetype is not None:
                self.set_parameter(name, new_value, new_units, new_description, new_valuetype)
            else:
                print('Warnning! But keep going.')
                print('new parameters are set but without any values!')
                print('the new parameters is %s'%(name))
            
    def get_parameter_info(self, name):
        """
        Get the information of a parameter.

        Args:
            name (str): The name of the parameter.

        Returns:
            dict: A dictionary containing the value, units, and description of the parameter.
        """
        if name in self.parameters:
            return self.parameters[name]

    def get_parameter_value(self, name):
        """
        Get the value of a parameter.

        Args:
            name (str): The name of the parameter.

        Returns:
            float: the value of the parameter
        """
        if name in self.parameters:
            return self.parameters[name]['value']
         

    def read_parameters_from_csv(self, filename, directory='./', extension='', verbose=False):
        """
        return self.parameters[name]['value']
        Read parameter information from a CSV file.

        Args:
            directory (str): The directory where to save this CSV file.
            filename (str): The name of the CSV file.
            extension (str): The file extension, recommand using '.csv'.
        """
        # the parameters that will be a grid so it would not be in the main dat file
        # and they will be stored in stand alone files
        stand_alone_parameter_list = ['ratio_g2d', 'sigmad_ref']
        
        filename_full = os.path.join(directory, filename+extension)
        
        with open(filename_full, 'r') as file:
            for line in file:
                line = line.strip()

                if not line.startswith('#') and line:
                    # skip the lines that start with '#' or are empty
                    
                    data = line.split(',')
                    if verbose:
                        print(line)
                        print(data)
                    
                    name = data[0].strip()
                    value = data[1].strip()
                    units = data[2].strip()
                    valuetype = data[3].strip()
                    description = data[4].strip()
                    # print(value)
                    
                    if name not in stand_alone_parameter_list:
                        # then the value would be the value here
                        if valuetype in ['str', 'list']:
                            value = value
                        elif valuetype == 'int':
                            value = int(value)
                        else:
                            try:
                                value = np.float64(value)
                            except ValueError:
                                print('the input value is not str and cannot be converted to float64')
                                print('If it is a str or list, please identify in the valuetype column')
                                value = value
                        
                    else:
                        # if name in stand_alone_parameter_list
                        # the value would be stroed in a another file in the same directory
                        stand_alone_file_name = filename.replace('.', '_') + '_' + name + '.dat'
                        value = np.loadtxt(os.path.join(directory, stand_alone_file_name))
                        
                    self.set_parameter(name, value, units, description, valuetype=valuetype)
                        
        return 1


    def write_parameters_to_csv(self, filename, directory='./', extension=''):
        """
        Write parameter information to a CSV file.

        Args:
            filename (str): The name of the CSV file.
            directory (str): The directory where to save this CSV file.
            extension (str): The file extension, recommand using '.csv'.
        """
        # the parameters that will be a grid so it would not be in the main dat file
        # and they will be stored in stand alone files
        stand_alone_parameter_list = ['ratio_g2d', 'sigmad_ref']
        
        filename_full = os.path.join(directory, filename+extension)
        
        with open(filename_full, 'w') as file:
            file.write(f'# the model parameters of\n# %s \n'%(self.chemical_save_name))
            file.write(f'# name(str),value(str/float),units(str),valuetype(str),description(str)\n')
            for name, info in self.parameters.items():
                if name not in stand_alone_parameter_list: 
                    file.write(f"{name},{info['value']},{info['units']},{info['valuetype']},{info['description']}\n")
                elif name in stand_alone_parameter_list:
                    stand_alone_file_name = filename.replace('.', '_') + '_' + name + '.dat'
                    file.write(f"{name},{stand_alone_file_name},{info['units']},{info['valuetype']},{info['description']}\n")
                    np.savetxt(os.path.join(directory, stand_alone_file_name), info['value'])
                    
        return 1
                    

    @property
    def nphot(self):
        return int(self.get_parameter_value("nphot"))

    @nphot.setter
    def nphot(self, value):
        self.edit_parameter("nphot", new_value=value, new_units=None, new_description='Number of Photon Packages')

    @property
    def nthreads(self):
        return int(self.get_parameter_value("nthreads"))

    @nthreads.setter
    def nthreads(self, value):
        self.edit_parameter("nthreads", new_value=value, new_units=None, new_description='Number of Threads using in RADMC3D')

    @property
    def scattering_mode(self):
        return int(self.get_parameter_value("scattering_mode"))

    @scattering_mode.setter
    def scattering_mode(self, value):
        self.edit_parameter("scattering_mode", new_value=value, new_units=None, new_description='Scatter Mode in RADMC3D')

    @property
    def nr(self):
        return int(self.get_parameter_value("nr"))

    @nr.setter
    def nr(self, value):
        self.edit_parameter("nr", new_value=value, new_units=None, new_description='Number of Grid in r-direction')

    @property
    def ntheta(self):
        return int(self.get_parameter_value("ntheta"))

    @ntheta.setter
    def ntheta(self, value):
        self.edit_parameter("ntheta", new_value=value, new_units=None, new_description='Number of Grid in theta-direction')

    @property
    def nphi(self):
        return int(self.get_parameter_value("nphi"))

    @nphi.setter
    def nphi(self, value):
        self.edit_parameter("nphi", new_value=value, new_units=None, new_description='Number of Grid in phi-direction')

    @property
    def rin(self):
        return self.get_parameter_value("rin") * au

    @rin.setter
    def rin(self, value):
        self.edit_parameter("rin", new_value=value, new_units='au', new_description='Inner edge of the r grid -- needs to push to dust sublimation radius to solve VHSE correctly')

    @property
    def rout(self):
        return self.get_parameter_value("rout") * au

    @rout.setter
    def rout(self, value):
        self.edit_parameter("rout", new_value=value, new_units='au', new_description='Outer edge of the r grid -- needs to be large enough to include all the dust and CO gas')

    @property
    def thetaup(self):
        return self.get_parameter_value("thetaup")

    @thetaup.setter
    def thetaup(self, value):
        self.edit_parameter("thetaup", new_value=value, new_units=None, new_description='Upper edge of theta grid -- needs to be large enough to include all possible emission')

    @property
    def r_step_max(self):
        return self.get_parameter_value("r_step_max") * au

    @r_step_max.setter
    def r_step_max(self, value):
        self.edit_parameter("r_step_max", new_value=value, new_units='au', new_description='Maximum allowable step in the grid')
    
    @property
    def fmodel_filename(self):
        return self.get_parameter_value("fmodel_filename")

    @fmodel_filename.setter
    def fmodel_filename(self, value):
        self.edit_parameter("fmodel_filename", new_value=value, new_units=None, new_description='File name for the stellar spectrum in RADMC3D format', new_valuetype='str')

    @property
    def tstar(self):
        return self.get_parameter_value("tstar")

    @tstar.setter
    def tstar(self, value):
        self.edit_parameter("tstar", new_value=value, new_units='K', new_description='Stellar temperature -- not used now since the stellar spectrum is required')

    @property
    def mstar(self):
        return self.get_parameter_value("mstar") * ms

    @mstar.setter
    def mstar(self, value):
        self.edit_parameter("mstar", new_value=value, new_units='solar masses', new_description='Stellar mass')

    @property
    def rstar(self):
        return self.get_parameter_value("rstar") * 0.999999*rs

    @rstar.setter
    def rstar(self, value):
        self.edit_parameter("rstar", new_value=value, new_units='solar radius', new_description='Stellar radius')

    #  @property
    #  def pstar(self):
        #  return self.get_parameter_value("pstar")
#
    #  @pstar.setter
    #  def pstar(self, value):
        #  self.edit_parameter("pstar", new_value=value, new_units=None, new_description='Position of the star in the RADMC3D')

    @property
    def incl(self):
        return self.get_parameter_value("incl")

    @incl.setter
    def incl(self, value):
        self.edit_parameter("incl", new_value=value, new_units='degree', new_description='Inclination angle of the disk -- refer to RADMC3D defination')
    
    @property
    def phi(self):
        return self.get_parameter_value("phi")

    @phi.setter
    def phi(self, value):
        self.edit_parameter("phi", new_value=value, new_units=None, new_description='Observation position angle of the disk -- refer to RADMC3D defination')
   
    @property
    def mdotacc(self):
        return self.get_parameter_value("mdotacc") * ms/(31557600)

    @mdotacc.setter
    def mdotacc(self, value):
        self.edit_parameter("mdotacc", new_value=value, new_units='solar masses per year', new_description='Disk mass accretion rate')

    @property
    def mdiskd(self):
        return self.get_parameter_value("mdiskd") * ms

    @mdiskd.setter
    def mdiskd(self, value):
        self.edit_parameter("mdiskd", new_value=value, new_units='solar masses', new_description='Disk dust mass')
        
    @property
    def mdiskd_2(self):
        return self.get_parameter_value("mdiskd_2") * ms

    @mdiskd_2.setter
    def mdiskd_2(self, value):
        self.edit_parameter("mdiskd_2", new_value=value, new_units='solar masses', new_description='Secondary Disk dust mass')

    @property
    def ratio_g2d_global(self):
        return self.get_parameter_value("ratio_g2d_global")

    @ratio_g2d_global.setter
    def ratio_g2d_global(self, value):
        self.edit_parameter("ratio_g2d_global", new_value=value, new_units=None, new_description='Disk global gas-to-dust mass ratio -- this is the overall ratio represented as only one number')

    @property
    def ratio_g2d(self):
        return self.get_parameter_value("ratio_g2d")

    @ratio_g2d.setter
    def ratio_g2d(self, value):
        self.edit_parameter("ratio_g2d", new_value=value, new_units=None, new_description='Disk gas-to-dust mass ratio -- can be a number or a 2D or 3D array representing the spatial varying ratio')

    @property
    def sigmad_ref(self):
        return self.get_parameter_value("sigmad_ref")

    @sigmad_ref.setter
    def sigmad_ref(self, value):
        self.edit_parameter("sigmad_ref", new_value=value, new_units='g cm**(-2)', new_description='Disk reference sigmad -- an 2D array; radius_incm = sigmad_ref[:, 0] and sigmad_ref_ingcm-2 = sigmad_ref[:, 1]. If not set up; it uses pls and pltap to setup a disk surface density distribution with a tappering off radius')
    
    @property
    def pl_sufdens(self):
        return self.get_parameter_value("pl_sufdens")

    @pl_sufdens.setter
    def pl_sufdens(self, value):
        self.edit_parameter("pl_sufdens", new_value=value, new_units=None, new_description='Power law index for disk surface density distribution')
        
    @property
    def pl_sufdens_2(self):
        return self.get_parameter_value("pl_sufdens_2")

    @pl_sufdens_2.setter
    def pl_sufdens_2(self, value):
        self.edit_parameter("pl_sufdens_2", new_value=value, new_units=None, new_description='Power law index for Secondary disk surface density distribution')

    @property
    def pl_tapoff(self):
        return self.get_parameter_value("pl_tapoff")

    @pl_tapoff.setter
    def pl_tapoff(self, value):
        self.edit_parameter("pl_tapoff", new_value=value, new_units=None, new_description='Power law index for disk tappring off outer edge')

    @property
    def rdisk_in(self):
        return self.get_parameter_value("rdisk_in") * au

    @rdisk_in.setter
    def rdisk_in(self, value):
        self.edit_parameter("rdisk_in", new_value=value, new_units=None, new_description='the sharp-cut radius of the disk inner radius')

    @property
    def Rtap(self):
        return self.get_parameter_value("Rtap") * au

    @Rtap.setter
    def Rtap(self, value):
        self.edit_parameter("Rtap", new_value=value, new_units=None, new_description='Disk tappering off radius')
        
    @property
    def rdisk_out_2(self):
        return self.get_parameter_value("rdisk_out_2") * au

    @rdisk_out_2.setter
    def rdisk_out_2(self, value):
        self.edit_parameter("rdisk_out_2", new_value=value, new_units=None, new_description='Secondary Disk outer radius')

    @property
    def scaleheight_index(self):
        return self.get_parameter_value("scaleheight_index")

    @scaleheight_index.setter
    def scaleheight_index(self, value):
        self.edit_parameter("scaleheight_index", new_value=value, new_units=None, new_description='Different index tells the model to use different method to set up the scale height')
        
    @property
    def scaleheight_index_2(self):
        return self.get_parameter_value("scaleheight_index_2")

    @scaleheight_index_2.setter
    def scaleheight_index_2(self, value):
        self.edit_parameter("scaleheight_index_2", new_value=value, new_units=None, new_description='For Secondary disk: Different index tells the model to use different method to set up the scale height')

    @property
    def hp100(self):
        return self.get_parameter_value("hp100") * au

    @hp100.setter
    def hp100(self, value):
        self.edit_parameter("hp100", new_value=value, new_units='au', new_description='Pressure scale height at 100 au')
        
    @property
    def hp100_2(self):
        return self.get_parameter_value("hp100_2") * au

    @hp100_2.setter
    def hp100_2(self, value):
        self.edit_parameter("hp100_2", new_value=value, new_units='au', new_description='For Secondary disk: Pressure scale height at 100 au')

    @property
    def hprcr(self):
        return self.get_parameter_value("hprcr")

    @hprcr.setter
    def hprcr(self, value):
        self.edit_parameter("hprcr", new_value=value, new_units=None, new_description='Ratio of pressure scale height over tappering off radius at tappering off radius')

    @property
    def plh(self):
        return self.get_parameter_value("plh")

    @plh.setter
    def plh(self, value):
        self.edit_parameter("plh", new_value=value, new_units=None, new_description='Power law index of the pressure scale height radial dependence')
    
    @property
    def plh_2(self):
        return self.get_parameter_value("plh_2")

    @plh_2.setter
    def plh_2(self, value):
        self.edit_parameter("plh_2", new_value=value, new_units=None, new_description='For Secondary disk: Power law index of the pressure scale height radial dependence')

    @property
    def dustopacname_1(self):
        return str(self.get_parameter_value("dustopacname_1"))

    @dustopacname_1.setter
    def dustopacname_1(self, value):
        self.edit_parameter("dustopacname_1", new_value=value, new_units=None, new_description='Opacity file name for the dust opacity 1', new_valuetype='str')

    @property
    def dustopacname_2(self):
        return str(self.get_parameter_value("dustopacname_2"))

    @dustopacname_2.setter
    def dustopacname_2(self, value):
        self.edit_parameter("dustopacname_2", new_value=value, new_units=None, new_description='Opacity file name for the dust opacity 2', new_valuetype='str')

    @property
    def nr_dust_1(self):
        return int(self.get_parameter_value("nr_dust_1"))

    @nr_dust_1.setter
    def nr_dust_1(self, value):
        self.edit_parameter("nr_dust_1", new_value=value, new_units=None, new_description='Number of the dust species with the dust opacity 1')

    @property
    def nr_dust_2(self):
        return int(self.get_parameter_value("nr_dust_2"))

    @nr_dust_2.setter
    def nr_dust_2(self, value):
        self.edit_parameter("nr_dust_2", new_value=value, new_units=None, new_description='Number of the dust species with the dust opacity 2')

    @property
    def dust_spec_nr(self):
        return int(self.get_parameter_value("dust_spec_nr"))

    @dust_spec_nr.setter
    def dust_spec_nr(self, value):
        self.edit_parameter("dust_spec_nr", new_value=value, new_units=None, new_description='Total number of dust species')

    @property
    def amin_1(self):
        return self.get_parameter_value("amin_1")

    @amin_1.setter
    def amin_1(self, value):
        self.edit_parameter("amin_1", new_value=value, new_units='cm', new_description='Minimum grain size for dust with opacity 1')

    @property
    def amin_2(self):
        return self.get_parameter_value("amin_2")

    @amin_2.setter
    def amin_2(self, value):
        self.edit_parameter("amin_2", new_value=value, new_units='cm', new_description='Minimum grain size for dust with opacity 2')

    @property
    def amax_1(self):
        return self.get_parameter_value("amax_1")

    @amax_1.setter
    def amax_1(self, value):
        self.edit_parameter("amax_1", new_value=value, new_units='cm', new_description='Maximum grain size for dust with opacity 1')

    @property
    def amax_2(self):
        return self.get_parameter_value("amax_2")

    @amax_2.setter
    def amax_2(self, value):
        self.edit_parameter("amax_2", new_value=value, new_units='cm', new_description='Maximum grain size for dust with opacity 2')

    @property
    def amin_all(self):
        return self.get_parameter_value("amin_all")

    @amin_all.setter
    def amin_all(self, value):
        self.edit_parameter("amin_all", new_value=value, new_units='cm', new_description='Overall minimum grain size for dust')

    @property
    def amax_all(self):
        return self.get_parameter_value("amax_all")

    @amax_all.setter
    def amax_all(self, value):
        self.edit_parameter("amax_all", new_value=value, new_units='cm', new_description='Overall maximum grain size for dust')

    @property
    def pla_dustsize(self):
        return self.get_parameter_value("pla_dustsize")

    @pla_dustsize.setter
    def pla_dustsize(self, value):
        self.edit_parameter("pla_dustsize", new_value=value, new_units=None, new_description='Power law index for dust grain size distribution slope')

    @property
    def rhobulk(self):
        return self.get_parameter_value("rhobulk")

    @rhobulk.setter
    def rhobulk(self, value):
        self.edit_parameter("rhobulk", new_value=value, new_units='g cm^-3', new_description='Bulk density of dust grains')
        
    @property
    def visc_alpha(self):
        return self.get_parameter_value("visc_alpha")

    @visc_alpha.setter
    def visc_alpha(self, value):
        self.edit_parameter("visc_alpha", new_value=value, new_units=None, new_description='alpha parameter for dust settling (ref. Estrada+2016) and radial distribution')
        
    @property
    def vel_frag(self):
        return self.get_parameter_value("vel_frag")

    @vel_frag.setter
    def vel_frag(self, value):
        self.edit_parameter("vel_frag", new_value=value, new_units='cm s^-1', new_description='fragmentation threshold for dust grains')

    @property
    def radius_drift(self):
        return self.get_parameter_value("radius_drift") * au

    @radius_drift.setter
    def radius_drift(self, value):
        self.edit_parameter("radius_drift", new_value=value, new_units='au', new_description='The grains that needs to be drifted will be dirfted inside this radius')
    
    @property
    def a_drift(self):
        return self.get_parameter_value("a_drift")

    @a_drift.setter
    def a_drift(self, value):
        self.edit_parameter("a_drift", new_value=value, new_units='cm', new_description='The grains that are larger than this size will be drifted in')
        
    @property
    def r_wall_surf(self):
        return self.get_parameter_value("r_wall_surf") * au

    @r_wall_surf.setter
    def r_wall_surf(self, value):
        self.edit_parameter("r_wall_surf", new_value=value, new_units='au', new_description='the surface of the disk inner wall (rim) as well as the outer radii')
    
    @property
    def factor_wall_height(self):
        return self.get_parameter_value("factor_wall_height")

    @factor_wall_height.setter
    def factor_wall_height(self, value):
        self.edit_parameter("factor_wall_height", new_value=value, new_units='', new_description='the factor of how much the disk inner wall will be elevated compared to hp')
    
    #  @property
    #  def fracs(self):
        #  return list(self.get_parameter_value("fracs"))
#
    #  @fracs.setter
    #  def fracs(self, value):
        #  self.edit_parameter("fracs", new_value=value, new_units=None, new_description=None)
#
    #  @property
    #  def ndsd(self):
        #  return list(self.get_parameter_value("ndsd"))
#
    #  @ndsd.setter
    #  def ndsd(self, value):
        #  self.edit_parameter("ndsd", new_value=value, new_units=None, new_description=None)

    @property
    def chemical_save_dir(self):
        return str(self.get_parameter_value("chemical_save_dir"))

    @chemical_save_dir.setter
    def chemical_save_dir(self, value):
        self.edit_parameter("chemical_save_dir", new_value=value, new_units=None, new_description='Directory to save the chemical model', new_valuetype='str')

    @property
    def chemical_save_name(self):
        return str(self.get_parameter_value("chemical_save_name"))

    @chemical_save_name.setter
    def chemical_save_name(self, value):
        self.edit_parameter("chemical_save_name", new_value=value, new_units=None, new_description='Name of this model (and also the chemical part)', new_valuetype='str')

    @property
    def nr_cyl_LIME(self):
        return int(self.get_parameter_value("nr_cyl_LIME"))

    @nr_cyl_LIME.setter
    def nr_cyl_LIME(self, value):
        self.edit_parameter("nr_cyl_LIME", new_value=value, new_units=None, new_description='Number of r grid that will be used in LIME')

    @property
    def nz_cyl_LIME(self):
        return int(self.get_parameter_value("nz_cyl_LIME"))

    @nz_cyl_LIME.setter
    def nz_cyl_LIME(self, value):
        self.edit_parameter("nz_cyl_LIME", new_value=value, new_units=None, new_description='Number of z grid that will be used in LIME')

    """
    (devloping quantities)
    """
    @property
    def R_temp_trans(self):
        return self.get_parameter_value("R_temp_trans") * au

    @R_temp_trans.setter
    def R_temp_trans(self, value):
        self.edit_parameter("R_temp_trans", new_value=value, new_units='au', new_description='The temperature transition radius where the gas temperature deviates (get hotter) from the dust temperature within it')
        
    @property
    def fact_Tgas_2_Tdust(self):
        return self.get_parameter_value("fact_Tgas_2_Tdust")

    @fact_Tgas_2_Tdust.setter
    def fact_Tgas_2_Tdust(self, value):
        self.edit_parameter("fact_Tgas_2_Tdust", new_value=value, new_units=None, new_description='the factor to which the gas temperature is larger than the dust temperature within the temperature transition radius')
        
    @property
    def G0Hab_set(self):
        return self.get_parameter_value("G0Hab_set")

    @G0Hab_set.setter
    def G0Hab_set(self, value):
        self.edit_parameter("G0Hab_set", new_value=value, new_units=None, new_description='The FUV field G0 in Habbin unit at the stellar radius')
    
    #  @property
    #  def dummy(self):
        #  return self.get_parameter_value("dummy")

    #  @dummy.setter
    #  def dummy(self, value):
        #  self.edit_parameter("dummy", new_value=value, new_units=None, new_description=None)

    def set_default_parameters(self):
        """
        ##########################################
        # set up the default required parameters #
        ##########################################
        
        This function sets default values for a wide range of parameters used in the DiskMINT model. These parameters cover various aspects of the model, including DiskMINT model settings, RADMC3D setup, stellar parameters, disk parameters, dust kappa setup, and chemistry setup. The settings include numerical values, file paths, booleans, and array data. This function essentially acts as a comprehensive initializer for a DiskMINT model.
        
        The initial parameters we set up is for modeling RU Lup with a constant gas-to-dust mass ratio of 100 everywhere in the disk.
        We note that this is not the best-fit model on RU Lup but just set up here as an example.
        
        Remember that any parameter value set by set_default_parameters can be overridden by reading a parameter value from an external source, such as a CSV file, by the function of read_parameters_from_csv. Or it can be directly changed by changing the model class attributes before running the model.
        """
        
        """
        # 0. DiskMINT Model Settings
        # now all DiskMINT model setting is directly
        # set into the Mint class that for the model
        #
        # radmc3d running directory
        # currently can only work in this root direcotry
        #  self.file_dir = './'
        self.file_dir = os.getcwd()+'/'
        #
        # whether to start over from beginning
        # if we only want to change the gas-to-dust ratio
        # it will be fine to start from previous GS init setup to save time
        self.bool_startover  = True
        # if the startover set to True, all the following parameters would
        # be valid and need to be set up properly
        self.data_dir_readin = ''
        # the read_in_type can be 'gaussian' # or 'vhse'
        self.read_in_type    = ''
        self.prev_setup_name = ''
        #
        # whether created new Dust Kappa File
        # so that want to calculate dust parameters
        # such as fracs and aave instead of reading in
        self.bool_MakeDustKappa = False
        #
        # whether to compute the SED
        self.bool_SED           = True
        #
        # whether to use the VHSE distribution
        self.bool_VHSE          = True
        #
        # whether to run the chemical network
        self.bool_chemistry     = True
        #
        # whether to save the radmc3d model results in another folder
        self.bool_savemodel     = True
        """
        #################################################
        """
        1. RADMC3D Setup
        the basic setups for using RADMC3D, r.f., RADMC3D-2.0 manual
        including:
        radmc3d input
        radmc3d grid
        """
        #
        # RADMC3D input
        #
        # Number of Photon Packages
        self.nphot             = 10000000 # 10M photons
        #
        # Number of Threads using in radmc3d
        self.nthreads          = 8
        #
        # Scattering?
        self.scattering_mode   = 0 # as not including scattering effect
        # Grid parameters
        #
        # The Grid numbers at each axis (here polar axises in a sphere)
        self.nr                = 170 # prev 80; adopted 170 for normal models
        self.ntheta            = 200 # prev 150; adopted 200 for normal models
        self.nphi              = 1   # since now it is a 2d model, nphi = 1 fixed
        #
        # the edge of the grid
        #
        # 0.035 * au to push to dust sublimation radius for RU Lup
        self.rin               = 0.035 # * au
        #
        # rout needs to be large enough to include all the dust in the outer radius
        # 500*au to include as much disk as possible
        self.rout              = 500 # * au
        #
        # maximum allowable step in the grid
        # if 0, nr will be used and it will be an log space r grid
        # if used, once the diff_r is larger than this,
        # the r-grid changes to linear scale
        self.r_step_max        = 0. # * au
        #
        # thetaup needs to be large enough to include all the dust in the higher layers
        # especially for the very puffed up Gaussian disk
        # 1.2; 0.7 # adopted 0.7 for VHSE disks
        # if it is a more puffed up GS disk, it needs to be smaller to reach higher -- 0.2
        self.thetaup           = 0.7

        #################################################
        """
        2. Stellar Parameters
        """
        #
        # stellar spectrum
        # set fmodel_filename = '', if the stellar spectrum is not provided, and then tstar must be provoided
        # Caution: currently the tstar option is not working as chemistry need to read in UV spectrum.
        self.fmodel_filename   = 'RULupSpecModel_wuv.inp' # the stellar spectrum
        #
        # the Temperature of Star, not used if stellar spectrum is provided
        self.tstar             = 4060
        #
        # the Mass of Star 0.67 as from parameters Google Docs
        self.mstar             = 0.7 # * ms
        #
        # the Radius of Star
        # here 0.999999 is used so that the rin could be set as the rstar
        self.rstar             = 2.5 # * 0.999999 *rs
        #
        # the Position of Star
        #  self.pstar             = np.array([0.,0.,0.])
        #
        # disk accretion rate, only used when turn on the viscous heating
        # ms/year where 1 year = 31557600 seconds
        self.mdotacc           = 1e-7 # *ms/(31557600)

        #################################################
        """
        3. Disk Parameters
        """
        #
        # 3.1 Masses
        #
        # Disk dust mass
        self.mdiskd            = 4.0e-4 # * ms
        #
        # Disk gas mass
        # Gas Mass is set by the Gas-to-dust Mass Ratio
        # method 1: simply setup the ratio of g2d
        self.ratio_g2d_global  = 100.
        #
        # method 2: it could be determined from setting up Mgas
        # self.mdiskg            = 4.0e-2 * ms
        # self.ratio_g2d_global  = self.mdiskg/self.mdiskd
        #
        # method 3: set the ratio of g2d reference values as an array
        # with this method, the g2d can be set as a function of radius
        # **TO DO** set up the global ratio g2d accordingly to the radio_g2d
        # read_absolute_dir = '/home/dingshandeng/github/DiskModeling/6-model_IMLup/'
        # ratio_g2d         = np.loadtxt(read_absolute_dir+'ratio_g2d_reference.dat')
        # ratio_g2d_global  = 100. # the global value still needs to be identified manually
        #
        # 3.2 Powerlaw of the surface density
        #
        # set up the reference surface density distribution directly
        self.sigmad_ref        = []
        # surface density power index; usually as -1.0
        # (NOTE) in some codes it is called the column density,
        # surface density is used here to avoid potential confusion
        self.pl_sufdens        = -1.0
        #
        # tapering off power index:
        # the power index of the exponential decay after critical radius (rc = Rc = Rtap)
        # 1.0 as default value -- people usually used
        # from viscous evolving solution, it should be
        # 2+pl_sufdens
        self.pl_tapoff         = 2+self.pl_sufdens
        #
        # tappering Radius Rtap = rc
        self.Rtap              = 63.0 # * au
        #
        # disk inner edge, if not 0, disk has a sharp cut inner edge
        self.rdisk_in          = 0.0 # * au
        #
        # 3.3 Flaring of the disk
        # (Not Important if VHSE turned on)
        #
        # the index refering to how to setup the scaleheight
        # if scaleheight_index == 1, using hp100 as DIANA project
        # hp = hp100 * (r/rc)**plh
        # if scaleheight_index == 2, using hprcr = hp/Rc as DALI model
        # hp = hprcr * rc * (r/rc)**plh
        self.scaleheight_index = 1
        #
        # scale height at 100 au
        self.hp100             = 30 # * au
        #
        # or the scale height could be setup by the scale height at critical radius
        # hprcr: hp over rc ratio
        # hprcr = hp/Rc at Rc (Rtap)
        self.hprcr             = 0.1
        #
        # flaring index of the disk
        self.plh               = 1.10
        #
        # dust settling parameter (ref. Estrada+2016)
        # default is 1e-2 for minor settling
        # or 1e-3 and even 1e-4 for stronger settling
        self.visc_alpha        = 1.0e-2
        #
        # the fragmentation threshold for dust grains
        self.vel_frag          = 100
        #
        # the grain radial drifting size
        self.a_drift           = 5e-3
        # 
        # the grain radial drifting radius
        self.radius_drift      = 300 # au
        #
        # disk inner rim (inner disk wall)
        self.r_wall_surf        = 0.0
        self.factor_wall_height = 1.0
        #
        #
        # 3.4 geometry
        # 
        # disk inclination angle
        self.incl              = 0.0 # 18.8
        #
        # disk position angle
        self.pos               = 0.0 # -121.
        #
        # which angle to observe the disk
        # used in radmc3d
        self.phi               = 0.
        #
        # 2nd disk
        # [**Dev**]
        self.mdiskd_2          = 0.
        self.pl_sufdens_2      = self.pl_sufdens
        self.rdisk_out_2       = 0.
        self.hp100_2           = self.hp100
        self.plh_2             = self.plh

        #################################################
        """
        4. Dust Kappa Setup
        the dust opacity is now setup using dsharp_opac

        **TODO**
        a. put several pre-setup dust opacity files to read
        in and then to be used.
        b. add the feature to set up the opacity using optool
        c. maybe set up the opacity as a stand along object
        """

        #
        # DustKappa File Names and Number of Dust Species
        #
        # file name of dust species
        # dustopacity 1 is setup for small grains, 2 is for large as default
        # but these two can also be changed to other formats
        self.dustopacname_1    = ''
        self.dustopacname_2    = ''
        #
        # number of dust species
        # to setup our VHSE models
        # [NOTE] here we support two different dust species
        # which can be further divided into different subgroups
        # with different sizes.
        self.nr_dust_1         = 20
        self.nr_dust_2         = 0
        # automatically set the total nr of dust
        self.dust_spec_nr      = self.nr_dust_1 + self.nr_dust_2

        #
        # use following parameters to make dust kappa files
        #
        # here a is the radius of dust grain in the unit of cm.
        #
        # minimum size of the dust
        # use 100 AA as 1.0e-6 so that dust is not as small as molecule
        self.amin_1             = 1.0e-6
        self.amin_2             = 1.0e-6
        #
        # maximum size of the dust
        self.amax_1             = 0.3
        self.amax_2             = 0.3
        #
        # automatically set the overall smallest
        # and largest sizes among the two
        self.amin_all           = np.min([self.amin_1, self.amin_2])
        self.amax_all           = np.max([self.amax_1, self.amax_2])
        #
        # power law index of the dust size distribution
        self.pla_dustsize       = 3.5

        #
        # bulk density
        # the bulk density should be larger instead of this value
        # but it seems like it is going to be cancelled so double check the code
        # the rhobulk is used in chemical network
        # prev used 2.5 (too small)
        #
        self.rhobulk            = 3.224   # assume in the unit of [g cm^-3]
        
        """
        **TODO**
        to add the feature of set up fracs and ndsd and all other neccessary
        dust info from the parameter file. But for now, they are either calculated 
        from the function in dust_density.py or from input files .inp
        """
        # mass fraction of dust species
        # it must match with the number of dust species
        # and sum up as 1 to avoid any potential issue
        # set as empty list when it is set up automatically
        #  self.fracs              = []
        # [NOTE] if fracs is setup manually, then ndsd needs to be calculated manually as well

        """
        **TODO**
        add the function of computing ndsd accordingly with fracs
        
        **NOTE**
        the ndsd is now be setup at the dust setup section
        """
        # ndsd = number density * a**2
        #  self.ndsd               = []

        """
        **TODO**
        add the option to change the dust compositions (to be added)
        
        **NOTE**
        the dust composition is now changed in the main file with the utils
        to call the optool
        """
        #
        # Set additional Composition for Troilite
        #
        # fv_troilite = 0.0258 # actually this is currently muted in the wrapper

        #################################################
        """
        5. Chemistry Setup
        """
        #
        # The Name List
        # Tracking all the models have been done
        #

        self.chemical_save_dir  = ''
        self.chemical_save_name = ''

        #
        # The FUV field
        # if 0.0, it will be calculated from the UV spectra
        # carefull here because normally the UV spectra is tricky
        # and likely to be underestimated due to not-enough resolution
        # 
        
        self.G0Hab_set = 0.0 # None
        
        #
        # Grids in Cylindrical Coordinate and
        # For LIME Running and VHSE computation
        #
        """
        #
        # Now nr and nz are hard coded LIME wrapper, so don't change the part here
        #

        # the nr used here will be the same as used in RADMC3D
        # but nz will be slightly increased to match the
        # grid point in ntheta in RADMC3D
        """

        self.nr_cyl_LIME        = 100
        self.nz_cyl_LIME        = 200

        """
        NOTE: to change some parameters of the chemistry part,
        we need to change it directly inside the chemistry/.../*.py files
        For now
        """

        """
        **TO DO**
        have the function of readin and writeout parameters
        from and into the description dat file so that it would be
        easy to change parameters if people prefer to do so
        and the description will be easy as well be easy as well
        """
        
        """
        (developing parameters)
        """
        self.R_temp_trans = 0.
        self.fact_Tgas_2_Tdust = 1.
        
        return 1

"""
para = parameters()
para.read_parameters_from_csv(filename='../test_parameters.dat')
print(para.parameters)
print(para.mdiskd)
"""

class Mint:
    def __init__(self, parameters, file_dir=os.getcwd()):
        """
        # DiskMINT Model Settings
        """
        # radmc3d running directory
        # currently can only work in this root direcotry
        #  self.file_dir = './'
        self.file_dir = file_dir 
        #
        # whether to start over from beginning
        # if we only want to change the gas-to-dust ratio
        # it will be fine to start from previous GS init setup to save time
        self.bool_startover  = True
        # if the startover set to True, all the following parameters would
        # be valid and need to be set up properly
        self.data_dir_readin = ''
        # the read_in_type can be 'gaussian' # or 'vhse'
        self.read_in_type    = ''
        # the name (the name for chemical_save_name) of previous model
        self.prev_setup_name = ''
        # the previous setup gas-to-dust mass ratio
        # self.ratio_g2d_prev  = None
        #
        # whether created new Dust Kappa File
        # so that want to calculate dust parameters
        # such as fracs and aave instead of reading in
        self.bool_MakeDustKappa = False
        #
        # whether to compute the SED
        self.bool_SED           = True
        #
        # whether to use the VHSE distribution
        self.bool_VHSE          = True
        # the number of loops it tries to go through to solve the VHSE
        self.n_vhse_loop        = int(20) 
        #
        # whether to run the chemical network
        self.bool_chemistry     = True
        #
        # where is the executable code for chemical network
        """**TODO** put the chemical network to a bin file"""
        self.chem_code_dir      = '/home/dingshandeng/github/DiskModeling/0-wrapper/chemistry/with_CO2_ice/'
        #
        # whether to save the radmc3d model results in another folder
        self.bool_savemodel     = True
        
        """
        **NOTE**
        below are the bools that are still under development   
        """
        """START (DEV)"""
        # decouple the gas and dust in the inner disk
        self.bool_temp_decouple = False
        
        # 
        # whether to enalble dust settling
        # (Developing)
        self.bool_dust_settling = False
        
        #
        # whether to clean up the noise of dust temperature at midplane
        #
        self.bool_clean_dusttemp_midplane = True
        
        #
        # whether to relocate the dust (apply dust dynamics)
        # (Developing)
        # we have two options: fragmentation and drifting
        self.bool_dust_fragmentation = False
        self.bool_dust_radial_drifting = False
        #
        # whether to add a disk inner rim
        self.bool_dust_inner_rim = False
        #
        # whether to have the same rc grid for LIME and chemistry
        self.bool_same_rc_as_radmc3d = False
        
        """END (DEV)"""
        
        #
        # copy the parameters we setup to the model
        #
        para = self.parameters = parameters
        #  para = self.parameters

        # Set up distributions for density, temperature, and grid
        self.setup_radmc3d_grid(para)
        self.setup_dust_density_initial_gaussian(para, bool_initial_setup=True)
        
        # If a Secondary disk is added
        # if para.mdiskd_2 > 0:
        self.setup_secondary_dust_density_initial_gaussian(para, bool_initial_setup=True)
            
        self.setup_gas_density(para, bool_initial_setup=True)

        print('input primary mdiskd: %.3e [ms]'%(para.mdiskd/ms))
        print('setup primary mdiskd: %.3e [ms]'%(self.mdiskd_setup_1/ms))
        
        # if para.mdiskd_2 > 0:
        print('input secondary mdiskd_2: %.3e [ms]'%(para.mdiskd_2/ms))
        print('setup secondary mdiskd_2: %.3e [ms]'%(self.mdiskd_setup_2/ms))

        print('setup total mdiskd: %.3e [ms]'%(self.mdiskd_setup/ms))

        print('setup mdiskg: %.3e [ms]'%(self.mdiskg_setup/ms))
        print('setup gtd: %.3e'%(self.mdiskg_setup/self.mdiskd_setup))
        
        # set up wavelength
        self.setup_wavelengths(para)

        """
        #  the dust info can also be directly read in from inp files
        if self.bool_MakeDustKappa:
            self.setup_dust_info(para)
        #  elif len(para.ndsd)
        else:
            #  then the dust info needs to in read from .inp files
            self.readin_dust_info(file_dir=self.file_dir)
        """

    def setup_radmc3d_grid(self, para):
        # Use the parameters object to set up the radmc3d model grid

        ########################
        # Make the coordinates #
        ########################
        #
        # Griding can always be a issue
        #

        # what we really want is the grids based on collumn density
        self.ri, self.nr = dd.grid_radial_profile(para.rin, para.rout, para.nr, para.r_step_max)

        if hasattr(para, 'vertical_ncoarse'):
            #  print("Parameter 'vertical_ncoarse' exists in parameters.")
            pass
        else:
            #  print("Parameter 'vertical_ncoarse' does not exist in parameters.")
            para.vertical_ncoarse = int(0.4 * para.ntheta) # 20 # cells in the coarse regions, all the rest in fine region
        """
        **TO DO**
        add the vertical_ncoarse to default parameters
        also have the ability to set up the vertical emitting layer in the parameters
        """
        self.thetai  = dd.grid_vertical_layer(para.thetaup, para.ntheta, para.vertical_ncoarse)

        #
        # Phi (not important)
        #
        self.phii     = np.linspace(0.e0,np.pi*2.e0,para.nphi+1)

        #
        # Set the center values now
        #
        self.rc       = 0.5 * ( self.ri[1:] + self.ri[:-1] )
        self.thetac   = 0.5 * ( self.thetai[1:] + self.thetai[:-1] )
        self.phic     = 0.5 * ( self.phii[1:] + self.phii[:-1] )

        #
        # Renew the numbers of grids at each direction
        #""
        self.nr       = len(self.rc)
        self.ntheta   = len(self.thetac)
        self.nphi     = len(self.phic)

        #
        # Make the grid
        #
        self.qq       = np.meshgrid(self.rc,self.thetac,self.phic,indexing='ij')
        self.rr       = self.qq[0]
        self.tt       = self.qq[1] # theta -- vertical
        self.zr       = np.pi/2.e0 - self.qq[1]

        ####################################################################
        #### Setup the Coordinate Converting Info for VHSE Computation #####
        ####################################################################

        #  print('grid size: nr, ntheta = ', nr, ntheta)
        self.nr_cyl_insitu = self.nr * 3
        self.ntheta_cyl_insitu = self.ntheta * 3
        #  print('insitu cylindrical grid size: nr, ntheta = ', nr_cyl_insitu, ntheta_cyl_insitu)

        #############################################################
        #### Setup the Coordinate Converting Info for LIME grid #####
        #############################################################

        # the LIME coordinate info is set up in the
        # Parameter class, it can be modifed here, but be cautious
        # on doing so
        # [NOTE] this nr and nz are in cylindrical coordinate
        # and it needs to be matched with the LIME input file
        self.nr_cyl_LIME = para.nr_cyl_LIME
        self.nz_cyl_LIME = para.nz_cyl_LIME
        
        return 1
    

    def calculate_volume(self):
        """
        # Calculate the volume
        # calculate the volume for each cell in spherical coordinate
        """

        v_grid = self.v_grid = np.meshgrid(self.ri, self.thetai, self.phii[0], indexing = 'ij')
        vrr         = v_grid[0]
        vtt         = v_grid[1]
        dvrr        = (vrr[1:self.nr+1] - vrr[0:self.nr])[:, 0:self.ntheta]
        dvtt        = (vtt[:, 1:self.ntheta+1] - vtt[:, 0:self.ntheta])[0:self.nr, :]

        # rr is the mid value of the grid of vrr, which should be better than adopting vrr
        # the total vol should be what we are calculating now times 2
        # since only half of the disk (upper part) is calculated

        self.vol = dvrr * self.rr * np.sin(self.tt) * 2*np.pi * self.rr * dvtt
        
        return 1
    

    def setup_dust_density_initial_gaussian(self, para, bool_initial_setup=False):
        """
        set up the initial gaussian vertical distribution
        
        Args:
            bool_initial_setup (bool): if true, put the sigmad_setup to the class
        """
        
        #
        # Disk Parameters
        #
        # the sigmag0 is muted because we set up sigmad0 first,
        # then scale the sigmad0 to match the total dust mass
        # then scale the sigmad(r) with gtd(r) to have the setup sigmag(r)
        # the sigmad0 here is a dummy parameter which will be then scaled
        # to the set up disk dust mass (mdiskd)
        sigmad0  = 1e1 # sigmag0 / ratio_g2d  # Sigma dust at 1 AU
        plsig    = para.pl_sufdens # Powerlaw of the surface density
        plsig2   = para.pl_tapoff # tapering off power index

        scaleheight_index = para.scaleheight_index # which method to setup scale height
        hp100 = para.hp100
        hprcr = para.hprcr
        plh   = para.plh
        Rtap  = para.Rtap

        sigmad_reference = para.sigmad_ref
        #
        # Make Dust Density model according to Woitke2016
        #
        if len(sigmad_reference) == 0:
            print('Used the powerlow for surface density and tapering-off: %.3f; %.3f'%(plsig, plsig2))
            self.sigmad = sigmad0 * (self.rr*np.cos(self.zr)/au)**plsig * np.exp(-(self.rr*np.cos(self.zr)/Rtap)**(plsig2))
        else:
            print('Use the input reference surface density distribution to build the sigmad')
            self.sigmad = dd.place_surface_density(sigmad_reference, rc=self.rc, thetac=self.thetac, phic=self.phic, coord='sph')
        
        if para.rdisk_in > 0.0:
            self.sigmad[self.rr < para.rdisk_in] = 0.0
        
        if scaleheight_index == 1:
            print('Used the Method 1 to setup Hp, with Hp100 = %.3f'%(hp100/au))
            self.hp   = hp100 * (self.rr*np.cos(self.zr)/(100*au))**plh
        elif scaleheight_index == 2:
            print('Used the Method 2 to setup Hp, with hprcr = %.3f at Rc = %.3f'%(hprcr, Rtap/au))
            self.hp   = hprcr * Rtap * (self.rr*np.cos(self.zr)/(Rtap))**plh
        else:
            print('WRONG INPUT for the ScaleHeight Index, \n now switching to default value of using hp100 = 10 AU at 100AU')
            self.hp   = 10 * au * (self.rr*np.cos(self.zr)/(100*au))**plh

        if hasattr(self, 'vol'):
            pass
        else:
            self.calculate_volume()

        self.rhod = self.sigmad/(np.sqrt(2*np.pi) * self.hp) * np.exp(-((self.rr*np.sin(self.zr))**2 / ( 2 * self.hp**2)))
        
        mass = (self.rhod*self.vol).sum(0).sum(0).sum(0)
        self.sigmad_1 = self.sigmad = self.sigmad * (para.mdiskd/2/mass)
        self.rhod_1   = self.rhod   = self.rhod * (para.mdiskd/2/mass)
        
        if bool_initial_setup:
            self.sigmad_setup = self.sigmad_1_setup = self.sigmad
            self.mdiskd_setup_half_1 = self.mdiskd_setup_half = (self.rhod*self.vol).sum(0).sum(0).sum(0)
            self.mdiskd_setup_1      = self.mdiskd_setup      = 2 * self.mdiskd_setup_half   
        
        return 1
        
        
    def setup_secondary_dust_density_initial_gaussian(self, para, bool_initial_setup=False):
        """
        **Devloping**
        set up the initial gaussian vertical distribution
        for Secondary dust disk Zone
        this is ideal for modeling transition disks and build Two-Zone model
        
        the secondary dust disk is normally the disk with smaller disk mass, thus the inner disk of a transit disk.
        
        Args:
            bool_initial_setup (bool): if true, put the sigmad_setup to the class
        """
        
        #
        # Disk Parameters
        #
        # the sigmag0 is muted because we set up sigmad0 first,
        # then scale the sigmad0 to match the total dust mass
        # then scale the sigmad(r) with gtd(r) to have the setup sigmag(r)
        # the sigmad0 here is a dummy parameter which will be then scaled
        # to the set up disk dust mass (mdiskd)
        
        #
        # set up the secondary secondary dust disk when it is not 0
        # 
        if para.mdiskd_2 > 0.0:
            sigmad0  = 1e1 # sigmag0 / ratio_g2d  # Sigma dust at 1 AU
            plsig    = para.pl_sufdens_2 # Powerlaw of the surface density
            # plsig2   = para.pl_tapoff_2 # tapering off power index

            scaleheight_index = para.scaleheight_index_2 # which method to setup scale height
            hp100 = para.hp100_2
            plh   = para.plh_2
            Rout  = para.rdisk_out_2

            #
            # Make Dust Density model according to Woitke2016
            #
            print('Used the powerlow for surface density: %.3f;'%(plsig))
            self.sigmad_2 = sigmad0 * (self.rr*np.cos(self.zr)/au)**plsig
            self.sigmad_2[self.rr > Rout] = 0.0 # sharp-cut the outer rim of the secondary dust disk.
            
            if scaleheight_index == 1:
                print('Used the Method 1 to setup Hp, with Hp100 = %.3f'%(hp100/au))
                self.hp   = hp100 * (self.rr*np.cos(self.zr)/(100*au))**plh
            else:
                print('WRONG INPUT for the ScaleHeight Index, \n now switching to default value of using hp100 = 10 AU at 100AU')
                self.hp   = 10 * au * (self.rr*np.cos(self.zr)/(100*au))**plh

            if hasattr(self, 'vol'):
                pass
            else:
                self.calculate_volume()

            self.rhod_2 = self.sigmad_2/(np.sqrt(2*np.pi) * self.hp) * np.exp(-((self.rr*np.sin(self.zr))**2 / ( 2 * self.hp**2)))
            mass = (self.rhod_2*self.vol).sum(0).sum(0).sum(0)
            self.sigmad_2 = self.sigmad_2 * (para.mdiskd_2/2/mass)
            #  sigmad_f = sigmad.copy()
            self.rhod_2 = self.rhod_2 * (para.mdiskd_2/2/mass)

            if bool_initial_setup:
                self.mdiskd_setup_half_2 = (self.rhod_2*self.vol).sum(0).sum(0).sum(0)
                self.mdiskd_setup_2      = 2 * self.mdiskd_setup_half_2
        
        else:
            self.rhod_2   = self.rhod.copy() * 0.0
            self.sigmad_2 = self.sigmad.copy() * 0.0
            if bool_initial_setup:
                self.mdiskd_setup_half_2 = 0.0
                self.mdiskd_setup_2      = 0.0
        
        # save a copy of the primary disk
        self.rhod_1 = self.rhod.copy()
        self.sigmad_1 = self.sigmad.copy()
        
        if bool_initial_setup:
            self.mdiskd_setup_1 = self.mdiskd_setup.copy()
            self.mdiskd_setup_half_1 = self.mdiskd_setup_half.copy()
        
        # combine the primary with the secondary disk
        self.rhod = self.rhod + self.rhod_2 
        self.sigmad = self.sigmad + self.sigmad_2
        if bool_initial_setup:
            self.sigmad_2_setup = self.sigmad_2
            self.sigmad_setup = self.sigmad
            self.mdiskd_setup = self.mdiskd_setup + self.mdiskd_setup_2
            self.mdiskd_setup_half = self.mdiskd_setup_half + self.mdiskd_setup_half_2
        
        return 1
    

    def setup_g2d_grid(self, para, bool_initial_setup=False):
        """
        # set up a grid of g2d in the shape of the model grid
        
        Args:
            bool_initial_setup (bool): if true, put the ratio_g2d_setup to the class
        """
        
        if hasattr(para, 'ratio_g2d') and para.ratio_g2d is not None:
            #  print('case 1')
            ratio_g2d_set = para.ratio_g2d
        else:
            #  print('case 2')
            ratio_g2d_set = para.ratio_g2d_global
        
        """
        dev
        put up a grid of g2d in the DiskMINT model
        """
        # ratio_g2d_set = copy.deepcopy(ratio_g2d)
        #  print(ratio_g2d_set)
        # ratio_g2d     = dd.place_ratio_g2d(ratio_g2d_set, self.sigmad.shape, self.rc)
        ratio_g2d_grid = dd.place_ratio_g2d(ratio_g2d_set, rc=self.rc, thetac=self.thetac, phic=self.phic, coord='sph', ratio_g2d_outgrid=para.ratio_g2d_global, rc_cyl_set=self.rc)
        
        self.ratio_g2d = ratio_g2d_grid
        
        if bool_initial_setup:
            self.ratio_g2d_setup = self.ratio_g2d
            self.sigmagas_setup = self.sigmad_setup * self.ratio_g2d_setup
        
        return 1
    

    def setup_gas_density(self, para, bool_initial_setup=False):
        """
        # set up the gas density accordingly to the dust disk
        
        Args:
            bool_initial_setup (bool): if true, put the rhogas_setup (not used now) and the mdiskg_setup to the code 

        Returns:
            add rhogas
            and mdiskg_setup to the mint class
            
        Description:
            TBA
        """
        
        if hasattr(self, 'ratio_g2d'):
            pass
        else:
            self.setup_g2d_grid(para, bool_initial_setup=bool_initial_setup)
        
        self.rhogas = self.rhod * self.ratio_g2d

        if hasattr(self, 'vol'):
            pass
        else:
            self.calculate_volume()
        
        if bool_initial_setup:
            # # record the initial rhogas setup
            # self.rhogas_setup = self.rhogas
            # calculate the total disk gas mass in the volume
            self.mdiskg_setup = (self.rhogas*self.vol).sum(0).sum(0).sum(0) * 2
        
        return 1
    

    def setup_wavelengths(self, para):
        """
        set up the wavelengths range for computation in RADMC3D SED
        """
        
        ######################################
        #### Wavelengths for Calculation #####
        ######################################
        
        if para.fmodel_filename != '':
            fmodel = np.loadtxt(os.path.join(self.file_dir,para.fmodel_filename))
            self.lam = fmodel[:,0]
            self.Fnu = fmodel[:,1]
            self.nlam = len(fmodel[:,0])
        else:
            #
            # Write the wavelength_micron.inp file
            # Using the default wavelengths setup number
            #
            # There are three ranges of lambda and they are connected together
            lam1     = 0.1e0
            # lam2     = 7.0e0
            # lam3     = 25.e0
            lam4     = 1.0e4
            # n12      = 20
            # n23      = 100
            # n34      = 30
            n14      = 100
            # lam12    = np.logspace(np.log10(lam1),np.log10(lam2),n12,endpoint=False)
            # lam23    = np.logspace(np.log10(lam2),np.log10(lam3),n23,endpoint=False)
            # lam34    = np.logspace(np.log10(lam3),np.log10(lam4),n34,endpoint=True)
            lam14    = np.logspace(np.log10(lam1),np.log10(lam4),n14,endpoint=True)
            # lam      = np.concatenate([lam12,lam23,lam34])
            self.lam      = lam14
            self.nlam     = self.lam.size
            
        return 1
    

    def setup_dust_info(self, para, bool_savefile=False):
        """
        # calculate the neccessary dust information here
        # it can also be setup by directly reading dust parameters file
        # the info that is needed including
        # fracs: mass fractions of dust species
        # a_ave: average dust radius of dust dpecies
        # fracs_nb: the number fractions of dust species
        # ndsd: nd(a) * pi * a**2
        # all these parameters should be setup when generating the dust
        # opacity files to prevent any mistakes, which we recommand to do
        # but we also keep the option of calculating them here
        """

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

        pla     = para.pla_dustsize
        # q_h     = para.pla_dustsize
        rhobulk = para.rhobulk

        #
        # Save fracs and a_ave Part
        #
        a_ave    = []
        fracs_nb = []
        """**TODO** [DEV] enalbling the feature of set the fracs in the parameters file"""
        fracs    = [] # para.fracs # fracs can be setup in the parameters file
        ndsd     = [] # para.fracs # ndsd can be setup in the parameters file;
        # if the fracs is provided, then ndsd needs to be provided manually

        # ai_sm = np.logspace(np.log10(para.amin_1), np.log10(para.amax_1), para.nr_dust_1 + 1)
        # ai_bg = np.logspace(np.log10(para.amin_2), np.log10(para.amax_2), para.nr_dust_2 + 1)
        # ai    = np.hstack([ai_sm[:-1], ai_bg])
        ai_1 = np.logspace(np.log10(para.amin_1), np.log10(para.amax_1), para.nr_dust_1 + 1)
        ai_2 = np.logspace(np.log10(para.amin_2), np.log10(para.amax_2), para.nr_dust_2 + 1)
        """
        **TODO**
        in current version, only one uniform dust composition
        is supported, will add more versatile treatment in the future
        """
        # ai = ai_sm.copy()

        #  print('input paras:', amin_1, amax_1, amin_2, amax_2, nr_dust_1, nr_dust_2, pla)

        
        if len(fracs) == 0:
            print('Computing Frations of Dust Species Assuming in the same Slope')
            
            for ai in [ai_1, ai_2]:
                mass_dustfrac_all = massdust_i(np.max(ai), np.min(ai), rhobulk = rhobulk, pla=pla)
                numb_dustfrac_all = numbdust_i(np.max(ai), np.min(ai), rhobulk = rhobulk, pla=pla)

                # set up the fractions for each dust sub-species
                for i in range(len(ai) - 1):
                    ai_max = ai[i+1]
                    ai_min = ai[i]

                    a_ave.append(a_average(ai_max, ai_min, pla=pla))
                    fracs.append(massdust_i(ai_max, ai_min, rhobulk = rhobulk, pla=pla)/mass_dustfrac_all)
                    fracs_nb.append(numbdust_i(ai_max, ai_min, rhobulk = rhobulk, pla=pla)/numb_dustfrac_all)
                    ndsd.append(ndsd_i(ai_max, ai_min, rhobulk = rhobulk, pla=pla))

        else:
            """
            **TODO** the part below might not work as designed so far.
            """
            print('Using User Provided Dust Mass Fractions as:', fracs)
            
            for ai in [ai_1, ai_2]:
                # set up the fractions for each dust sub-species
                for i in range(len(ai) - 1):
                    ai_max = ai[i+1]
                    ai_min = ai[i]
                    a_ave_t = a_average(ai_max, ai_min, pla=pla)
                    a_ave.append(a_ave_t)
                    fracs_nb.append(fracs[i]/(4/3*np.pi*rhobulk*a_ave_t**3))

            fracs_nb = np.array(fracs_nb)/np.sum(np.array(fracs_nb))

        self.a_ave    = a_ave.copy()
        self.fracs    = fracs.copy()
        self.fracs_nb = fracs_nb.copy()
        self.ndsd     = ndsd.copy()
        
        self.dust_spec_nr = len(fracs)

        if bool_savefile:
            np.savetxt('aave.inp', a_ave)
            np.savetxt('fracs.inp', fracs)
            np.savetxt('fracs_numb.inp', fracs_nb)
            np.savetxt('ndsd.inp', ndsd)
            
        return 1
    

    def readin_dust_info(self, file_dir, \
                         file_aave='aave.inp', file_fracs='fracs.inp',\
                         file_fracs_numb='fracs_numb.inp', file_ndsd='ndsd.inp'):
        """
        # readin the neccessary dust information here
        # it can also be setup by directly reading dust parameters file
        # the info that is needed including
        # fracs: mass fractions of dust species
        # a_ave: average dust radius of dust dpecies
        # fracs_nb: the number fractions of dust species
        # ndsd: nd(a) * pi * a**2
        # all these parameters should be setup when generating the dust
        # opacity files to prevent any mistakes, which we recommand to do
        """
        
        #################################
        #### Read Dust Species Info #####
        #################################

        a_ave    = np.loadtxt(os.path.join(file_dir,file_aave)) # mass fracs of dust species
        fracs    = np.loadtxt(os.path.join(file_dir,file_fracs)) # average dust radius of dust species
        fracs_nb = np.loadtxt(os.path.join(file_dir,file_fracs_numb)) # number fracs of dust species
        ndsd     = np.loadtxt(os.path.join(file_dir,file_ndsd)) # n_d(a) pi a^2

        if len(fracs.shape) == 0:
            fracs = np.array([fracs])
            
        if len(a_ave.shape) == 0:
            a_ave = np.array([a_ave])
            
        if len(fracs_nb.shape) == 0:
            fracs_nb = np.array([fracs_nb])
            
        if len(ndsd.shape) == 0:
            ndsd = np.array([ndsd])

        self.a_ave    = a_ave.copy()
        self.fracs    = fracs.copy()
        self.fracs_nb = fracs_nb.copy()
        self.ndsd     = ndsd.copy()

        self.dust_spec_nr = len(fracs)

        return 1
