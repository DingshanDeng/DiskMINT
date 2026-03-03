#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author: Dingshan Deng @ University of Arizona
# contact: dingshandeng@arizona.edu

"""
Refactored DiskMINT model runner for HPC array jobs
"""

import os, sys, copy, shutil, glob
import numpy as np
import platform
import time
import argparse

def main():
    parser = argparse.ArgumentParser(description='Run DiskMINT model with specified parameters')
    parser.add_argument('--mstar', type=float, required=True, help='Stellar mass in solar masses')
    parser.add_argument('--ratio_mdisk', type=float, required=True, help='Disk-to-star mass ratio')
    parser.add_argument('--ratio_g2d', type=float, required=True, help='Gas-to-dust ratio')
    parser.add_argument('--rtap', type=float, required=True, help='Tapering radius in AU')
    parser.add_argument('--model_date', type=str, required=True, help='Model date identifier')
    parser.add_argument('--model_mark', type=str, required=True, help='Model mark/label')
    parser.add_argument('--model_id', type=int, required=True, help='Model ID number')
    parser.add_argument('--param_file', type=str, default=None, help='Base parameter CSV file')
    parser.add_argument('--model_dir', type=str, default=None, help='Model output directory where data/ and output/ will be created and the model will be run. If not specified, it will be created in the current working directory as an example.')
    
    args = parser.parse_args()
    
    time0 = time.time()
    
    # Package setup
    if platform.node() == 'garnet':
        package_position = "/home/dingshandeng/github/DiskMINT"
    else: # this is on hpc
        package_position = "/home/u24/dingshandeng/github/DiskMINT"
        
    sys.path.append(os.path.join(package_position, 'src'))
    import diskmint.constants as const
    import diskmint.model as model
    import diskmint.disk_density as dd
    import diskmint.execute as exe
    
    # Directory setup
    # workding_dir is where we stored the parameter_grid.csv, the input/ and the parameter.csv
    # it is also where the job is submitted
    working_dir = os.getcwd()
    read_dir = os.path.join(working_dir, 'input')
    # model_dir is where the current model will be running
    model_dir = args.model_dir if args.model_dir is not None else os.path.join(working_dir, 'diskmint_radmc3d_model_example')
    data_dir = os.path.join(model_dir, 'data')
    save_dir = os.path.join(model_dir, 'output')
    
    # Create necessary directories
    for dir_path in [model_dir, data_dir, save_dir]:
        if not os.path.exists(dir_path):
            print(f"Directory {dir_path} does not exist. Creating it.")
            os.makedirs(dir_path, exist_ok=True)
            print(f"Directory created: {dir_path}")
        else:
            print(f"Directory already exists: {dir_path}")
            print(f"Caution: this run is overwriting existing data inside /data.")
            print(f"Caution: make sure this is intended to avoid data loss.")
    
    # Model options
    bool_MakeDustKappa = True
    bool_SED = False
    bool_VHSE = True
    n_vhse_loop = 20
    bool_dust_settling = True
    bool_clean_dusttemp_midplane = True
    bool_chemistry = False
    bool_savemodel = True
    
    # Load base parameters
    para = model.Parameters()
    # the paramter csvs are saved in the model dir
    # the format is 
    # model_dir/
    # model_parameters_template.csv
    # -- data/ # where the current model is running
    # -- output/ # save the finished models
    # -- -- model_name/
    # -- -- -- all the files for this model
    if args.param_file is not None:
        para.read_parameters_from_csv(filename=args.param_file, directory=working_dir, extension='')
    else:
        para.read_parameters_from_csv(filename=("model_parameters_%.1fMsolar.csv" % args.mstar).replace('.', 'p'), directory=working_dir, extension='')
    
    # Update parameters
    para.edit_parameter("chemical_save_dir", new_value=save_dir)
    para.edit_parameter("nphot", new_value=1e7)
    para.edit_parameter("nthreads", new_value=94)
    
    # Calculate and set mass parameters
    mdisk_gas = args.ratio_mdisk * args.mstar
    mdisk_dust = mdisk_gas / args.ratio_g2d
    para.edit_parameter("mdiskd", new_value=mdisk_dust)
    para.edit_parameter("ratio_g2d_global", new_value=args.ratio_g2d)
    para.edit_parameter("Rtap", new_value=args.rtap)
    
    # Set G0Hab based on stellar mass (you can extend this)
    g0hab_values = {
        0.1: 2.136e+09,
        0.3: 7.395e+09,
        0.5: 1.822e+10,
        0.7: 3.480e+10,
        1.0: 6.343e+10,
        2.0: 2.081e+11
    }
    para.edit_parameter("G0Hab_set", new_value=g0hab_values.get(args.mstar, 1.0e+10))
    
    # Construct model name
    model_mark = args.model_mark
    if 'DSHARP' in para.dustopacname_1:
        model_mark += '_DSHARPdust'
    if bool_dust_settling:
        model_mark += '_dustset'
    
    name_of_this_model = (
        f"diskmintmodel{args.model_date}_t{args.model_id}_{model_mark}_"
        f"mgas{mdisk_gas:.2e}ms_mdust{mdisk_dust:.2e}ms_"
        f"gtd{args.ratio_g2d:.0f}_rc{args.rtap:.1f}_"
        f"alphav{para.get_parameter_value('visc_alpha'):.1e}"
    ).replace('.', 'p')
    para.edit_parameter("chemical_save_name", new_value=name_of_this_model)
    
    print('=' * 30)
    print(f'Model #{args.model_id}: {name_of_this_model}')
    print(f'Mstar={args.mstar} Msun, Mdisk/Mstar={args.ratio_mdisk}, g2d={args.ratio_g2d}')
    print('=' * 30)
    
    # Change to data directory
    # this is where the model will be run and all the input files should be copied to this directory
    os.chdir(data_dir)
    
    # Copy input files
    chem_code_dir = os.path.join(package_position, 'chemistry', 'reducedRGH22')
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
    
    # Set up and run model
    mint = model.Mint(para, file_dir=data_dir)
    mint.bool_MakeDustKappa = bool_MakeDustKappa
    mint.bool_SED = bool_SED
    mint.bool_VHSE = bool_VHSE
    mint.n_vhse_loop = n_vhse_loop
    mint.bool_dust_settling = bool_dust_settling
    mint.bool_clean_dusttemp_midplane = bool_clean_dusttemp_midplane
    mint.bool_chemistry = bool_chemistry
    mint.bool_savemodel = bool_savemodel
    mint.chem_code_dir = chem_code_dir
    
    mint.bool_temp_decouple = False
    mint.bool_dust_fragmentation = False
    mint.bool_dust_radial_drifting = False
    mint.bool_dust_inner_rim = False
    mint.bool_same_rc_as_radmc3d = True
    
    # Run the model
    exe.runmodel(mint, test_alliteration=False)
    
    print('Model run completed. Saving results...')
    print('Results are saving to the output/ directory inside the model directory.')
    
    # Save CO data files
    data_save_dir = os.path.join(save_dir, para.chemical_save_name)
    os.makedirs(data_save_dir, exist_ok=True)
    
    for file_pattern in ['COinitgrid-*.dat', 'COendgrid-*.chem']:
        for f in glob.glob(os.path.join(data_dir, file_pattern)):
            shutil.move(f, data_save_dir)
    
       
    elapsed_time = time.time() - time0
    print('=' * 30)
    print(f'Model #{args.model_id} completed')
    print(f'Wall time: {elapsed_time/3600:.2f} hours')
    print(f'CPU hours: {elapsed_time * para.nthreads / 3600:.2f}')
    print('=' * 30)
    
    # clean up the data directory to save space
    def empty_dir(path: str) -> None:
        path = os.path.abspath(path)
        if not os.path.isdir(path):
            raise NotADirectoryError(path)

        # Delete everything inside, but keep the folder itself
        for name in os.listdir(path):
            p = os.path.join(path, name)
            try:
                if os.path.isdir(p) and not os.path.islink(p):
                    shutil.rmtree(p)
                else:
                    os.remove(p)  # files + symlinks
            except FileNotFoundError:
                # In case something else deletes it concurrently
                pass
    
    os.chdir(model_dir)
    empty_dir(data_dir)
    print(f"Emptied contents of: {data_dir}")

if __name__ == '__main__':
    main()
