#!/usr/bin/env python
"""
Generate parameter grid CSV files for DiskMINT model runs
"""

import pandas as pd
import numpy as np
import os
import argparse

def generate_grid(output_file, grid_type='default', mstar="unknown", model_dir="models_example"):
    """
    Generate parameter grid CSV file
    
    grid_type options:
    - 'default': Original 12-model grid (3 masses × 4 g2d)
    - 'extended': Larger grid with more combinations
    - 'custom': Define your own ranges
    """
    
    grid_setup_dict = None
    ratio_mdisk_values = None
    ratio_g2d_values = None
    rtap_values = None
    
    if grid_type == 'default':
        # Your original grid
        ratio_mdisk_values = [1e-4, 1e-3]
        ratio_g2d_values = [10., 30., 100., 300.]
        rtap_values = [100.]
        
    elif grid_type == 'extended':
        # Extended grid
        ratio_mdisk_values = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
        ratio_g2d_values = [10., 30., 100., 300.]
        rtap_values = [10., 30., 100., 300.]
        
    elif grid_type == 'test':
        # Quick test grid
        ratio_mdisk_values = [1e-4, 5e-3]
        ratio_g2d_values = [100.]
        rtap_values = [100.]
        
    elif grid_type == 'supplementary':
        # # Quick test grid
        # ratio_mdisk_values = [1e-4, 1e-2]
        # ratio_g2d_values = [100.]
        # rtap_values = [100.]
        
        print('Building up the grid using the supplementary grid')
        
        # # for 0.3 ms on 02/20/2026
        # grid_setup_dict = {
        #     "ratio_mdisk_values": [1e-4, 1e-2],
        #     "ratio_g2d_values": [300, 100],
        #     "rtap_values": [10, 10]
        # }
        
        # for 0.7 ms on 02/22/2026
        # mdisk/mstar: 1.00e-03, gtd: 10, rc: 300 AU
        grid_setup_dict = {
            "ratio_mdisk_values": [1e-3],
            "ratio_g2d_values": [10],
            "rtap_values": [300]
        }
    
    def _get_model_mark_short(mstar, ratio_mdisk, ratio_g2d, rtap):
        return f"diskmint_mstar{mstar}ms_mdisktmstar{ratio_mdisk:.2e}_gtd{int(ratio_g2d)}_rtap{int(rtap)}au".replace('.', 'p')
    
    grid_data = []
    if grid_setup_dict is not None:
        for i in range(len(grid_setup_dict["ratio_mdisk_values"])):
            
            ratio_mdisk = grid_setup_dict["ratio_mdisk_values"][i]
            ratio_g2d = grid_setup_dict["ratio_g2d_values"][i]
            rtap = grid_setup_dict["rtap_values"][i]
            
            model_mark_short_t = _get_model_mark_short(mstar, ratio_mdisk, ratio_g2d, rtap)
            model_directory_t = os.path.join(model_dir, model_mark_short_t)
            
            grid_data.append({
                'ratio_mdisk': ratio_mdisk,
                'ratio_g2d': ratio_g2d,
                'rtap': rtap,
                'model_mark_short': model_mark_short_t,
                'model_directory': model_directory_t
            })
    
    else: 
        # Generate all combinations
        for ratio_mdisk in ratio_mdisk_values:
            for ratio_g2d in ratio_g2d_values:
                for rtap in rtap_values:
                    
                    model_mark_short_t = _get_model_mark_short(mstar, ratio_mdisk, ratio_g2d, rtap) 
                    model_directory_t = os.path.join(model_dir, model_mark_short_t)
                    
                    grid_data.append({
                        'ratio_mdisk': ratio_mdisk,
                        'ratio_g2d': ratio_g2d,
                        'rtap': rtap,
                        'model_mark_short': model_mark_short_t,
                        'model_directory': model_directory_t
                    })
    
    df = pd.DataFrame(grid_data)
    df.to_csv(output_file, index=False)
    
    print(f"Generated grid with {len(df)} models")
    print(f"Saved to: {output_file}")
    print("\nGrid summary:")
    
    if grid_setup_dict is not None:
        print('build with a dict grid')
        print(f"  Mdisk/Mstar: {grid_setup_dict['ratio_mdisk_values']}")
        print(f"  g2d ratios: {grid_setup_dict['ratio_g2d_values']}")
        print(f"  Rtap values: {grid_setup_dict['rtap_values']}")
    
    if ratio_mdisk_values is not None:
        print('build with the mesh grid')
        print(f"  Mdisk/Mstar: {ratio_mdisk_values}")
        print(f"  g2d ratios: {ratio_g2d_values}")
        print(f"  Rtap values: {rtap_values}")
    
    return len(df)

def generate_file_trees(model_dir, grid_csv):
    """
    Generate file trees for each model based on the grid CSV file
    
    For each model defined in the grid CSV, create a directory structure like:
    model_dir/
        diskmint_mstar{mstar}ms_mdisk{ratio_mdisk:.2e}ms_gtd{int(ratio_g2d)}_rtap{int(rtap)}/
            data/
            output/
        diskmint_xxx/
            data/
            output/
        ...
    
    This function reads the grid CSV to determine how many models there are and creates the corresponding directories.
    """
    # Read the grid CSV to determine the number of models
    df = pd.read_csv(grid_csv)
    num_models = len(df)
    df.index = np.arange(num_models)  # Ensure index starts from 0 to num_models-1
    
    for i in df.index:
        if 'model_directory' in df.columns:
            model_path = df.loc[i, 'model_directory']  # Use the model_directory if it exists
        else:
            model_t = df.loc[i, 'model_mark_short']  # Use the model_mark_short as the unique identifier
            model_path = os.path.join(model_dir, model_t)
        
        # Create data and output directories for each model
        data_path = os.path.join(model_path, 'data')
        output_path = os.path.join(model_path, 'output')
        
        os.makedirs(data_path, exist_ok=True)
        os.makedirs(output_path, exist_ok=True)
    
    print(f"Generated file trees for {num_models} models in {model_dir}") 

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default=None, help='Output CSV file')
    parser.add_argument('--grid', default='test', choices=['default', 'extended', 'test', 'supplementary'],
                       help='Grid type to generate')
    parser.add_argument('--mstar', default='unknown', help='Stellar mass for model marking')
    parser.add_argument('--model_dir', default=None, help='Base model directory')
    parser.add_argument('--generate_file_trees', action='store_true', default=False, help='Generate file trees for each model')
    args = parser.parse_args()
    
    if args.output is None:
        mstar_str = "%.1f"%(float(args.mstar)) if args.mstar != 'unknown' else 'unknown'
        output_filename = f"grids/grid_{mstar_str}ms_{args.grid}".replace('.', 'p') + ".csv"
        print(f'Specify the default output file as {output_filename}')
        args.output = output_filename
    
    generate_grid(args.output, args.grid, mstar=args.mstar, model_dir=args.model_dir)
    
    if args.model_dir is not None and args.generate_file_trees:
        generate_file_trees(args.model_dir, args.output)
    