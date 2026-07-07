#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author: Dingshan Deng @ University of Arizona
# contact: dingshandeng@arizona.edu
# created: 07/07/2026

"""
Example: fit continuum and C18O radial profiles by changing surface density
and gas-to-dust ratio.

Run from this directory:
python -u 1-model_radial_profile_fit_advanced.py 2>&1 | tee -a output_radial_profile_fit.log
"""

import argparse
import glob
import os
import shutil
import time

import numpy as np

import diskmint.execute as exe
import diskmint.model as model


def copy_input_files(input_dir, data_dir, para):
    os.makedirs(data_dir, exist_ok=True)
    patterns = [
        para.fmodel_filename,
        f"dustkappa_{para.dustopacname_1}_*.inp",
        f"dustkappa_{para.dustopacname_2}_*.inp",
        "aave.inp",
        "fracs.inp",
        "fracs_numb.inp",
        "ndsd.inp",
    ]

    copied = 0
    for pattern in dict.fromkeys(patterns):
        for source_path in glob.glob(os.path.join(input_dir, pattern)):
            if os.path.isfile(source_path):
                shutil.copy2(source_path, data_dir)
                copied += 1
    print(f"Copied {copied} input files to {data_dir}.")


def load_reference_table(path):
    table = np.loadtxt(path)
    if table.ndim != 2 or table.shape[1] < 2:
        raise ValueError(f"Expected a two-column reference table: {path}")
    return table[:, :2]


def configure_model(para, working_dir, mode):
    inputs_dir = os.path.join(working_dir, "model_inputs")
    sigmad_ref = load_reference_table(
        os.path.join(inputs_dir, "sigma_reference_template.dat")
    )

    para.edit_parameter("sigmad_ref", new_value=sigmad_ref)
    para.edit_parameter("mdiskd", new_value=4.0e-4)
    para.edit_parameter("pl_sufdens", new_value=-1.0)
    para.edit_parameter("pl_tapoff", new_value=1.0)
    para.edit_parameter("Rtap", new_value=100.0)

    if mode == "model_A_constant_g2d":
        para.edit_parameter("ratio_g2d_global", new_value=100.0)
        para.edit_parameter("chemical_save_name", new_value="example_radial_profile_fit_model_A_constant_g2d")
    elif mode == "model_B_radial_profile_fit":
        ratio_g2d_ref = load_reference_table(
            os.path.join(inputs_dir, "ratio_g2d_reference_template.dat")
        )
        para.edit_parameter("ratio_g2d_global", new_value=100.0)
        para.edit_parameter("ratio_g2d", new_value=ratio_g2d_ref)
        para.edit_parameter("chemical_save_name", new_value="example_radial_profile_fit_model_B_radial_g2d")
    else:
        raise ValueError(f"Unknown mode: {mode}")

    return para


def build_mint(para, data_dir, args):
    mint = model.Mint(para, file_dir=data_dir)
    mint.bool_MakeDustKappa = True
    mint.bool_SED = True
    mint.bool_VHSE = True
    mint.n_vhse_loop = args.n_vhse_loop
    mint.bool_dust_settling = args.dust_settling
    mint.bool_chemistry = not args.skip_chemistry
    mint.bool_savemodel = True
    mint.chem_code_dir = "reducedRGH22"
    return mint


def main():
    parser = argparse.ArgumentParser(
        description="Run the DiskMINT radial profile fitting example."
    )
    parser.add_argument(
        "--mode",
        choices=["model_A_constant_g2d", "model_B_radial_profile_fit"],
        default="model_A_constant_g2d",
        help="Choose the baseline constant-g2d model or the radial fitting template.",
    )
    parser.add_argument(
        "--n-vhse-loop",
        type=int,
        default=10,
        help="Maximum VHSE iterations.",
    )
    parser.add_argument(
        "--dust-settling",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Enable stable DiskMINT dust settling (default: enabled; use --no-dust-settling to disable).",
    )
    parser.add_argument(
        "--skip-chemistry",
        action="store_true",
        help="Skip chemistry for a faster thermal/density test run.",
    )
    args = parser.parse_args()

    time0 = time.time()
    working_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(working_dir, "data")
    input_dir = os.path.join(working_dir, "input")
    output_dir = os.path.join(working_dir, "output")

    para = model.Parameters()
    para.read_parameters_from_csv(
        filename="radial_profile_fit_parameters",
        directory=working_dir,
        extension=".csv",
    )
    para.edit_parameter("chemical_save_dir", new_value=output_dir)
    para = configure_model(para, working_dir, args.mode)

    copy_input_files(input_dir, data_dir, para)
    para.write_parameters_to_csv(
        para.chemical_save_name + "_parameters_setup",
        directory=data_dir,
        extension=".csv",
    )

    os.chdir(data_dir)
    mint = build_mint(para, data_dir, args)
    print(f"Running mode: {args.mode}")
    print(f"Model name: {para.chemical_save_name}")
    print(f"Global gas-to-dust ratio: {para.ratio_g2d_global:.1f}")
    print(f"Dust settling: {mint.bool_dust_settling}")
    print(f"Chemistry: {mint.bool_chemistry}")

    exe.runmodel(mint, test_alliteration=False)
    print(f"Finished in {(time.time() - time0) / 3600.0:.2f} hours.")


if __name__ == "__main__":
    main()
