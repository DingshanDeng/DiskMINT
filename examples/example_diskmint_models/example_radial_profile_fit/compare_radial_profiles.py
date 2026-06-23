#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Lightweight plotting helper for the radial profile fitting example.

This script plots the bundled target profiles and the input fitting controls.
After a full model run, extend this file with project-specific image or line
radiative transfer products for direct model-vs-observation comparisons.
"""

import os

import matplotlib.pyplot as plt
import numpy as np


def read_csv(path):
    return np.genfromtxt(path, delimiter=",", names=True)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    model_inputs = os.path.join(here, "model_inputs")
    target_profiles = os.path.join(here, "target_profiles")
    output_dir = os.path.join(here, "profile_comparison")
    os.makedirs(output_dir, exist_ok=True)

    sigmad = np.loadtxt(os.path.join(model_inputs, "sigma_reference_template.dat"))
    ratio_g2d = np.loadtxt(os.path.join(model_inputs, "ratio_g2d_reference_template.dat"))
    cont = read_csv(os.path.join(target_profiles, "example_continuum_profile.csv"))
    c18o = read_csv(os.path.join(target_profiles, "example_c18o_profile.csv"))

    au = 1.495978707e13
    radius_sigmad_au = sigmad[:, 0] / au
    radius_g2d_au = ratio_g2d[:, 0] / au

    fig, axes = plt.subplots(2, 2, figsize=(8, 6), constrained_layout=True)

    axes[0, 0].plot(radius_sigmad_au, sigmad[:, 1], color="tab:blue")
    axes[0, 0].set_xscale("log")
    axes[0, 0].set_yscale("log")
    axes[0, 0].set_xlabel("Radius [au]")
    axes[0, 0].set_ylabel("Dust surface density [g cm$^{-2}$]")

    axes[0, 1].plot(radius_g2d_au, ratio_g2d[:, 1], color="tab:orange")
    axes[0, 1].set_xlabel("Radius [au]")
    axes[0, 1].set_ylabel("Gas-to-dust ratio")

    axes[1, 0].errorbar(
        cont["radius_au"],
        cont["intensity_norm"],
        yerr=cont["uncertainty_norm"],
        marker="o",
        color="tab:green",
    )
    axes[1, 0].set_xlabel("Radius [au]")
    axes[1, 0].set_ylabel("Continuum intensity [normalized]")

    axes[1, 1].errorbar(
        c18o["radius_au"],
        c18o["integrated_intensity_norm"],
        yerr=c18o["uncertainty_norm"],
        marker="o",
        color="tab:red",
    )
    axes[1, 1].set_xlabel("Radius [au]")
    axes[1, 1].set_ylabel("C18O integrated intensity [normalized]")

    fig.savefig(os.path.join(output_dir, "radial_profile_fit_inputs.png"), dpi=200)
    print(f"Saved {os.path.join(output_dir, 'radial_profile_fit_inputs.png')}")


if __name__ == "__main__":
    main()
