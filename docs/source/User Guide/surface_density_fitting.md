# Surface Density Fitting Guide

This page summarizes how DiskMINT can be used to infer radial surface density profiles from resolved observations. A fully reproducible worked example will be added in a future update.

---

DiskMINT supports a direct profile-fitting workflow in addition to analytic parameterized surface density models. In this mode, the dust and gas surface densities are treated separately and are iteratively adjusted to match resolved observations. The approach was developed and applied in [Deng et al. (2025), ApJ 995, 98](https://ui.adsabs.harvard.edu/abs/2025ApJ...995...98D) for IM Lup.

## Goal

The goal is to derive:

- $\Sigma_\mathrm{dust}(r)$ from the resolved dust continuum radial profile
- $\Sigma_\mathrm{gas}(r)$ from the resolved C$^{18}$O radial profile

This allows the dust-to-gas ratio to vary with radius instead of being fixed to a single global value.

## Basic assumptions

This workflow assumes:

- the outer-disk dust continuum is approximately optically thin, so intensity per area roughly traces dust surface density
- the C$^{18}$O radial profile can be interpreted with the thermochemical model as a tracer of gas surface density
- the disk inclination and position angle are known well enough to extract radial profiles in deprojected elliptical annuli

As with any profile-based inference, the least constrained regions are the optically thick inner disk and the outermost disk beyond the measured gas or dust radii.

## High-level workflow

1. Start from an initial analytic model for the disk structure.
2. Generate synthetic continuum and line images with DiskMINT.
3. Measure radial profiles from both the synthetic images and the observations using elliptical annuli matched to the observed geometry.
4. Update $\Sigma_\mathrm{dust}(r)$ using the ratio between observed and modeled continuum radial intensity.
5. Update $\Sigma_\mathrm{gas}(r)$ using the ratio between observed and modeled C$^{18}$O radial intensity.
6. Re-run DiskMINT and repeat until the radial profiles converge.

This procedure allows the dust and gas radial structures to be fit separately, which is important for disks where the dust is more compact than the gas.

## Example

*This section is under construction*

The complete worked example is planned for a future documentation update and will include:

- expected observational inputs
- how the radial profiles are extracted
- how custom surface density profiles are passed into the modeling workflow
- which scripts or notebooks run the fitting loop
- common caveats such as optical depth, beam smearing, and unconstrained inner or outer radii
  
The IM Lup analysis in Deng et al. (2025) is the primary reference for the method.
