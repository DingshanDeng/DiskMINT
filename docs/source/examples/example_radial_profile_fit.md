# Example: Fitting Continuum and C18O Radial Profiles

This example shows how to set up a DiskMINT model for iterative radial-profile fitting.
It is target-agnostic: the bundled values are templates for learning the workflow, not a published best-fit model.

The example files are located at:

```text
examples/example_diskmint_models/example_radial_profile_fit/
```

## What This Example Demonstrates

The continuum radial profile mainly constrains the dust surface density and opacity assumptions.
C18O emission also responds to gas mass, gas-to-dust ratio, temperature structure, UV field, and chemistry.

The intended workflow is:

1. Start with `model_A_constant_g2d`.
2. Compare the continuum and C18O radial profiles.
3. Adjust `model_inputs/sigma_reference_template.dat` to change the dust surface density.
4. Adjust `ratio_g2d_global`, or use `model_inputs/ratio_g2d_reference_template.dat`, to change the gas mass distribution.
5. Rerun and compare again.

## Directory Structure

```text
example_radial_profile_fit/
├── example_model_radial_profile_fit.py
├── radial_profile_fit_parameters.csv
├── compare_radial_profiles.py
├── input/
│   ├── BTSettl_1p0Msolar_1pc_um_ergpcm2hzs.inp
│   ├── dustkappa_optool_20250122_*.inp
│   ├── aave.inp
│   ├── fracs.inp
│   ├── fracs_numb.inp
│   └── ndsd.inp
├── model_inputs/
│   ├── sigma_reference_template.dat
│   └── ratio_g2d_reference_template.dat
└── target_profiles/
    ├── example_continuum_profile.csv
    └── example_c18o_profile.csv
```

## Running the Example

Run the baseline model with a constant gas-to-dust ratio:

```bash
python -u example_model_radial_profile_fit.py --mode model_A_constant_g2d 2>&1 | tee -a output_radial_profile_fit.log
```

Run the radial fitting template, which uses both a reference dust surface density and a radial gas-to-dust table:

```bash
python -u example_model_radial_profile_fit.py --mode model_B_radial_profile_fit 2>&1 | tee -a output_radial_profile_fit.log
```

For a quicker density/thermal test before running chemistry:

```bash
python -u example_model_radial_profile_fit.py --mode model_B_radial_profile_fit --skip-chemistry 2>&1 | tee -a output_radial_profile_fit.log
```

## Key Controls

The reference dust surface density is read from:

```text
model_inputs/sigma_reference_template.dat
```

It is a two-column table:

```text
radius_cm  sigma_dust_g_cm-2
```

The optional radial gas-to-dust ratio is read from:

```text
model_inputs/ratio_g2d_reference_template.dat
```

It is also a two-column table:

```text
radius_cm  gas_to_dust_ratio
```

The example uses only stable DiskMINT runtime flags: VHSE, SED, dust settling, chemistry, model saving, and dust-kappa metadata handling.
Experimental dust fragmentation, radial drift, temperature decoupling, inner-rim, and same-chemistry-grid flags are intentionally not used here.

## Plotting the Bundled Profiles

Before or after a model run, plot the template profiles with:

```bash
python compare_radial_profiles.py
```

The helper writes:

```text
profile_comparison/radial_profile_fit_inputs.png
```

Extend this plotting helper with your own continuum image or C18O line-radiative-transfer products for direct model-vs-observation comparisons.
