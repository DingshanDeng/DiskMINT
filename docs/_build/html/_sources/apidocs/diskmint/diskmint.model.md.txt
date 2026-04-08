# {py:mod}`diskmint.model`

```{py:module} diskmint.model
```

```{autodoc2-docstring} diskmint.model
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`Parameters <diskmint.model.Parameters>`
  - ```{autodoc2-docstring} diskmint.model.Parameters
    :summary:
    ```
* - {py:obj}`Mint <diskmint.model.Mint>`
  - ```{autodoc2-docstring} diskmint.model.Mint
    :summary:
    ```
````

### API

`````{py:class} Parameters()
:canonical: diskmint.model.Parameters

```{autodoc2-docstring} diskmint.model.Parameters
```

```{rubric} Initialization
```

```{autodoc2-docstring} diskmint.model.Parameters.__init__
```

````{py:method} set_parameter(name, value, units=None, description=None, valuetype='float64')
:canonical: diskmint.model.Parameters.set_parameter

```{autodoc2-docstring} diskmint.model.Parameters.set_parameter
```

````

````{py:method} edit_parameter(name, new_value=None, new_units=None, new_description=None, new_valuetype=None)
:canonical: diskmint.model.Parameters.edit_parameter

```{autodoc2-docstring} diskmint.model.Parameters.edit_parameter
```

````

````{py:method} get_parameter_info(name)
:canonical: diskmint.model.Parameters.get_parameter_info

```{autodoc2-docstring} diskmint.model.Parameters.get_parameter_info
```

````

````{py:method} get_parameter_value(name)
:canonical: diskmint.model.Parameters.get_parameter_value

```{autodoc2-docstring} diskmint.model.Parameters.get_parameter_value
```

````

````{py:method} read_parameters_from_csv(filename, directory='./', extension='', verbose=False)
:canonical: diskmint.model.Parameters.read_parameters_from_csv

```{autodoc2-docstring} diskmint.model.Parameters.read_parameters_from_csv
```

````

````{py:method} write_parameters_to_csv(filename, directory='./', extension='')
:canonical: diskmint.model.Parameters.write_parameters_to_csv

```{autodoc2-docstring} diskmint.model.Parameters.write_parameters_to_csv
```

````

````{py:property} nphot
:canonical: diskmint.model.Parameters.nphot

```{autodoc2-docstring} diskmint.model.Parameters.nphot
```

````

````{py:property} nthreads
:canonical: diskmint.model.Parameters.nthreads

```{autodoc2-docstring} diskmint.model.Parameters.nthreads
```

````

````{py:property} scattering_mode
:canonical: diskmint.model.Parameters.scattering_mode

```{autodoc2-docstring} diskmint.model.Parameters.scattering_mode
```

````

````{py:property} nr
:canonical: diskmint.model.Parameters.nr

```{autodoc2-docstring} diskmint.model.Parameters.nr
```

````

````{py:property} ntheta
:canonical: diskmint.model.Parameters.ntheta

```{autodoc2-docstring} diskmint.model.Parameters.ntheta
```

````

````{py:property} nphi
:canonical: diskmint.model.Parameters.nphi

```{autodoc2-docstring} diskmint.model.Parameters.nphi
```

````

````{py:property} rin
:canonical: diskmint.model.Parameters.rin

```{autodoc2-docstring} diskmint.model.Parameters.rin
```

````

````{py:property} rout
:canonical: diskmint.model.Parameters.rout

```{autodoc2-docstring} diskmint.model.Parameters.rout
```

````

````{py:property} thetaup
:canonical: diskmint.model.Parameters.thetaup

```{autodoc2-docstring} diskmint.model.Parameters.thetaup
```

````

````{py:property} r_step_max
:canonical: diskmint.model.Parameters.r_step_max

```{autodoc2-docstring} diskmint.model.Parameters.r_step_max
```

````

````{py:property} fmodel_filename
:canonical: diskmint.model.Parameters.fmodel_filename

```{autodoc2-docstring} diskmint.model.Parameters.fmodel_filename
```

````

````{py:property} tstar
:canonical: diskmint.model.Parameters.tstar

```{autodoc2-docstring} diskmint.model.Parameters.tstar
```

````

````{py:property} mstar
:canonical: diskmint.model.Parameters.mstar

```{autodoc2-docstring} diskmint.model.Parameters.mstar
```

````

````{py:property} rstar
:canonical: diskmint.model.Parameters.rstar

```{autodoc2-docstring} diskmint.model.Parameters.rstar
```

````

````{py:property} incl
:canonical: diskmint.model.Parameters.incl

```{autodoc2-docstring} diskmint.model.Parameters.incl
```

````

````{py:property} phi
:canonical: diskmint.model.Parameters.phi

```{autodoc2-docstring} diskmint.model.Parameters.phi
```

````

````{py:property} mdotacc
:canonical: diskmint.model.Parameters.mdotacc

```{autodoc2-docstring} diskmint.model.Parameters.mdotacc
```

````

````{py:property} mdiskd
:canonical: diskmint.model.Parameters.mdiskd

```{autodoc2-docstring} diskmint.model.Parameters.mdiskd
```

````

````{py:property} mdiskd_2
:canonical: diskmint.model.Parameters.mdiskd_2

```{autodoc2-docstring} diskmint.model.Parameters.mdiskd_2
```

````

````{py:property} ratio_g2d_global
:canonical: diskmint.model.Parameters.ratio_g2d_global

```{autodoc2-docstring} diskmint.model.Parameters.ratio_g2d_global
```

````

````{py:property} ratio_g2d
:canonical: diskmint.model.Parameters.ratio_g2d

```{autodoc2-docstring} diskmint.model.Parameters.ratio_g2d
```

````

````{py:property} sigmad_ref
:canonical: diskmint.model.Parameters.sigmad_ref

```{autodoc2-docstring} diskmint.model.Parameters.sigmad_ref
```

````

````{py:property} pl_sufdens
:canonical: diskmint.model.Parameters.pl_sufdens

```{autodoc2-docstring} diskmint.model.Parameters.pl_sufdens
```

````

````{py:property} pl_sufdens_2
:canonical: diskmint.model.Parameters.pl_sufdens_2

```{autodoc2-docstring} diskmint.model.Parameters.pl_sufdens_2
```

````

````{py:property} pl_tapoff
:canonical: diskmint.model.Parameters.pl_tapoff

```{autodoc2-docstring} diskmint.model.Parameters.pl_tapoff
```

````

````{py:property} rdisk_in
:canonical: diskmint.model.Parameters.rdisk_in

```{autodoc2-docstring} diskmint.model.Parameters.rdisk_in
```

````

````{py:property} Rtap
:canonical: diskmint.model.Parameters.Rtap

```{autodoc2-docstring} diskmint.model.Parameters.Rtap
```

````

````{py:property} rdisk_out_2
:canonical: diskmint.model.Parameters.rdisk_out_2

```{autodoc2-docstring} diskmint.model.Parameters.rdisk_out_2
```

````

````{py:property} scaleheight_index
:canonical: diskmint.model.Parameters.scaleheight_index

```{autodoc2-docstring} diskmint.model.Parameters.scaleheight_index
```

````

````{py:property} scaleheight_index_2
:canonical: diskmint.model.Parameters.scaleheight_index_2

```{autodoc2-docstring} diskmint.model.Parameters.scaleheight_index_2
```

````

````{py:property} hp100
:canonical: diskmint.model.Parameters.hp100

```{autodoc2-docstring} diskmint.model.Parameters.hp100
```

````

````{py:property} hp100_2
:canonical: diskmint.model.Parameters.hp100_2

```{autodoc2-docstring} diskmint.model.Parameters.hp100_2
```

````

````{py:property} hprcr
:canonical: diskmint.model.Parameters.hprcr

```{autodoc2-docstring} diskmint.model.Parameters.hprcr
```

````

````{py:property} plh
:canonical: diskmint.model.Parameters.plh

```{autodoc2-docstring} diskmint.model.Parameters.plh
```

````

````{py:property} plh_2
:canonical: diskmint.model.Parameters.plh_2

```{autodoc2-docstring} diskmint.model.Parameters.plh_2
```

````

````{py:property} dustopacname_1
:canonical: diskmint.model.Parameters.dustopacname_1

```{autodoc2-docstring} diskmint.model.Parameters.dustopacname_1
```

````

````{py:property} dustopacname_2
:canonical: diskmint.model.Parameters.dustopacname_2

```{autodoc2-docstring} diskmint.model.Parameters.dustopacname_2
```

````

````{py:property} nr_dust_1
:canonical: diskmint.model.Parameters.nr_dust_1

```{autodoc2-docstring} diskmint.model.Parameters.nr_dust_1
```

````

````{py:property} nr_dust_2
:canonical: diskmint.model.Parameters.nr_dust_2

```{autodoc2-docstring} diskmint.model.Parameters.nr_dust_2
```

````

````{py:property} dust_spec_nr
:canonical: diskmint.model.Parameters.dust_spec_nr

```{autodoc2-docstring} diskmint.model.Parameters.dust_spec_nr
```

````

````{py:property} amin_1
:canonical: diskmint.model.Parameters.amin_1

```{autodoc2-docstring} diskmint.model.Parameters.amin_1
```

````

````{py:property} amin_2
:canonical: diskmint.model.Parameters.amin_2

```{autodoc2-docstring} diskmint.model.Parameters.amin_2
```

````

````{py:property} amax_1
:canonical: diskmint.model.Parameters.amax_1

```{autodoc2-docstring} diskmint.model.Parameters.amax_1
```

````

````{py:property} amax_2
:canonical: diskmint.model.Parameters.amax_2

```{autodoc2-docstring} diskmint.model.Parameters.amax_2
```

````

````{py:property} amin_all
:canonical: diskmint.model.Parameters.amin_all

```{autodoc2-docstring} diskmint.model.Parameters.amin_all
```

````

````{py:property} amax_all
:canonical: diskmint.model.Parameters.amax_all

```{autodoc2-docstring} diskmint.model.Parameters.amax_all
```

````

````{py:property} pla_dustsize
:canonical: diskmint.model.Parameters.pla_dustsize

```{autodoc2-docstring} diskmint.model.Parameters.pla_dustsize
```

````

````{py:property} rhobulk
:canonical: diskmint.model.Parameters.rhobulk

```{autodoc2-docstring} diskmint.model.Parameters.rhobulk
```

````

````{py:property} visc_alpha
:canonical: diskmint.model.Parameters.visc_alpha

```{autodoc2-docstring} diskmint.model.Parameters.visc_alpha
```

````

````{py:property} vel_frag
:canonical: diskmint.model.Parameters.vel_frag

```{autodoc2-docstring} diskmint.model.Parameters.vel_frag
```

````

````{py:property} radius_drift
:canonical: diskmint.model.Parameters.radius_drift

```{autodoc2-docstring} diskmint.model.Parameters.radius_drift
```

````

````{py:property} a_drift
:canonical: diskmint.model.Parameters.a_drift

```{autodoc2-docstring} diskmint.model.Parameters.a_drift
```

````

````{py:property} r_wall_surf
:canonical: diskmint.model.Parameters.r_wall_surf

```{autodoc2-docstring} diskmint.model.Parameters.r_wall_surf
```

````

````{py:property} factor_wall_height
:canonical: diskmint.model.Parameters.factor_wall_height

```{autodoc2-docstring} diskmint.model.Parameters.factor_wall_height
```

````

````{py:property} chemical_save_dir
:canonical: diskmint.model.Parameters.chemical_save_dir

```{autodoc2-docstring} diskmint.model.Parameters.chemical_save_dir
```

````

````{py:property} chemical_save_name
:canonical: diskmint.model.Parameters.chemical_save_name

```{autodoc2-docstring} diskmint.model.Parameters.chemical_save_name
```

````

````{py:property} nr_cyl_LIME
:canonical: diskmint.model.Parameters.nr_cyl_LIME

```{autodoc2-docstring} diskmint.model.Parameters.nr_cyl_LIME
```

````

````{py:property} nz_cyl_LIME
:canonical: diskmint.model.Parameters.nz_cyl_LIME

```{autodoc2-docstring} diskmint.model.Parameters.nz_cyl_LIME
```

````

````{py:property} R_temp_trans
:canonical: diskmint.model.Parameters.R_temp_trans

```{autodoc2-docstring} diskmint.model.Parameters.R_temp_trans
```

````

````{py:property} fact_Tgas_2_Tdust
:canonical: diskmint.model.Parameters.fact_Tgas_2_Tdust

```{autodoc2-docstring} diskmint.model.Parameters.fact_Tgas_2_Tdust
```

````

````{py:property} G0Hab_set
:canonical: diskmint.model.Parameters.G0Hab_set

```{autodoc2-docstring} diskmint.model.Parameters.G0Hab_set
```

````

````{py:method} set_default_parameters()
:canonical: diskmint.model.Parameters.set_default_parameters

```{autodoc2-docstring} diskmint.model.Parameters.set_default_parameters
```

````

`````

`````{py:class} Mint(parameters, file_dir=os.getcwd())
:canonical: diskmint.model.Mint

```{autodoc2-docstring} diskmint.model.Mint
```

```{rubric} Initialization
```

```{autodoc2-docstring} diskmint.model.Mint.__init__
```

````{py:attribute} bool_savemodel
:canonical: diskmint.model.Mint.bool_savemodel
:value: >
   True

```{autodoc2-docstring} diskmint.model.Mint.bool_savemodel
```

````

````{py:attribute} bool_dust_inner_rim
:canonical: diskmint.model.Mint.bool_dust_inner_rim
:value: >
   False

```{autodoc2-docstring} diskmint.model.Mint.bool_dust_inner_rim
```

````

````{py:method} setup_radmc3d_grid(para)
:canonical: diskmint.model.Mint.setup_radmc3d_grid

```{autodoc2-docstring} diskmint.model.Mint.setup_radmc3d_grid
```

````

````{py:method} calculate_volume()
:canonical: diskmint.model.Mint.calculate_volume

```{autodoc2-docstring} diskmint.model.Mint.calculate_volume
```

````

````{py:method} setup_dust_density_initial_gaussian(para, bool_initial_setup=False)
:canonical: diskmint.model.Mint.setup_dust_density_initial_gaussian

```{autodoc2-docstring} diskmint.model.Mint.setup_dust_density_initial_gaussian
```

````

````{py:method} setup_secondary_dust_density_initial_gaussian(para, bool_initial_setup=False)
:canonical: diskmint.model.Mint.setup_secondary_dust_density_initial_gaussian

```{autodoc2-docstring} diskmint.model.Mint.setup_secondary_dust_density_initial_gaussian
```

````

````{py:method} setup_g2d_grid(para, bool_initial_setup=False)
:canonical: diskmint.model.Mint.setup_g2d_grid

```{autodoc2-docstring} diskmint.model.Mint.setup_g2d_grid
```

````

````{py:method} setup_gas_density(para, bool_initial_setup=False)
:canonical: diskmint.model.Mint.setup_gas_density

```{autodoc2-docstring} diskmint.model.Mint.setup_gas_density
```

````

````{py:method} setup_wavelengths(para)
:canonical: diskmint.model.Mint.setup_wavelengths

```{autodoc2-docstring} diskmint.model.Mint.setup_wavelengths
```

````

````{py:method} setup_dust_info(para, bool_savefile=False)
:canonical: diskmint.model.Mint.setup_dust_info

```{autodoc2-docstring} diskmint.model.Mint.setup_dust_info
```

````

````{py:method} readin_dust_info(file_dir, file_aave='aave.inp', file_fracs='fracs.inp', file_fracs_numb='fracs_numb.inp', file_ndsd='ndsd.inp')
:canonical: diskmint.model.Mint.readin_dust_info

```{autodoc2-docstring} diskmint.model.Mint.readin_dust_info
```

````

`````
