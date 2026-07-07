# {py:mod}`diskmint.garden.infer`

```{py:module} diskmint.garden.infer
```

```{autodoc2-docstring} diskmint.garden.infer
:allowtitles:
```

## Module Contents

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`get_model_path <diskmint.garden.infer.get_model_path>`
  - ```{autodoc2-docstring} diskmint.garden.infer.get_model_path
    :summary:
    ```
* - {py:obj}`load_model <diskmint.garden.infer.load_model>`
  - ```{autodoc2-docstring} diskmint.garden.infer.load_model
    :summary:
    ```
* - {py:obj}`from_dataframe <diskmint.garden.infer.from_dataframe>`
  - ```{autodoc2-docstring} diskmint.garden.infer.from_dataframe
    :summary:
    ```
* - {py:obj}`from_observations <diskmint.garden.infer.from_observations>`
  - ```{autodoc2-docstring} diskmint.garden.infer.from_observations
    :summary:
    ```
````

### API

`````{py:exception} GardenDependencyError()
:canonical: diskmint.garden.infer.GardenDependencyError

Bases: {py:obj}`ImportError`

```{autodoc2-docstring} diskmint.garden.infer.GardenDependencyError
```

```{rubric} Initialization
```

```{autodoc2-docstring} diskmint.garden.infer.GardenDependencyError.__init__
```

````{py:method} add_note()
:canonical: diskmint.garden.infer.GardenDependencyError.add_note

````

```{py:class} args
:canonical: diskmint.garden.infer.GardenDependencyError.args

```

```{py:class} msg
:canonical: diskmint.garden.infer.GardenDependencyError.msg

```

```{py:class} name
:canonical: diskmint.garden.infer.GardenDependencyError.name

```

```{py:class} path
:canonical: diskmint.garden.infer.GardenDependencyError.path

```

````{py:method} with_traceback()
:canonical: diskmint.garden.infer.GardenDependencyError.with_traceback

````

`````

````{py:function} get_model_path(band: str = 'band6')
:canonical: diskmint.garden.infer.get_model_path

```{autodoc2-docstring} diskmint.garden.infer.get_model_path
```
````

````{py:function} load_model(band: str = 'band6') -> dict[str, typing.Any]
:canonical: diskmint.garden.infer.load_model

```{autodoc2-docstring} diskmint.garden.infer.load_model
```
````

````{py:function} from_dataframe(observations, *, band: str = 'band6', model: dict[str, typing.Any] | None = None, continuum_flux_col: str = 'flux_mm_mjy', c18o_flux_col: str = 'flux_c18o_mjy_kms', distance_col: str = 'distance_pc', mstar_col: str = 'mstar_msun', rdust_col: str = 'Rdust_90_au', wavelength=None, frequency=None, rest_frequency_ghz: float | None = None)
:canonical: diskmint.garden.infer.from_dataframe

```{autodoc2-docstring} diskmint.garden.infer.from_dataframe
```
````

````{py:function} from_observations(*, flux_mm, flux_c18o, distance, mstar, rdust_90, band: str = 'band6', wavelength=None, frequency=None, rest_frequency_ghz: float | None = None) -> dict[str, typing.Any]
:canonical: diskmint.garden.infer.from_observations

```{autodoc2-docstring} diskmint.garden.infer.from_observations
```
````
