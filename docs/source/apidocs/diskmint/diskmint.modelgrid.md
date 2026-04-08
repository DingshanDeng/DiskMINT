# {py:mod}`diskmint.modelgrid`

```{py:module} diskmint.modelgrid
```

```{autodoc2-docstring} diskmint.modelgrid
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`radmc3dGrid <diskmint.modelgrid.radmc3dGrid>`
  - ```{autodoc2-docstring} diskmint.modelgrid.radmc3dGrid
    :summary:
    ```
* - {py:obj}`radmc3dData <diskmint.modelgrid.radmc3dData>`
  - ```{autodoc2-docstring} diskmint.modelgrid.radmc3dData
    :summary:
    ```
````

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`readGrid <diskmint.modelgrid.readGrid>`
  - ```{autodoc2-docstring} diskmint.modelgrid.readGrid
    :summary:
    ```
* - {py:obj}`readData <diskmint.modelgrid.readData>`
  - ```{autodoc2-docstring} diskmint.modelgrid.readData
    :summary:
    ```
````

### API

`````{py:class} radmc3dGrid()
:canonical: diskmint.modelgrid.radmc3dGrid

Bases: {py:obj}`object`

```{autodoc2-docstring} diskmint.modelgrid.radmc3dGrid
```

```{rubric} Initialization
```

```{autodoc2-docstring} diskmint.modelgrid.radmc3dGrid.__init__
```

````{py:method} makeWavelengthGrid(wbound=None, nw=None, ppar=None)
:canonical: diskmint.modelgrid.radmc3dGrid.makeWavelengthGrid

```{autodoc2-docstring} diskmint.modelgrid.radmc3dGrid.makeWavelengthGrid
```

````

````{py:method} writeWavelengthGrid(fname='', old=False)
:canonical: diskmint.modelgrid.radmc3dGrid.writeWavelengthGrid

```{autodoc2-docstring} diskmint.modelgrid.radmc3dGrid.writeWavelengthGrid
```

````

````{py:method} makeSpatialGrid(crd_sys=None, xbound=None, ybound=None, zbound=None, nxi=None, nyi=None, nzi=None, ppar=None)
:canonical: diskmint.modelgrid.radmc3dGrid.makeSpatialGrid

```{autodoc2-docstring} diskmint.modelgrid.radmc3dGrid.makeSpatialGrid
```

````

````{py:method} writeSpatialGrid(fname='', old=False)
:canonical: diskmint.modelgrid.radmc3dGrid.writeSpatialGrid

```{autodoc2-docstring} diskmint.modelgrid.radmc3dGrid.writeSpatialGrid
```

````

````{py:method} readWavelengthGrid(fname=None, old=False)
:canonical: diskmint.modelgrid.radmc3dGrid.readWavelengthGrid

```{autodoc2-docstring} diskmint.modelgrid.radmc3dGrid.readWavelengthGrid
```

````

````{py:method} readSpatialGrid(fname='', old=False)
:canonical: diskmint.modelgrid.radmc3dGrid.readSpatialGrid

```{autodoc2-docstring} diskmint.modelgrid.radmc3dGrid.readSpatialGrid
```

````

````{py:method} readGrid(old=False)
:canonical: diskmint.modelgrid.radmc3dGrid.readGrid

```{autodoc2-docstring} diskmint.modelgrid.radmc3dGrid.readGrid
```

````

````{py:method} getCellVolume()
:canonical: diskmint.modelgrid.radmc3dGrid.getCellVolume

```{autodoc2-docstring} diskmint.modelgrid.radmc3dGrid.getCellVolume
```

````

`````

````{py:function} readGrid(sgrid=True, wgrid=True, sgrid_fname=None, wgrid_fname=None, old=False)
:canonical: diskmint.modelgrid.readGrid

```{autodoc2-docstring} diskmint.modelgrid.readGrid
```
````

`````{py:class} radmc3dData(grid=None)
:canonical: diskmint.modelgrid.radmc3dData

Bases: {py:obj}`object`

```{autodoc2-docstring} diskmint.modelgrid.radmc3dData
```

```{rubric} Initialization
```

```{autodoc2-docstring} diskmint.modelgrid.radmc3dData.__init__
```

````{py:method} getTauOneDust(idust=0, axis='', kappa=0.0)
:canonical: diskmint.modelgrid.radmc3dData.getTauOneDust

```{autodoc2-docstring} diskmint.modelgrid.radmc3dData.getTauOneDust
```

````

````{py:method} getTau(idust=None, axis='xy', wav=0.0, kappa=None, old=False)
:canonical: diskmint.modelgrid.radmc3dData.getTau

```{autodoc2-docstring} diskmint.modelgrid.radmc3dData.getTau
```

````

````{py:method} readDustDens(fname=None, old=False)
:canonical: diskmint.modelgrid.radmc3dData.readDustDens

```{autodoc2-docstring} diskmint.modelgrid.radmc3dData.readDustDens
```

````

````{py:method} readDustTemp(fname=None, old=False)
:canonical: diskmint.modelgrid.radmc3dData.readDustTemp

```{autodoc2-docstring} diskmint.modelgrid.radmc3dData.readDustTemp
```

````

````{py:method} readGasVel(fname=None)
:canonical: diskmint.modelgrid.radmc3dData.readGasVel

```{autodoc2-docstring} diskmint.modelgrid.radmc3dData.readGasVel
```

````

````{py:method} readVTurb(fname=None)
:canonical: diskmint.modelgrid.radmc3dData.readVTurb

```{autodoc2-docstring} diskmint.modelgrid.radmc3dData.readVTurb
```

````

````{py:method} readGasDens(fname=None, ispec='')
:canonical: diskmint.modelgrid.radmc3dData.readGasDens

```{autodoc2-docstring} diskmint.modelgrid.radmc3dData.readGasDens
```

````

````{py:method} readGasTemp(fname=None)
:canonical: diskmint.modelgrid.radmc3dData.readGasTemp

```{autodoc2-docstring} diskmint.modelgrid.radmc3dData.readGasTemp
```

````

````{py:method} writeDustDens(fname='', binary=True, old=False)
:canonical: diskmint.modelgrid.radmc3dData.writeDustDens

```{autodoc2-docstring} diskmint.modelgrid.radmc3dData.writeDustDens
```

````

````{py:method} writeDustTemp(fname='', binary=True)
:canonical: diskmint.modelgrid.radmc3dData.writeDustTemp

```{autodoc2-docstring} diskmint.modelgrid.radmc3dData.writeDustTemp
```

````

````{py:method} writeGasDens(fname=None, ispec='', binary=True)
:canonical: diskmint.modelgrid.radmc3dData.writeGasDens

```{autodoc2-docstring} diskmint.modelgrid.radmc3dData.writeGasDens
```

````

````{py:method} writeGasTemp(fname='', binary=True)
:canonical: diskmint.modelgrid.radmc3dData.writeGasTemp

```{autodoc2-docstring} diskmint.modelgrid.radmc3dData.writeGasTemp
```

````

````{py:method} writeGasVel(fname='', binary=True)
:canonical: diskmint.modelgrid.radmc3dData.writeGasVel

```{autodoc2-docstring} diskmint.modelgrid.radmc3dData.writeGasVel
```

````

````{py:method} writeVTurb(fname='', binary=True)
:canonical: diskmint.modelgrid.radmc3dData.writeVTurb

```{autodoc2-docstring} diskmint.modelgrid.radmc3dData.writeVTurb
```

````

````{py:method} getSigmaDust(idust=-1)
:canonical: diskmint.modelgrid.radmc3dData.getSigmaDust

```{autodoc2-docstring} diskmint.modelgrid.radmc3dData.getSigmaDust
```

````

````{py:method} getSigmaGas()
:canonical: diskmint.modelgrid.radmc3dData.getSigmaGas

```{autodoc2-docstring} diskmint.modelgrid.radmc3dData.getSigmaGas
```

````

`````

````{py:function} readData(ddens=False, dtemp=False, gdens=False, gtemp=False, gvel=False, ispec=None, vturb=False, grid=None, old=False, binary=None)
:canonical: diskmint.modelgrid.readData

```{autodoc2-docstring} diskmint.modelgrid.readData
```
````
