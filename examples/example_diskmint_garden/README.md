# DiskMINT-GARDEN Quickstart

This example shows the `diskmint.garden` inference API. It does not run a
full DiskMINT thermochemical model; it loads the bundled DiskMINT-GARDEN
surrogate model and predicts disk properties from observed fluxes.

Install the optional dependencies first:

```bash
pip install "diskmint[garden]"
```

Run the smoke example:

```bash
python examples/example_diskmint_garden/run_garden_example.py
```

Or open `diskmint_garden_quickstart.ipynb` for an annotated notebook version.
