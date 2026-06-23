# DiskMINT-GARDEN Quickstart

DiskMINT-GARDEN is the beta inference API for estimating disk properties from
observed millimeter continuum and $\mathrm{C^{18}O}$ fluxes. Unlike the other
examples in this section, it does not run RADMC-3D or the chemistry network.

Install the optional dependencies:

```bash
pip install "diskmint[garden]"
```

Run the repository smoke example:

```bash
python examples/example_diskmint_garden/run_garden_example.py
```

The same workflow is available as an annotated notebook:

```text
examples/example_diskmint_garden/diskmint_garden_quickstart.ipynb
```

Minimal single-target use:

```python
import diskmint.garden.infer as garden

result = garden.from_observations(
    flux_mm=120.0,       # continuum flux density [mJy]
    flux_c18o=850.0,     # C18O integrated flux [mJy km/s]
    distance=140.0,      # distance [pc]
    mstar=0.8,           # stellar mass [Msun]
    rdust_90=80.0,       # 90 percent dust radius [au]
    band="band6",        # band6: C18O J=2-1; band7: C18O J=3-2
)

print(result["Mgas_pred_Msun"])
print(result["gtd_pred_dimless"])
print(result["Rc_pred_au"])
```

Always check `is_outside_grid` and `outside_features` in the output. A flagged
prediction means the observed inputs fall outside the model training domain.
