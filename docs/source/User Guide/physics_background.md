# Physics Background

A brief summary of the physical models underlying DiskMINT. Each section points to the relevant parameters and the primary literature for depth.

---

## Vertical Hydrostatic Equilibrium (VHSE)

DiskMINT solves for a self-consistent 2D disk structure by iterating between two steps:

1. **Radiative transfer** — RADMC-3D computes the dust (and gas) temperature by propagating stellar photons through the density structure via Monte Carlo (controlled by `nphot`, `nthreads`).
2. **Vertical hydrostatic equilibrium** — given the temperature, the vertical pressure gradient is balanced against gravity to update the density structure:

   $$\frac{\partial P}{\partial z} = -\rho_\mathrm{gas} \frac{GM_\star z}{(r^2 + z^2)^{3/2}}$$

   where $P = \rho_\mathrm{gas} c_s^2$ with the local isothermal sound speed $c_s = \sqrt{k_B T / \mu m_H}$.

The loop runs until the relative change in gas density between iterations falls below 5 % (`delta_rho/rho_max < 0.05`), or until `n_vhse_loop` iterations are reached. Roughly 10–20 iterations are required to converge from a Gaussian initial condition.

**Relevant parameters:** `n_vhse_loop`, `nphot`, `nthreads`, `nr`, `ntheta`

**Primary reference:** [Deng et al. (2025), ApJ 995, 98](https://ui.adsabs.harvard.edu/abs/2025ApJ...995...98D)

---

## Dust Settling

When `bool_dust_settling = True`, DiskMINT applies the Dubrulle prescription for turbulent dust settling within each VHSE iteration. The dust scale height $h_d$ for a grain of size $a$ is smaller than the gas scale height $h_g$ by:

$$\frac{h_d}{h_g} = \left(1 + \frac{\mathrm{St}}{\alpha_v}\right)^{-1/2}$$

where St is the Stokes number and $\alpha_v$ is the turbulence parameter (`visc_alpha`). Larger grains settle closer to the midplane; smaller grains remain well-mixed with the gas.

This produces a multi-species density structure where each size bin has its own vertical profile — important for correctly computing the dust temperature in an optically thick midplane.

**Relevant parameters:** `visc_alpha`, `nr_dust_1`, `amin_1`, `amax_1`, `pla_dustsize`

**Reference:** Dubrulle et al. (1995), A&A 309, 209; see also [Estrada et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016ApJ...818..200E)

---

## Grain Size Distribution

DiskMINT models the grain population as a power law in grain radius $a$:

$$n(a) \propto a^{-p}$$

for $a_\mathrm{min} \le a \le a_\mathrm{max}$, discretised into `nr_dust_1` logarithmically spaced size bins. The MRN distribution (Mathis, Rumpl & Nordsieck 1977) uses $p = 3.5$ (`pla_dustsize`). Each bin has its own opacity table (`dustkappa_{name}_{i}.inp`).

The code computes cross-section-weighted averages (`ndsd.inp`), mass fractions (`fracs.inp`), and number fractions (`fracs_numb.inp`) automatically when `bool_MakeDustKappa = True`.

**Relevant parameters:** `amin_1`, `amax_1`, `nr_dust_1`, `pla_dustsize`, `rhobulk`

---

## Surface Density Profile

DiskMINT uses a tapered power law for the surface density:

$$\Sigma(r) \propto \left(\frac{r}{R_\mathrm{tap}}\right)^{\gamma} \exp\!\left[-\left(\frac{r}{R_\mathrm{tap}}\right)^{\xi}\right]$$

where $\gamma =$ `pl_sufdens` and $\xi =$ `pl_tapoff`. This is the self-similar solution from viscous disk evolution (Lynden-Bell & Pringle 1974), which gives $\xi = 2 + \gamma$ (i.e., `pl_tapoff = 2 + pl_sufdens` for a physically motivated model). The parameter `Rtap` is the characteristic radius — roughly the outer edge of the main disk body.

**Relevant parameters:** `pl_sufdens`, `pl_tapoff`, `Rtap`, `mdiskd`, `ratio_g2d_global`

---

## Chemical Network

The chemical network is a reduced version of the full network of [Ruaud, Gorti & Hollenbach (2022)](https://ui.adsabs.harvard.edu/abs/2022ApJ...925...49R/abstract), focusing on the processes that control C¹⁸O abundance:

1. **Isotopologue-selective photodissociation** — UV photons from the star preferentially dissociate ¹²CO and ¹³CO over C¹⁸O due to self-shielding differences. This is the dominant mechanism depleting CO in the upper disk layers.
2. **Grain-surface chemistry** — CO freezes onto grain surfaces at $T < 20$–30 K and converts to CO₂ ice via grain-surface reactions. This depletes gas-phase CO in the cold midplane.
3. **Photodesorption** — UV photons desorb CO ice back to the gas phase, creating a warm molecular layer between the hot surface and cold midplane.

The network is solved on a cylindrical grid (`nr_cyl_LIME × nz_cyl_LIME`) after VHSE convergence. The stellar UV field strength at 1 au is set by `G0Hab_set` (in Habing units at the stellar surface).

**Reliability:** C¹⁸O emission from the warm molecular layer (intermediate heights, $r \sim 50$–300 au) is well-reproduced. The cold midplane CO abundance is more uncertain; it depends sensitively on the assumed ice binding energy and grain surface area.

**Relevant parameters:** `G0Hab_set`, `nr_cyl_LIME`, `nz_cyl_LIME`, `chemical_save_name`

**Primary reference:** [Ruaud, Gorti & Hollenbach (2022), ApJ 925, 49](https://ui.adsabs.harvard.edu/abs/2022ApJ...925...49R)

---

## UV Propagation and the `chi` / `G0Hab` Parameter

The UV radiation field drives photodissociation. DiskMINT parameterises the stellar UV with `G0Hab_set` — the flux at the stellar surface in Habing units ($G_0 = 1.6 \times 10^{-3}$ erg cm⁻² s⁻¹). The field then scales as $r^{-2}$ outward from the star.

Representative values by stellar mass (from the BT-Settl models used in the HPC example):

| Stellar mass (M☉) | `G0Hab_set` |
|---|---|
| 0.1 | 2.1 × 10⁹ |
| 0.3 | 7.4 × 10⁹ |
| 0.5 | 1.8 × 10¹⁰ |
| 0.7 | 3.5 × 10¹⁰ |
| 1.0 | 6.3 × 10¹⁰ |
| 2.0 | 2.1 × 10¹¹ |

For targets with UV excess (e.g., T Tauri stars with accretion), the UV contribution from accretion shocks should be added to the photospheric value. Setting `G0Hab_set` an order of magnitude higher than the photospheric value is a reasonable first estimate for actively accreting sources.

---

## Key Papers

| Paper | What it covers |
|---|---|
| [Deng et al. (2023), ApJ 954, 165](https://ui.adsabs.harvard.edu/abs/2023ApJ...954..165D) | Original DiskMINT method: VHSE + chemistry for CO isotopologue masses |
| [Deng et al. (2025), ApJ 995, 98](https://ui.adsabs.harvard.edu/abs/2025ApJ...995...98D) | Application to IM Lup; radially varying g2d ratio |
| [Ruaud, Gorti & Hollenbach (2022), ApJ 925, 49](https://ui.adsabs.harvard.edu/abs/2022ApJ...925...49R) | Reduced chemical network used for C¹⁸O |
| [Dubrulle et al. (1995), A&A 309, 209](https://ui.adsabs.harvard.edu/abs/1995A%26A...309..209D) | Turbulent dust settling prescription |
| [Mathis, Rumpl & Nordsieck (1977), ApJ 217, 425](https://ui.adsabs.harvard.edu/abs/1977ApJ...217..425M) | MRN grain size distribution |
