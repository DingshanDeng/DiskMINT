# Physics Background

A brief summary of the physical models underlying DiskMINT. Each section points to the relevant parameters and the primary literature for depth.

For a dedicated overview of the three core thermochemical ingredients behind DiskMINT, see {doc}`thermochemical_processes`.

---

## Vertical Hydrostatic Equilibrium (VHSE)

DiskMINT solves for a self-consistent 2D disk structure by iterating between two steps:

1. **Radiative transfer** — RADMC-3D computes the dust (and gas) temperature by propagating stellar photons through the density structure via Monte Carlo (controlled by `nphot`, `nthreads`).
2. **Vertical hydrostatic equilibrium** — given the temperature, the vertical pressure gradient is balanced against gravity to update the density structure:

   $$\frac{\partial P}{\partial z} = -\rho_\mathrm{gas} \frac{GM_\star z}{(r^2 + z^2)^{3/2}}$$

   where $P = \rho_\mathrm{gas} c_s^2$ with the local isothermal sound speed $c_s = \sqrt{k_B T / \mu m_H}$.

The loop runs until the relative change in gas density between iterations falls below 5 % (`delta_rho/rho_max < 0.05`), or until `n_vhse_loop` iterations are reached. Roughly 10–20 iterations are required to converge from a Gaussian initial condition.

**Relevant parameters:** `n_vhse_loop`, `nphot`, `nthreads`, `nr`, `ntheta`

**Primary reference:** [Deng et al. (2023), ApJ 954, 165](https://ui.adsabs.harvard.edu/abs/2023ApJ...954..165D)

---

## Dust Settling

When `bool_dust_settling = True`, DiskMINT computes dust settling by balancing gravitational settling toward the midplane against upward turbulent diffusion. In the current framework, this is not implemented as a simple post-processing scale-height correction. Instead, as descrbied in [Deng et al. (2025)](https://ui.adsabs.harvard.edu/abs/2025ApJ...995...98D) (Section 3.2.1, Figure 4), DiskMINT solves the steady-state vertical diffusion equation for the dust-to-gas density ratio:

$$\frac{\partial}{\partial z}\left(\ln{\frac{\rho_{\rm dust}}{\rho_{\rm gas}}}\right) = -\frac{\Omega^2 \tau_{\rm s}}{D}z$$

where $\Omega$ is the Keplerian angular frequency, $\tau_{\rm s}$ is the dust stopping time, and $D$ is the turbulent diffusion coefficient. The stopping time depends on the local grain size, gas density, and sound speed, so the settling strength is coupled to the local thermodynamic structure of the disk rather than being set by radius alone.

DiskMINT solves this dust diffusion equation inside the same iterative loop used to converge the vertical disk structure. Starting from an initially well-mixed gas and dust distribution, the code iterates through:

1. radiative transfer to estimate the dust temperature,
2. thermal balance to estimate the gas temperature,
3. vertical hydrostatic equilibrium to update the gas density,
4. vertical dust diffusion to update the density distribution of each grain size bin.

The vertically integrated dust surface density at each radius and grain size is kept fixed during this step; settling only redistributes dust vertically. The result is a self-consistent stratified structure in which larger grains become concentrated toward the midplane while smaller grains remain more vertically extended in the disk surface layers.

In the IM Lup application of Deng et al. (2025), this treatment improves the match to the observed SED, especially near the far-infrared, while having little effect on the inferred total gas mass. Figure 4 of that paper gives a useful flow chart of the settling workflow implemented in DiskMINT.

**Relevant parameters:** `visc_alpha`, `nr_dust_1`, `amin_1`, `amax_1`, `pla_dustsize`

**Reference:** [Deng et al. (2025), ApJ 995, 98](https://ui.adsabs.harvard.edu/abs/2025ApJ...995...98D); see also [Estrada et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016ApJ...818..200E)

---

## Grain Size Distribution

DiskMINT models the grain population as a power law in grain radius $a$:

$$n(a) \propto a^{-p}$$

for $a_\mathrm{min} \le a \le a_\mathrm{max}$, discretised into `nr_dust_1` logarithmically spaced size bins. The MRN distribution (Mathis, Rumpl & Nordsieck 1977) uses $p = 3.5$ (`pla_dustsize`). Each bin has its own opacity table (`dustkappa_{name}_{i}.inp`).

The code computes cross-section-weighted averages (`ndsd.inp`), mass fractions (`fracs.inp`), and number fractions (`fracs_numb.inp`) automatically when `bool_MakeDustKappa = True`.

**Relevant parameters:** `amin_1`, `amax_1`, `nr_dust_1`, `pla_dustsize`, `rhobulk`

---

## Surface Density Profile

DiskMINT supports more than one way of defining the radial surface density structure for gas and dust.

### Analytic surface density

The simplest option is an analytic parameterization. A commonly used choice is the tapered power-law profile:

$$\Sigma(r) \propto \left(\frac{r}{R_\mathrm{tap}}\right)^{\gamma} \exp\!\left[-\left(\frac{r}{R_\mathrm{tap}}\right)^{\xi}\right]$$

where $\gamma =$ `pl_sufdens` and $\xi =$ `pl_tapoff`. This is the self-similar solution from viscous disk evolution (Lynden-Bell & Pringle 1974), which gives $\xi = 2 + \gamma$ (i.e., `pl_tapoff = 2 + pl_sufdens` for a physically motivated viscous disk). The parameter `Rtap` is the characteristic radius (often written as $R_\mathrm{c}$) and sets the scale where the profile begins to taper. This analytic form is useful for forward modeling, parameter studies, and cases where no spatially resolved constraints are available.

### Observation-driven surface density fitting

DiskMINT can also use directly specified radial surface density profiles, including profiles derived from resolved observations. In the structured IM Lup model of [Deng et al. (2025), ApJ 995, 98](https://ui.adsabs.harvard.edu/abs/2025ApJ...995...98D), the dust and gas surface densities are treated separately as $\Sigma_\mathrm{dust}(r)$ and $\Sigma_\mathrm{gas}(r)$, so the dust-to-gas ratio becomes a function of radius rather than a single global constant.

In that workflow, $\Sigma_\mathrm{dust}(r)$ is iteratively adjusted to reproduce the observed dust continuum radial profile, while $\Sigma_\mathrm{gas}(r)$ is iteratively adjusted to reproduce the observed $\mathrm{C^{18}O}$ radial profile. The synthetic and observed radial profiles are measured in elliptical annuli using the observed disk geometry, and the next surface-density estimate is updated from the ratio of observation to model at each radius.

This observation-driven mode is useful when high-quality resolved continuum and line data are available and the gas and dust are radially decoupled. In practice, the analytic and direct-fitting approaches are complementary: the analytic profile is a natural starting point for general modeling, while direct fitting is the more flexible option for well-resolved disks. A worked example of this fitting workflow will be added to the documentation in a future update.

**Relevant parameters:** `pl_sufdens`, `pl_tapoff`, `Rtap`, `mdiskd`, `ratio_g2d_global`

---

## Chemical Network

The chemical network is a reduced version of the full CO chemistry framework developed in [Ruaud & Gorti (2019)](https://ui.adsabs.harvard.edu/abs/2019ApJ...885..146R/abstract) and later adapted for $\mathrm{C^{18}O}$ modeling in [Ruaud, Gorti & Hollenbach (2022)](https://ui.adsabs.harvard.edu/abs/2022ApJ...925...49R/abstract). It focuses on the processes that control $\mathrm{C^{18}O}$ abundance:

- **Isotope-selective photodissociation** — UV photons from the star preferentially dissociate $\mathrm{C^{18}O}$ and $\mathrm{{}^{13}CO}$ over $\mathrm{{}^{12}CO}$ due to self-shielding differences. This is the dominant mechanism depleting $\mathrm{C^{18}O}$ in the upper disk layers, following [Visser, van Dishoeck & Black (2009)](https://www.aanda.org/10.1051/0004-6361/200912129).
- **Grain-surface chemistry** — CO freezes onto grain surfaces at $T < 20$–30 K and converts to $\mathrm{CO_2}$ ice via grain-surface reactions. This depletes gas-phase CO in the cold midplane. Then it is balanced with **Photodesorption** — UV photons desorb CO ice back to the gas phase, creating a warm molecular layer between the hot surface and cold midplane.

The network is solved on a cylindrical grid (`nr_cyl_LIME x nz_cyl_LIME`) after VHSE convergence. 

There are two ways of setting of the UV field. If the UV spectra is included in the input stellar spectra, the UV filed is solved from the input, and in that case, the stellar UV field strength at 1 au is set by `G0Hab_set` (in Habing units at the stellar surface) with `G0Hab_set = None`. 
If the input UV spectra is not available, you may set up the `G0Hab_set` value to the value you prefer.

**Relevant parameters:** `G0Hab_set`, `nr_cyl_LIME`, `nz_cyl_LIME`, `chemical_save_name`

**Primary references:** [Visser, van Dishoeck & Black (2009), A&A 503, 323](https://www.aanda.org/10.1051/0004-6361/200912129); [Ruaud & Gorti (2019), ApJ 885, 146](https://ui.adsabs.harvard.edu/abs/2019ApJ...885..146R); [Ruaud, Gorti & Hollenbach (2022), ApJ 925, 49](https://ui.adsabs.harvard.edu/abs/2022ApJ...925...49R)

---

## Key Papers

| Paper | What it covers |
|---|---|
| [Mathis, Rumpl & Nordsieck (1977), ApJ 217, 425](https://ui.adsabs.harvard.edu/abs/1977ApJ...217..425M) | MRN grain size distribution |
| [Visser, van Dishoeck & Black (2009), A&A 503, 323](https://www.aanda.org/10.1051/0004-6361/200912129) | Original isotope-selective CO photodissociation model |
| [Estrada et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016ApJ...818..200E) | Find the self-consistent disk veritcal structure |
| [Ruaud & Gorti (2019), ApJ 885, 146](https://ui.adsabs.harvard.edu/abs/2019ApJ...885..146R) | Original CO chemistry model underlying the reduced $\mathrm{C^{18}O}$ network |
| [Ruaud, Gorti & Hollenbach (2022), ApJ 925, 49](https://ui.adsabs.harvard.edu/abs/2022ApJ...925...49R) | Reduced chemical network used for $\mathrm{C^{18}O}$ |
| [Deng et al. (2023), ApJ 954, 165](https://ui.adsabs.harvard.edu/abs/2023ApJ...954..165D) | Original DiskMINT method: VHSE + chemistry for CO isotopologue masses |
| [Deng et al. (2025), ApJ 995, 98](https://ui.adsabs.harvard.edu/abs/2025ApJ...995...98D) | Application to IM Lup; radially varying g2d ratio |
