# Thermochemical Processes in DiskMINT

This page summarizes the three key thermochemical ingredients that form the physical backbone of DiskMINT:

1. vertical hydrostatic equilibrium (VHSE)
2. isotopologue-selective photodissociation
3. grain-surface chemistry, especially CO freeze-out and conversion to CO$_2$ ice

These processes are central to how DiskMINT interprets C$^{18}$O emission and infers the gas structure of a disk. Together, they are grounded in three key DiskMINT references:

- [Deng et al. (2023), ApJ 954, 165](https://ui.adsabs.harvard.edu/abs/2023ApJ...954..165D) for the original DiskMINT framework
- [Deng et al. (2025), ApJ 995, 98](https://ui.adsabs.harvard.edu/abs/2025ApJ...995...98D) for the IM Lup application and the extended structured-model workflow
- [Ruaud, Gorti & Hollenbach (2022), ApJ 925, 49](https://ui.adsabs.harvard.edu/abs/2022ApJ...925...49R) for the reduced chemical network used to model C$^{18}$O

---

## Why thermochemical modeling matters

Weak CO isotopologue emission does not always mean that a disk is gas-poor. The emission can also be reduced because CO is selectively photodissociated, frozen onto grains, chemically converted into other ice species, or shifted into a different vertical emitting layer by the disk structure itself.

DiskMINT is designed to model these effects self-consistently instead of relying on a fixed CO abundance or a simple global depletion factor. In practice, the C$^{18}$O emission depends on where the gas is dense enough, warm or cold enough, and chemically protected enough for the molecule to survive and emit.

## 1. Vertical Hydrostatic Equilibrium (VHSE)

DiskMINT solves for a self-consistent vertical gas structure by iterating between density and temperature. Rather than imposing a fixed Gaussian vertical profile, the code updates the gas density so that the vertical pressure gradient balances the stellar gravitational force.

This matters because the line-emitting layer depends on the density and temperature structure of the disk. If the vertical structure is inconsistent with the thermal structure, the predicted CO isotopologue emission can be biased.

In the current DiskMINT framework:

- radiative transfer determines the dust temperature
- the gas temperature is computed from thermal coupling with dust grains
- the gas density is then updated by solving vertical hydrostatic equilibrium

This iterative treatment is one of the core physical ingredients of the model and is essential for locating the emitting layers of the rarer CO isotopologues. The original DiskMINT implementation of this self-consistent framework is presented in [Deng et al. (2023), ApJ 954, 165](https://ui.adsabs.harvard.edu/abs/2023ApJ...954..165D), and its application to the structured IM Lup modeling is developed further in [Deng et al. (2025), ApJ 995, 98](https://ui.adsabs.harvard.edu/abs/2025ApJ...995...98D).

## 2. Selective Chemistry: Isotopologue-selective Photodissociation

DiskMINT does not assume that C$^{18}$O follows a fixed isotopic scaling from $^{12}$CO. Instead, it includes isotopologue-selective photodissociation in the reduced chemistry network.

The reason is that the more abundant isotopologues self-shield more effectively against UV radiation. Rarer isotopologues such as C$^{18}$O are less protected and can be dissociated deeper into the disk. This changes both the radial and vertical abundance structure of the gas that actually emits in C$^{18}$O.

This effect is especially important when using C$^{18}$O as a gas-mass tracer. Without selective photodissociation, the modeled C$^{18}$O abundance can be overestimated, which can in turn bias the inferred gas surface density or total gas mass. In DiskMINT, this chemistry is inherited from the reduced network of [Ruaud, Gorti & Hollenbach (2022), ApJ 925, 49](https://ui.adsabs.harvard.edu/abs/2022ApJ...925...49R), which preserves the processes most important for C$^{18}$O.

## 3. Grain-surface Chemistry

DiskMINT also includes more than simple CO freeze-out. In cold parts of the disk, CO can freeze onto dust grains, but the chemistry does not stop there. On grain surfaces, CO can be converted into CO$_2$ ice and related ice reservoirs.

This matters because freeze-out with grain-surface conversion removes CO from the gas phase more effectively than a simple reversible freeze-out picture. It changes the location of the effective CO snow surface and reduces the amount of CO available to emit in the cold disk.

In DiskMINT, this process is one of the main ingredients needed to interpret weak CO isotopologue emission physically rather than replacing it with an arbitrary global depletion factor. This treatment follows the reduced chemistry framework of [Ruaud, Gorti & Hollenbach (2022), ApJ 925, 49](https://ui.adsabs.harvard.edu/abs/2022ApJ...925...49R) and is one of the ingredients emphasized in both [Deng et al. (2023), ApJ 954, 165](https://ui.adsabs.harvard.edu/abs/2023ApJ...954..165D) and [Deng et al. (2025), ApJ 995, 98](https://ui.adsabs.harvard.edu/abs/2025ApJ...995...98D).

## How these three processes work together

These ingredients are tightly connected:

- VHSE determines the density and temperature structure of the disk
- selective photodissociation determines where C$^{18}$O can survive under UV irradiation
- grain-surface chemistry determines where CO is removed from the gas and locked into ices

Together, they set the location and strength of the C$^{18}$O emitting layer. This is why DiskMINT treats them as a coupled thermochemical problem rather than as separate post-processing corrections.

## Relationship to the DiskMINT workflow

In practical terms:

- the Python modeling workflow solves the disk structure, including VHSE
- the reduced chemistry network is then run on the converged structure
- the chemistry is designed to retain the processes most important for C$^{18}$O modeling while remaining computationally practical for target-by-target studies

For a broader overview of the physical setup, see {doc}`physics_background`. For the original DiskMINT framework, see [Deng et al. (2023), ApJ 954, 165](https://ui.adsabs.harvard.edu/abs/2023ApJ...954..165D); for the reduced chemistry network, see [Ruaud, Gorti & Hollenbach (2022), ApJ 925, 49](https://ui.adsabs.harvard.edu/abs/2022ApJ...925...49R); and for the structured IM Lup application, see [Deng et al. (2025), ApJ 995, 98](https://ui.adsabs.harvard.edu/abs/2025ApJ...995...98D).
