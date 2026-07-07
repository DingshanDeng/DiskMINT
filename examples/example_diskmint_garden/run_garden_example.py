"""Smoke example for the DiskMINT-GARDEN inference API."""

from __future__ import annotations

import pandas as pd

import diskmint.garden.infer as garden


def main() -> None:
    one_target = garden.from_observations(
        flux_mm=120.0,
        flux_c18o=850.0,
        distance=140.0,
        mstar=0.8,
        rdust_90=80.0,
        band="band6",
    )
    print("Single-target prediction:")
    for key in ("Mgas_pred_Msun", "Mdust_pred_Msun", "gtd_pred_dimless", "Rc_pred_au"):
        print(f"  {key}: {one_target[key]:.4g}")
    print(f"  outside grid: {one_target['is_outside_grid']}")

    table = pd.DataFrame(
        {
            "name": ["demo-band6", "demo-band6-compact"],
            "flux_mm_mjy": [120.0, 45.0],
            "flux_c18o_mjy_kms": [850.0, 180.0],
            "distance_pc": [140.0, 160.0],
            "mstar_msun": [0.8, 0.4],
            "Rdust_90_au": [80.0, 35.0],
        }
    )
    predictions = garden.from_dataframe(table, band="band6")
    print("\nTable prediction:")
    print(
        predictions[
            [
                "name",
                "Mgas_pred_Msun",
                "Mdust_pred_Msun",
                "gtd_pred_dimless",
                "Rc_pred_au",
                "is_outside_grid",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
