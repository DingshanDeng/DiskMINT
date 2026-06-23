"""Predict-only DiskMINT-GARDEN inference API.

This module exposes the beta inference interface for the bundled DiskMINT-GARDEN
XGBoost surrogate models. Training and diagnostic plotting helpers remain
development-only in the research repository.
"""

from __future__ import annotations

from importlib import resources
from typing import Any

import numpy as np

_BAND_FILES = {
    "band6": "diskmint_xgb_band6.joblib",
    "b6": "diskmint_xgb_band6.joblib",
    "c18o21": "diskmint_xgb_band6.joblib",
    "2-1": "diskmint_xgb_band6.joblib",
    "band7": "diskmint_xgb_band7.joblib",
    "b7": "diskmint_xgb_band7.joblib",
    "c18o32": "diskmint_xgb_band7.joblib",
    "3-2": "diskmint_xgb_band7.joblib",
}

_C18O_REST_FREQUENCIES_GHZ = {
    "band6": 219.5603541,
    "b6": 219.5603541,
    "c18o21": 219.5603541,
    "2-1": 219.5603541,
    "band7": 329.3305525,
    "b7": 329.3305525,
    "c18o32": 329.3305525,
    "3-2": 329.3305525,
}


class GardenDependencyError(ImportError):
    """Raised when optional DiskMINT-GARDEN dependencies are unavailable."""


def _require_optional_dependencies() -> dict[str, Any]:
    missing = []
    modules: dict[str, Any] = {}
    for name in ("pandas", "astropy.units", "astropy.constants", "joblib"):
        try:
            if name == "pandas":
                import pandas as pd

                modules["pd"] = pd
            elif name == "astropy.units":
                import astropy.units as u

                modules["u"] = u
            elif name == "astropy.constants":
                import astropy.constants as const

                modules["const"] = const
            elif name == "joblib":
                import joblib

                modules["joblib"] = joblib
        except ImportError:
            missing.append(name.split(".", maxsplit=1)[0])
    if missing:
        deps = ", ".join(sorted(set(missing)))
        raise GardenDependencyError(
            "DiskMINT-GARDEN inference requires optional dependencies "
            f"({deps}). Install them with `pip install 'diskmint[garden]'`."
        )
    return modules


def _normalise_band(band: str) -> str:
    key = str(band).lower().replace("_", "").replace(" ", "")
    if key not in _BAND_FILES:
        valid = "band6, band7"
        raise ValueError(f"Unknown GARDEN band '{band}'. Use one of: {valid}.")
    return key


def get_model_path(band: str = "band6"):
    """Return the packaged path-like object for a bundled GARDEN model."""

    key = _normalise_band(band)
    return resources.files("diskmint.garden.data").joinpath(_BAND_FILES[key])


def load_model(band: str = "band6") -> dict[str, Any]:
    """Load one bundled DiskMINT-GARDEN surrogate model.

    Parameters
    ----------
    band
        `band6` for the C18O J=2-1 model, or `band7` for the C18O J=3-2 model.
    """

    modules = _require_optional_dependencies()
    try:
        import xgboost  # noqa: F401
    except ImportError as exc:
        raise GardenDependencyError(
            "DiskMINT-GARDEN bundled models require xgboost. "
            "Install it with `pip install 'diskmint[garden]'`."
        ) from exc
    path = get_model_path(band)
    with resources.as_file(path) as model_path:
        try:
            return modules["joblib"].load(model_path)
        except ImportError as exc:
            raise GardenDependencyError(
                "DiskMINT-GARDEN could not load the bundled model because an "
                "optional dependency is missing. Install with "
                "`pip install 'diskmint[garden]'`."
            ) from exc


def _safe_log10(x, eps: float = 1e-40):
    x = np.asarray(x, dtype=float)
    return np.log10(np.clip(x, eps, None))


def _flux2luminosity(flux_intg, dpc, restfreq_ghz):
    """Convert integrated line flux in mJy km/s to line luminosity in Lsun."""

    modules = _require_optional_dependencies()
    u = modules["u"]
    const = modules["const"]
    conv_factor = restfreq_ghz * u.GHz / const.c
    return (
        4 * np.pi * flux_intg * conv_factor * u.mJy * u.km / u.s * (dpc * u.pc) ** 2
    ).to(u.Lsun).value


def _luminosity_from_flux_density(
    flux_mjy,
    d_pc,
    *,
    wavelength=None,
    frequency=None,
):
    """Convert continuum flux density in mJy to monochromatic nu Lnu in Lsun."""

    modules = _require_optional_dependencies()
    u = modules["u"]
    const = modules["const"]
    if (wavelength is None) == (frequency is None):
        raise ValueError("Provide exactly one of `wavelength` or `frequency`.")

    flux_mjy = np.asanyarray(flux_mjy)
    s_nu = (flux_mjy * 1e-3) * u.Jy
    distance = (d_pc * u.pc).to(u.m)
    if wavelength is not None:
        wave = wavelength if isinstance(wavelength, u.Quantity) else wavelength * u.mm
        nu = (const.c / wave.to(u.m)).to(u.Hz)
    else:
        freq = frequency if isinstance(frequency, u.Quantity) else frequency * u.GHz
        nu = freq.to(u.Hz)
    lnu_si = (4 * np.pi * distance**2 * s_nu.to(u.W / u.m**2 / u.Hz)).to(u.W / u.Hz)
    return (nu * lnu_si).to(u.Lsun).value


def _domain_check(meta: dict[str, Any], x_new, feature_names: list[str]) -> list[dict[str, Any]]:
    modules = _require_optional_dependencies()
    try:
        from sklearn.neighbors import NearestNeighbors
    except ImportError as exc:
        raise GardenDependencyError(
            "DiskMINT-GARDEN inference requires scikit-learn. "
            "Install it with `pip install 'diskmint[garden]'`."
        ) from exc

    x_new = np.asarray(x_new, dtype=float)
    below = x_new < meta["X_min"]
    above = x_new > meta["X_max"]
    nn_dist = np.full(x_new.shape[0], np.nan)
    if "X_train_domain" in meta:
        nn = NearestNeighbors(n_neighbors=1).fit(meta["X_train_domain"])
        nn_dist, _ = nn.kneighbors(x_new)
        nn_dist = nn_dist[:, 0]

    rows = []
    for i in range(x_new.shape[0]):
        outside = below[i] | above[i]
        bad_idx = np.where(outside)[0]
        rows.append(
            {
                "is_outside_grid": bool(outside.any()),
                "outside_features": ", ".join(feature_names[j] for j in bad_idx),
                "nn_dist": float(nn_dist[i]) if np.isfinite(nn_dist[i]) else np.nan,
            }
        )
    return rows


def _predict_from_luminosities(
    model: dict[str, Any],
    df,
    *,
    lmm_col: str,
    lc18o_col: str,
    rdust_col: str,
    mstar_col: str,
):
    meta = model["meta"]
    missing = [col for col in (lmm_col, lc18o_col, rdust_col, mstar_col) if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required column(s): {missing}")

    x_list = []
    feature_names = []
    for col in meta["obs_cols"]:
        if col == "Lmm_model":
            arr = df[lmm_col].to_numpy(dtype=float)
        elif col == "LC18O_model":
            arr = df[lc18o_col].to_numpy(dtype=float)
        elif col == "R90dust_model":
            arr = df[rdust_col].to_numpy(dtype=float)
        elif col == "R90CO_model":
            raise ValueError("This beta API does not support GARDEN models requiring R90CO.")
        else:
            raise ValueError(f"Unsupported GARDEN feature in model metadata: {col}")
        x_list.append(_safe_log10(arr) if meta["log_obs"] else arr)
        feature_names.append(col)

    if meta["include_mstar"]:
        x_list.append(df[mstar_col].to_numpy(dtype=float))
        feature_names.append("mstar")

    x_new = np.column_stack(x_list)
    flags = _domain_check(meta, x_new, feature_names)
    if meta.get("clip_X_to_grid", True):
        x_new = np.clip(x_new, meta["X_min"], meta["X_max"])

    y_pred = model["reg"].predict(x_new)
    if meta.get("clip_Y_to_grid", True):
        lo = meta.get("Y_q01", meta["Y_min"])
        hi = meta.get("Y_q99", meta["Y_max"])
        y_pred = np.clip(y_pred, lo, hi)

    log_flags = meta.get(
        "param_log_flags",
        {param: meta.get("log_params", True) for param in meta["param_cols"]},
    )
    out = df.copy()
    for j, pname in enumerate(meta["param_cols"]):
        values = y_pred[:, j]
        out[f"{pname}_pred"] = 10.0**values if log_flags.get(pname, False) else values

    out["is_outside_grid"] = [row["is_outside_grid"] for row in flags]
    out["outside_features"] = [row["outside_features"] for row in flags]
    out["nn_dist"] = [row["nn_dist"] for row in flags]
    if "mdust_pred" in out.columns:
        out["Mdust_pred_Msun"] = out["mdust_pred"]
    if "mdisk_pred" in out.columns:
        out["Mgas_pred_Msun"] = out["mdisk_pred"]
    if {"mdisk_pred", "mdust_pred"}.issubset(out.columns):
        out["gtd_pred_dimless"] = out["mdisk_pred"] / out["mdust_pred"]
    if "rc_pred" in out.columns:
        out["Rc_pred_au"] = out["rc_pred"]
    return out


def from_dataframe(
    observations,
    *,
    band: str = "band6",
    model: dict[str, Any] | None = None,
    continuum_flux_col: str = "flux_mm_mjy",
    c18o_flux_col: str = "flux_c18o_mjy_kms",
    distance_col: str = "distance_pc",
    mstar_col: str = "mstar_msun",
    rdust_col: str = "Rdust_90_au",
    wavelength=None,
    frequency=None,
    rest_frequency_ghz: float | None = None,
):
    """Run DiskMINT-GARDEN inference for an observation table.

    Parameters
    ----------
    observations
        A pandas DataFrame containing continuum flux, C18O integrated flux,
        distance, stellar mass, and 90 percent dust radius columns.
    band
        `band6` for C18O J=2-1, or `band7` for C18O J=3-2.
    wavelength, frequency
        Continuum observing wavelength or frequency. If neither is given, the
        default is 1.3 mm.
    rest_frequency_ghz
        Optional C18O line rest frequency override in GHz.
    """

    modules = _require_optional_dependencies()
    pd = modules["pd"]
    u = modules["u"]
    key = _normalise_band(band)
    df = observations.copy()
    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame(df)
    required = [continuum_flux_col, c18o_flux_col, distance_col, mstar_col, rdust_col]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required column(s): {missing}")

    if wavelength is None and frequency is None:
        wavelength = 1.3 * u.mm
    elif wavelength is not None and frequency is not None:
        raise ValueError("Provide only one of `wavelength` or `frequency`.")
    rest = rest_frequency_ghz or _C18O_REST_FREQUENCIES_GHZ[key]

    df["Lmm_Lsun"] = _luminosity_from_flux_density(
        df[continuum_flux_col].to_numpy(dtype=float),
        df[distance_col].to_numpy(dtype=float),
        wavelength=wavelength,
        frequency=frequency,
    )
    df["LC18O_Lsun"] = _flux2luminosity(
        df[c18o_flux_col].to_numpy(dtype=float),
        df[distance_col].to_numpy(dtype=float),
        rest,
    )

    loaded_model = model if model is not None else load_model(key)
    return _predict_from_luminosities(
        loaded_model,
        df,
        lmm_col="Lmm_Lsun",
        lc18o_col="LC18O_Lsun",
        rdust_col=rdust_col,
        mstar_col=mstar_col,
    )


def from_observations(
    *,
    flux_mm,
    flux_c18o,
    distance,
    mstar,
    rdust_90,
    band: str = "band6",
    wavelength=None,
    frequency=None,
    rest_frequency_ghz: float | None = None,
) -> dict[str, Any]:
    """Run DiskMINT-GARDEN inference for one target.

    Flux density is in mJy, C18O integrated flux is in mJy km/s, distance is in
    pc, stellar mass is in Msun, and `rdust_90` is the 90 percent dust radius in
    au.
    """

    modules = _require_optional_dependencies()
    pd = modules["pd"]
    df = pd.DataFrame(
        {
            "flux_mm_mjy": [flux_mm],
            "flux_c18o_mjy_kms": [flux_c18o],
            "distance_pc": [distance],
            "mstar_msun": [mstar],
            "Rdust_90_au": [rdust_90],
        }
    )
    result = from_dataframe(
        df,
        band=band,
        wavelength=wavelength,
        frequency=frequency,
        rest_frequency_ghz=rest_frequency_ghz,
    )
    row = result.iloc[0].to_dict()
    for key, value in list(row.items()):
        if isinstance(value, np.generic):
            row[key] = value.item()
    return row
