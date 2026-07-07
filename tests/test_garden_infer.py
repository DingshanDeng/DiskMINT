from __future__ import annotations

import importlib.util

import pytest

from diskmint.garden import infer


def _has_garden_deps() -> bool:
    return all(
        importlib.util.find_spec(name) is not None
        for name in ("pandas", "astropy", "joblib", "sklearn", "xgboost")
    )


def test_bundled_model_paths_exist():
    assert infer.get_model_path("band6").name == "diskmint_xgb_band6.joblib"
    assert infer.get_model_path("band7").name == "diskmint_xgb_band7.joblib"


def test_invalid_band_name_raises():
    with pytest.raises(ValueError, match="Unknown GARDEN band"):
        infer.get_model_path("band8")


@pytest.mark.skipif(
    importlib.util.find_spec("xgboost") is not None,
    reason="only relevant when xgboost is absent",
)
def test_load_model_reports_missing_xgboost():
    with pytest.raises(infer.GardenDependencyError, match="diskmint\\[garden\\]"):
        infer.load_model("band6")


@pytest.mark.skipif(not _has_garden_deps(), reason="requires diskmint[garden] dependencies")
def test_loads_bundled_models():
    model_b6 = infer.load_model("band6")
    model_b7 = infer.load_model("band7")
    assert model_b6["method"] == "xgb"
    assert model_b7["method"] == "xgb"
    assert model_b6["meta"]["obs_cols"] == ("Lmm_model", "LC18O_model", "R90dust_model")
    assert model_b7["meta"]["obs_cols"] == ("Lmm_model", "LC18O_model", "R90dust_model")


@pytest.mark.skipif(not _has_garden_deps(), reason="requires diskmint[garden] dependencies")
def test_from_observations_smoke_band6():
    result = infer.from_observations(
        flux_mm=120.0,
        flux_c18o=850.0,
        distance=140.0,
        mstar=0.8,
        rdust_90=80.0,
        band="band6",
    )
    for key in ("Mgas_pred_Msun", "Mdust_pred_Msun", "gtd_pred_dimless", "Rc_pred_au"):
        assert key in result
        assert result[key] > 0
    assert "is_outside_grid" in result
    assert "outside_features" in result


@pytest.mark.skipif(not _has_garden_deps(), reason="requires diskmint[garden] dependencies")
def test_from_dataframe_requires_columns():
    import pandas as pd

    with pytest.raises(ValueError, match="Missing required column"):
        infer.from_dataframe(pd.DataFrame({"flux_mm_mjy": [120.0]}))
