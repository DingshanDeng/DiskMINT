from __future__ import annotations

import importlib.util
import inspect
import os
import shutil
import subprocess
import sys
from argparse import Namespace
from pathlib import Path

import numpy as np
import pytest

import diskmint.model as model


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE_DIR = ROOT / "examples" / "example_diskmint_models" / "example_radial_profile_fit"
SCRIPT_PATH = EXAMPLE_DIR / "example_model_radial_profile_fit.py"


def _load_example_module():
    spec = importlib.util.spec_from_file_location("example_model_radial_profile_fit", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_radial_profile_fit_parameter_file_loads():
    para = model.Parameters()
    para.read_parameters_from_csv(
        filename="radial_profile_fit_parameters",
        directory=str(EXAMPLE_DIR),
        extension=".csv",
    )

    assert para.fmodel_filename == "BTSettl_1p0Msolar_1pc_um_ergpcm2hzs.inp"
    assert para.nr == 30
    assert para.ntheta == 100
    assert para.ratio_g2d_global == 100.0


def test_radial_profile_fit_reference_tables_load():
    module = _load_example_module()
    model_inputs = EXAMPLE_DIR / "model_inputs"

    sigmad = module.load_reference_table(model_inputs / "sigma_reference_template.dat")
    ratio_g2d = module.load_reference_table(model_inputs / "ratio_g2d_reference_template.dat")

    assert sigmad.shape == (215, 2)
    assert ratio_g2d.shape == (27, 2)
    assert np.all(np.diff(sigmad[:, 0]) > 0)
    assert np.all(np.diff(ratio_g2d[:, 0]) > 0)
    assert np.all(sigmad[:, 1] > 0)
    assert np.all(ratio_g2d[:, 1] > 0)


def test_radial_profile_fit_modes_configure_expected_parameters():
    module = _load_example_module()

    para_a = model.Parameters()
    para_a.read_parameters_from_csv(
        filename="radial_profile_fit_parameters",
        directory=str(EXAMPLE_DIR),
        extension=".csv",
    )
    module.configure_model(para_a, str(EXAMPLE_DIR), "model_A_constant_g2d")
    assert para_a.chemical_save_name == "example_radial_profile_fit_model_A_constant_g2d"
    assert para_a.ratio_g2d_global == 100.0
    assert para_a.get_parameter_value("ratio_g2d") is None
    assert para_a.sigmad_ref.shape == (215, 2)

    para_b = model.Parameters()
    para_b.read_parameters_from_csv(
        filename="radial_profile_fit_parameters",
        directory=str(EXAMPLE_DIR),
        extension=".csv",
    )
    module.configure_model(para_b, str(EXAMPLE_DIR), "model_B_radial_profile_fit")
    assert para_b.chemical_save_name == "example_radial_profile_fit_model_B_radial_g2d"
    assert para_b.ratio_g2d_global == 100.0
    assert para_b.ratio_g2d.shape == (27, 2)
    assert para_b.sigmad_ref.shape == (215, 2)


def test_radial_profile_fit_uses_only_stable_runtime_flags(tmp_path):
    module = _load_example_module()
    para = model.Parameters()
    para.read_parameters_from_csv(
        filename="radial_profile_fit_parameters",
        directory=str(EXAMPLE_DIR),
        extension=".csv",
    )
    module.configure_model(para, str(EXAMPLE_DIR), "model_B_radial_profile_fit")

    args = Namespace(n_vhse_loop=3, dust_settling=True, skip_chemistry=True)
    module.copy_input_files(str(EXAMPLE_DIR / "input"), str(tmp_path), para)
    mint = module.build_mint(para, str(tmp_path), args)

    assert mint.bool_MakeDustKappa is True
    assert mint.bool_SED is True
    assert mint.bool_VHSE is True
    assert mint.n_vhse_loop == 3
    assert mint.bool_dust_settling is True
    assert mint.bool_chemistry is False
    assert mint.bool_savemodel is True
    assert mint.chem_code_dir == "reducedRGH22"

    build_mint_source = inspect.getsource(module.build_mint)
    for attr in (
        "bool_temp_decouple",
        "bool_dust_fragmentation",
        "bool_dust_radial_drifting",
        "bool_dust_inner_rim",
        "bool_same_rc_as_radmc3d",
    ):
        assert attr not in build_mint_source


def test_radial_profile_fit_copies_required_inputs(tmp_path):
    module = _load_example_module()
    para = model.Parameters()
    para.read_parameters_from_csv(
        filename="radial_profile_fit_parameters",
        directory=str(EXAMPLE_DIR),
        extension=".csv",
    )

    module.copy_input_files(str(EXAMPLE_DIR / "input"), str(tmp_path), para)

    assert (tmp_path / para.fmodel_filename).is_file()
    assert len(list(tmp_path.glob("dustkappa_optool_20250122_*.inp"))) == 20
    for filename in ("aave.inp", "fracs.inp", "fracs_numb.inp", "ndsd.inp"):
        assert (tmp_path / filename).is_file()


@pytest.mark.diskmint_example_run
def test_radial_profile_fit_short_model_run_when_enabled(pytestconfig, tmp_path):
    if not pytestconfig.getoption("--run-diskmint-examples"):
        pytest.skip("requires --run-diskmint-examples and a configured RADMC-3D environment")

    run_dir = tmp_path / "example_radial_profile_fit"
    shutil.copytree(EXAMPLE_DIR, run_dir)

    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src") + os.pathsep + env.get("PYTHONPATH", "")
    command = [
        sys.executable,
        "-u",
        "example_model_radial_profile_fit.py",
        "--mode",
        "model_A_constant_g2d",
        "--n-vhse-loop",
        "1",
        "--skip-chemistry",
    ]
    result = subprocess.run(
        command,
        cwd=run_dir,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=3600,
        check=False,
    )

    assert result.returncode == 0, result.stdout
    assert (run_dir / "data" / "example_radial_profile_fit_model_A_constant_g2d_parameters_setup.csv").is_file()
