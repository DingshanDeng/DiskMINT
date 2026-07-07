from __future__ import annotations

import sys
from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


def pytest_addoption(parser):
    parser.addoption(
        "--run-diskmint-examples",
        action="store_true",
        default=False,
        help="run slow DiskMINT example model smoke tests that require external binaries",
    )


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "diskmint_example_run: runs a DiskMINT example model and requires external binaries",
    )
