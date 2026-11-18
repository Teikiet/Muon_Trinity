"""
Tests for InteractionHistogram IO reader.

(c) Copyright 2025 CORSIKA Project, corsika8@kit.edu
"""

import os
import os.path as op
import tempfile

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
import yaml

from corsika8.io import Library
from corsika8.io.outputs import InteractionHistogram


def _write_minimal_lab(path: str, shower: int, pdg_dict: dict) -> None:
    """Build a minimal lab parquet matching writer schema (shower + edges + pdg_*)."""
    n = len(pdg_dict["energy_low_eV"])
    rows = [
        [shower] * n,
        pdg_dict["energy_low_eV"],
        pdg_dict["energy_high_eV"],
    ]
    names = ["shower", "energy_low_eV", "energy_high_eV"]
    arrays = [pa.array(rows[0]), pa.array(rows[1]), pa.array(rows[2])]
    for k, v in pdg_dict.items():
        if k.startswith("pdg_"):
            names.append(k)
            arrays.append(pa.array(v))
    table = pa.table(dict(zip(names, arrays)))
    pq.write_table(table, path)


def _write_minimal_cms(path: str, shower: int, pdg_dict: dict) -> None:
    n = len(pdg_dict["sqrt_s_low_eV"])
    names = ["shower", "sqrt_s_low_eV", "sqrt_s_high_eV"]
    arrays = [
        pa.array([shower] * n),
        pa.array(pdg_dict["sqrt_s_low_eV"]),
        pa.array(pdg_dict["sqrt_s_high_eV"]),
    ]
    for k, v in pdg_dict.items():
        if k.startswith("pdg_"):
            names.append(k)
            arrays.append(pa.array(v))
    table = pa.table(dict(zip(names, arrays)))
    pq.write_table(table, path)


@pytest.fixture
def interaction_hist_library() -> str:
    """Minimal valid library with one interaction_hist output and two showers."""
    tmp = tempfile.TemporaryDirectory()
    root = tmp.name

    with open(op.join(root, "config.yaml"), "w") as f:
        yaml.dump({"creator": "CORSIKA8"}, f)

    with open(op.join(root, "summary.yaml"), "w") as f:
        yaml.dump({"output_dirs": ["interaction_hist"]}, f)

    sub = op.join(root, "interaction_hist")
    os.makedirs(sub)

    with open(op.join(sub, "config.yaml"), "w") as f:
        yaml.dump(
            {
                "name": "interaction_hist",
                "type": "InteractionHistogram",
                "format": "parquet",
            },
            f,
        )

    summary_shower = {
        "rows": 2,
        "integrated_interactions_per_pdg": {"2212": 2.0},
        "projectile_pdgs_lab_energy": [2212],
        "projectile_pdgs_cms_energy": [2212],
    }
    with open(op.join(sub, "summary.yaml"), "w") as f:
        yaml.dump(
            {
                "shower_0": summary_shower,
                "shower_1": summary_shower,
                "total_rows": 8,
            },
            f,
        )

    # Shower 0: only pdg_2212
    lab0 = {
        "energy_low_eV": [1e3, 2e3],
        "energy_high_eV": [2e3, 3e3],
        "pdg_2212": [1.0, 0.0],
    }
    cms0 = {
        "sqrt_s_low_eV": [1e4, 2e4],
        "sqrt_s_high_eV": [2e4, 3e4],
        "pdg_2212": [1.0, 0.0],
    }
    _write_minimal_lab(op.join(sub, "lab_energy_0.parquet"), 0, lab0)
    _write_minimal_cms(op.join(sub, "cms_energy_0.parquet"), 0, cms0)

    # Shower 1: adds pdg_2112 in second bin only
    lab1 = {
        "energy_low_eV": [1e3, 2e3],
        "energy_high_eV": [2e3, 3e3],
        "pdg_2212": [0.0, 1.0],
        "pdg_2112": [0.0, 2.0],
    }
    cms1 = {
        "sqrt_s_low_eV": [1e4, 2e4],
        "sqrt_s_high_eV": [2e4, 3e4],
        "pdg_2212": [0.0, 1.0],
        "pdg_2112": [0.0, 2.0],
    }
    _write_minimal_lab(op.join(sub, "lab_energy_1.parquet"), 1, lab1)
    _write_minimal_cms(op.join(sub, "cms_energy_1.parquet"), 1, cms1)

    yield root
    tmp.cleanup()


def test_interaction_histogram_combined_and_average(
    interaction_hist_library: str,
) -> None:
    root = interaction_hist_library
    lib = Library(root)
    out = lib.get("interaction_hist")
    assert isinstance(out, InteractionHistogram)
    assert out.is_good()

    lab = out.combined_lab()
    assert len(lab) == 4
    assert "pdg_2212" in lab.columns and "pdg_2112" in lab.columns
    # Shower 0: pdg_2112 missing in file -> filled with 0
    assert lab.loc[lab["shower"] == 0, "pdg_2112"].tolist() == [0.0, 0.0]

    avg = out.average_lab_histogram()
    assert set(avg.columns) == {
        "energy_low_eV",
        "energy_high_eV",
        "pdg_2112",
        "pdg_2212",
    }
    assert len(avg) == 2
    row0 = avg[(avg["energy_low_eV"] == 1e3) & (avg["energy_high_eV"] == 2e3)].iloc[0]
    assert row0["pdg_2212"] == 0.5
    assert row0["pdg_2112"] == 0.0

    data = out.astype("pandas")
    assert "lab" in data and "cms" in data
    pd.testing.assert_frame_equal(data["lab"], lab)


def test_interaction_histogram_integrated(interaction_hist_library: str) -> None:
    root = interaction_hist_library
    lib = Library(root)
    out = lib.get("interaction_hist")
    assert isinstance(out, InteractionHistogram)
    assert out.integrated_interactions(0) == {"2212": 2.0}
