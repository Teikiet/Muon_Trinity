import numpy as np

import corsika8


def GetPrimaryDict() -> dict:
    prop_dict = {
        "x": 1,
        "y": 2,
        "z": 3,
        "nx": 1 / np.sqrt(3),
        "ny": 1 / np.sqrt(3),
        "nz": -1 / np.sqrt(3),
        "pdg": 2212,
        "total_energy": 1234.5678,
    }
    return prop_dict


def test_particle_init() -> None:
    prop_dict = GetPrimaryDict()
    prim = corsika8.io.outputs.Particle(prop_dict)

    assert prim.x
    assert prim.y
    assert prim.z

    assert prim.x == 1
    assert prim.y == 2
    assert prim.z == 3

    assert len(prim.position) == 3

    assert prim.nx
    assert prim.ny
    assert prim.nz

    assert len(prim.direction) == 3

    assert prim.pdg


def test_particle_init_extended() -> None:
    prop_dict = GetPrimaryDict()
    prop_dict["some_new_field"] = 77
    prim = corsika8.io.outputs.Particle(prop_dict)

    assert prim.some_new_field == 77
