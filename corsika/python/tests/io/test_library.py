"""
Tests for `corsika.io.outputs.energy_loss`

(c) Copyright 2018 CORSIKA Project, corsika8@kit.edu

This software is distributed under the terms of the 3-clause BSD license.
See file LICENSE for a full version of the license.
"""

import os
import os.path as op
import subprocess
import tarfile
import tempfile

import pytest
import yaml

from corsika8.io import Library, outputs

from .. import build_directory

bindir = op.join(build_directory, "tests/output")


def generate_data() -> str:
    """
    Generate a test with `testOutput`.

    Returns
    -------
    str
        The path to the generated data.
    """

    # Expected output directory that is made from the testOutput binary
    output_dir = op.join(os.getcwd(), "out_test/check")

    if not op.exists(output_dir):  # only make if not already run
        binary = op.join(bindir, "testOutput")

        # ensure that the binary exists (not trivial on the CI)
        if not op.exists(binary):
            msg = f"Could not find testOutput binary at {binary}\n"
            msg += f"Binary dir contains {os.listdir(bindir)}"
            raise RuntimeError(msg)

        subprocess.call([binary, "OutputManager"])

        # Check if it still doesn't exist
        if not op.exists(output_dir):
            msg = "After running binary, could not find expected"
            msg += f" output dir {output_dir}\n"
            msg += "The binary did not execute successfully or the"
            msg += " OutputManager tests have changed"
            raise RuntimeError(msg)

    return output_dir


def test_basic_Library() -> None:
    dir_to_test = generate_data()
    lib = Library(dir_to_test)

    assert lib is not None
    assert 0 == len(lib.names)
    assert isinstance(lib.summary, dict)
    assert isinstance(lib.config, dict)

    # Check what happens for an unknown output subdir
    assert lib.get("Does not exits") is None


def test_bad_Library() -> None:
    with pytest.raises(ValueError):
        Library("This does not exist")


# Dummy output class for testing
class DummyOutput(outputs.Output):
    def __init__(self, path):
        self._path = path

    def is_good(self):
        return True

    def astype(self, dtype: str, **kwargs):
        return super().astype(dtype, **kwargs)


# Inject DummyOutput into outputs namespace
outputs.MockOutput = DummyOutput  # type: ignore


@pytest.fixture
def valid_library_dir():
    temp_dir = tempfile.TemporaryDirectory()
    path = temp_dir.name

    config = {"creator": "CORSIKA8"}
    with open(os.path.join(path, "config.yaml"), "w") as f:
        yaml.dump(config, f)

    summary = {"output_dirs": ["output1"]}
    with open(os.path.join(path, "summary.yaml"), "w") as f:
        yaml.dump(summary, f)

    output_path = os.path.join(path, "output1")
    os.mkdir(output_path)
    output_config = {"name": "output1", "type": "MockOutput"}
    with open(os.path.join(output_path, "config.yaml"), "w") as f:
        yaml.dump(output_config, f)

    yield path
    temp_dir.cleanup()


def test_library_initialization(valid_library_dir):
    lib = Library(valid_library_dir)

    assert lib.config is not None
    assert lib.summary is not None
    assert "output1" in lib.names

    output = lib.get("output1")
    assert isinstance(output, DummyOutput)
    assert output.is_good()

    assert lib.get("DOES_NOT_EXIST") is None


def test_library_from_tar(valid_library_dir):
    # Create a tarball from the valid library directory
    with tempfile.NamedTemporaryFile(suffix=".tar") as tar_file:
        tar_path = tar_file.name
        base_name = os.path.basename(valid_library_dir)

        with tarfile.open(tar_path, "w") as tar:
            tar.add(valid_library_dir, arcname=base_name)

        # Create Library from tarball
        lib = Library(tar_path)

        assert lib.config is not None
        assert lib.summary is not None
        assert "output1" in lib.names

        output = lib.get("output1")
        assert isinstance(output, DummyOutput)
        assert output.is_good()

    # do it again but with missing parts
    def add_without(name):
        with tempfile.NamedTemporaryFile(suffix=".tar") as tar_file:
            with tarfile.open(tar_file.name, "w") as tar:
                for element in os.listdir(valid_library_dir):
                    if element == name:
                        continue
                    full_path = os.path.join(valid_library_dir, element)
                    tar.add(full_path, arcname=element)

            # Create Library from tarball
            lib = Library(tar_file.name)
            return lib

    lib = add_without("summary.yaml")
    assert lib.config is not None
    assert lib.summary is None
    with pytest.raises(ValueError):
        add_without("config.yaml")  # invalid
    lib = add_without("output1")
    assert lib.config is not None
    assert lib.summary is not None


def test_invalid_library_raises():
    with tempfile.TemporaryDirectory() as path:
        with pytest.raises(ValueError):
            Library(path)
