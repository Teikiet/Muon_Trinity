"""
Tests for `corsika`

(c) Copyright 2018 CORSIKA Project, corsika8@kit.edu

This software is distributed under the terms of the 3-clause BSD license.
See file LICENSE for a full version of the license.
"""

import os
import os.path as op

__all__ = ["build_directory"]


def find_build_directory() -> str:
    """
    Return the absolute path to the CORSIKA8 build directory.

    Returns
    -------
    str
        The absolute path to the build directory.

    Raises
    ------
    RuntimeError
        If the build directory cannot be found or it is empty.
    """

    # if we are running on Gitlab
    is_on_CI = os.getenv("CI_PROJECT_DIR") is not None
    if is_on_CI:
        build_dir = op.abspath(op.join(os.getenv("CI_PROJECT_DIR"), "build"))

    else:  # otherwise, we are running locally
        here = op.dirname(os.path.realpath(__file__))
        for name in ["build", "corsika-build"]:
            build_dir = op.realpath(
                op.join(
                    here,
                    op.pardir,
                    op.pardir,
                    op.pardir,
                    "corsika-build",
                )
            )
            if op.isdir(build_dir):
                break

    # check that the build directory contains 'CMakeCache.txt'
    if not op.exists(op.join(build_dir, "CMakeCache.txt")):
        msg = "Python tests cannot find C8 build dir.\n"
        msg += f"Is on CI: {is_on_CI}\n"
        msg += f"Path for CI: {os.getenv('CI_PROJECT_DIR')}\n"
        msg += f"Checked in {build_dir}\n"
        msg += f"Contents: {os.listdir(build_dir)}"
        raise RuntimeError(msg)

    return build_dir


# find the build_directory once
build_directory = find_build_directory()
