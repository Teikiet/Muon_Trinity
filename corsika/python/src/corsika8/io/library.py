"""
This file allows for reading/working with C8 libraries.

(c) Copyright 2020 CORSIKA Project, corsika8@kit.edu

This software is distributed under the terms of the 3-clause BSD license.
See file LICENSE for a full version of the license.
"""

import os
import os.path as op
import tarfile
import tempfile
from typing import Any, Dict, List, Optional, Union

import yaml

from . import outputs
from .logger import c8_logger


class Library:
    """
    Represents a library ("run") of showers produced by C8.
    """

    def __init__(self, path: str):
        """

        Parameters
        ----------
        path: str
            The path to the directory containing the library.


        Raises
        ------
        ValueError
            If `path` does not contain a valid CORSIKA8 library.
        """

        # If the path is a tarball, make a temp dir to open the file into
        if path.endswith(".tar"):
            path, self._temp_dir = self._extract_tar(path)

        # Check that this is a valid library
        if not self.__valid_library(path):
            raise ValueError(f"'{path}' does not contain a valid CORSIKA8 library.")

        # store the top-level path
        self.path = path

        output_dirs = None

        # load the config and summary files
        self.config = self.load_config(self.path)
        self.summary = self.load_summary(self.path)

        if self.summary is None:
            msg = f"Missing summary file in '{path}'."
            msg += " The simulation may not have finished. Will not load library"
            c8_logger.warning(msg)
            return

        if "output_dirs" in self.summary.keys():
            output_dirs = self.summary["output_dirs"]
            msg = f"Reading in sub-directories: {output_dirs}"
            c8_logger.debug(msg)
        else:
            msg = "Sub-directories not specified in summary.yaml file."
            msg += " Will find then dynamically"
            c8_logger.debug(msg)

        # build the list of outputs
        self.__outputs = self.__build_outputs(self.path, output_dirs)

    def _extract_tar(self, path: str) -> tuple[str, tempfile.TemporaryDirectory]:
        temp_dir = tempfile.TemporaryDirectory(prefix="c8_temp_dir_")
        msg = f"Unzipping corsika library {path} to {temp_dir.name}"
        c8_logger.info(msg)

        with tarfile.open(path, "r") as tar:
            # filter not available <py3.12, but add back one day
            # tar.extractall(temp_dir.name, filter="fully_trusted")
            tar.extractall(temp_dir.name)

        # Find the first directory inside the temp_dir
        extracted_items = os.listdir(temp_dir.name)
        if (
            "summary.yaml" in extracted_items or "config.yaml" in extracted_items
        ):  # does not include a directory as top level
            extracted_path = temp_dir.name
        else:
            subdirs = [
                item
                for item in extracted_items
                if op.isdir(op.join(temp_dir.name, item))
            ]
            if not subdirs:
                raise ValueError(f"No directory found inside tarball '{path}'")
            extracted_path = op.join(temp_dir.name, subdirs[0])

        return extracted_path, temp_dir

    @property
    def names(self) -> List[str]:
        """
        Return the list of registered outputs.
        """
        return list(self.__outputs.keys())

    def get(self, name: str) -> Optional[outputs.Output]:
        """
        Return the output with a given name.
        """
        if name in self.__outputs:
            return self.__outputs[name]
        else:
            msg = f"Output with name '{name}' not available in this library. Skipping."
            c8_logger.warning(msg)
            return None

    @staticmethod
    def __load_yaml(path: str, yaml_name: str) -> Optional[Dict[str, Any]]:
        """
        Load the yaml from a given library path.


        Parameters
        ----------
        path: str
            The path to the directory containing the library.

        yaml_name: str
            The name of the yaml file in `path`


        Returns
        -------
        dict or None:
            A dict of the yaml file, if valid directory, otherwise returns `None`

        """
        config_file = op.join(path, yaml_name)
        if op.exists(config_file):
            with open(config_file, "r") as f:
                return yaml.load(f, Loader=yaml.Loader)

        return None

    @staticmethod
    def load_config(path: str) -> Optional[Dict[str, Any]]:
        """
        Load the top-level config from a given library path.


        Parameters
        ----------
        path: str
            The path to the directory containing the library.

        Returns
        -------
        dict or None:
            The config if valid directory, otherwise returns None

        """

        return Library.__load_yaml(path, "config.yaml")

    @staticmethod
    def load_summary(path: str) -> Optional[Dict[str, Any]]:
        """
        Load the top-level summary from a given library path.


        Parameters
        ----------
        path: str
            The path to the directory containing the library.

        Returns
        -------
        dict or None:
            The config if valid directory, otherwise returns None

        """

        return Library.__load_yaml(path, "summary.yaml")

    @staticmethod
    def __valid_library(path: str) -> bool:
        """
        Check if the library pointed to by 'path' is a valid C8 library.

        Parameters
        ----------
        path: str
            The path to the directory containing the library.

        Returns
        -------
        bool:
            True if this is a valid C8 library.

        """

        # check that the config file exists
        if not op.exists(op.join(path, "config.yaml")):
            return False

        # the config file exists, we load it
        config = Library.load_config(path)
        if config is None:
            return False

        # and check that the config's "writer" key is correct
        return config.get("creator") == "CORSIKA8"

    @staticmethod
    def __build_outputs(
        path: str, dirs: Union[list, None]
    ) -> Dict[str, outputs.Output]:
        """
        Build the outputs contained in this library.

        This will print a warning message if a particular
        output is invalid but will continue to load additional
        outputs afterwards.

        Parameters
        ----------
        path: str
            The path to the directory containing this library.

        dirs: list[str]
            List of directory names that will be read in.
            If None, will attempt to find the directory names

        Returns
        -------
        Dict[str, Output]:
            A dictionary mapping names to initialized outputs.
        """

        if dirs is None:
            # if not supplied, get a list of the subdirectories in the library
            _, dirs, _ = next(os.walk(path))

        # this is the dictionary where we store our components
        components: Dict[str, Any] = {}

        # loop over the subdirectories
        for subdir in dirs:
            # read the config file for this output
            config = Library.load_config(op.join(path, subdir))
            if config is None:
                msg = "Could not find a configuration file in"
                msg += f" {op.join(path, subdir)}. Skipping sub-directory"
                c8_logger.warning(msg)
                continue

            # the name keyword is our unique identifier
            name = config.get("name")

            # get the "type" keyword to identify this component
            out_type = config.get("type")

            # if `out_type` was None, this is an invalid output
            if out_type is None or name is None:
                msg = (
                    f"'{subdir}' does not contain a valid config."
                    "Missing 'type' or 'name' keyword."
                )
                c8_logger.warning(msg)
                continue  # skip to the next output, don't error

            # we now have a valid component type, get the corresponding
            # type from the proccesses subdirectory
            try:
                # instantiate the output and store it in our dict
                component = getattr(outputs, out_type)(op.join(path, subdir))

                # check if the read failed
                if not component.is_good():
                    msg = (
                        f"'{name}' encountered an error while reading. "
                        "This process will be not be loaded."
                    )
                    c8_logger.warning(msg)
                else:
                    components[name] = component

            except AttributeError as e:
                msg = (
                    f"Unable to instantiate an instance of '{out_type}' "
                    f"for a process called '{name}'. Skipping '{subdir}'"
                )
                c8_logger.warning(msg)
                c8_logger.warning(e)
                continue  # skip to the next output, don't error

        # and we are done building - return the constructed outputs
        return components
