"""
Converter for the pyarrow data types to numpy ones

(c) Copyright 2020 CORSIKA Project, corsika8@kit.edu

This software is distributed under the terms of the 3-clause BSD license.
See file LICENSE for a full version of the license.
"""

import numpy as np
import pyarrow

# Type conversions for pyarrow data types to numpy ones
# https://arrow.apache.org/docs/python/data.html
# https://numpy.org/doc/stable/reference/arrays.dtypes.html#arrays-dtypes-constructing
PYARROW_TO_NUMPY = {
    pyarrow.int8().id: "int8",
    pyarrow.int16().id: "int16",
    pyarrow.int32().id: "int32",
    pyarrow.int64().id: "int64",
    pyarrow.uint8().id: "uint8",
    pyarrow.uint16().id: "uint16",
    pyarrow.uint32().id: "uint32",
    pyarrow.uint64().id: "uint64",
    pyarrow.float16().id: "float16",
    pyarrow.float32().id: "float32",
    pyarrow.float64().id: "float64",
}


def convert_to_numpy(pyarrow_table: pyarrow.Table) -> np.ndarray:
    """
    Converts a pyarrow Table to a numpy structured array

    Parameters
    ----------
    pyarrow_table: pyarrow.Table
        PyArrow table of any dimension to be sliced

    Returns
    -------
    np.ndarray:
        converted table with the same column labels and data types

    Raises
    ------
    TypeError
        If a column type is not supported.

    """

    column_names = pyarrow_table.column_names
    column_types = []

    for name in column_names:
        arrow_type_id = pyarrow_table[name].type.id
        try:
            numpy_type = PYARROW_TO_NUMPY[arrow_type_id]
        except KeyError:
            raise TypeError(
                "Unsupported PyArrow type for column"
                f" '{name}': {pyarrow_table[name].type}"
            )
        column_types.append(numpy_type)

    dtype = list(zip(column_names, column_types))
    np_table = np.zeros(pyarrow_table.num_rows, dtype=dtype)

    for name in column_names:
        np_table[name] = pyarrow_table[name]

    return np_table
