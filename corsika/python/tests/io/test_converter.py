import numpy as np
import pyarrow as pa
import pytest

from corsika8.io import converters


def test_convert_to_numpy_basic_types():
    # Create a PyArrow table with supported types
    data = {
        "int_col": pa.array([1, 2, 3], type=pa.int32()),
        "float_col": pa.array([1.1, 2.2, 3.3], type=pa.float64()),
        "uint_col": pa.array([10, 20, 30], type=pa.uint8()),
    }
    table = pa.table(data)

    # Convert to NumPy structured array
    result = converters.convert_to_numpy(table)

    # Check dtype
    expected_dtype = np.dtype(
        [
            ("int_col", "int32"),
            ("float_col", "float64"),
            ("uint_col", "uint8"),
        ]
    )
    assert result.dtype == expected_dtype

    # Check values
    assert np.array_equal(result["int_col"], np.array([1, 2, 3], dtype="int32"))
    assert np.allclose(result["float_col"], [1.1, 2.2, 3.3])
    assert np.array_equal(result["uint_col"], np.array([10, 20, 30], dtype="uint8"))


def test_convert_to_numpy_unsupported_type():
    # Create a table with an unsupported type (e.g., string)
    data = {"str_col": pa.array(["a", "b", "c"], type=pa.string())}
    table = pa.table(data)

    # Expect a TypeError
    with pytest.raises(TypeError, match="Unsupported PyArrow type"):
        converters.convert_to_numpy(table)
