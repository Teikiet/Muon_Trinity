# CORSIKA 8 - Python Library

The python library to read the CORSIKA 8 output can be installed using the `pip` command. It is recommended that both CORSIKA 8 and this library are installed within a virtual environment. To learn how to do this, see the virtual-environment installation instructions [here](../README.md).

## For general users

Install the library directly from the main branch (you will not be able to make any changes to them).
**Note:** if you are NOT using a virtual environment, you may want to use the `pip install --user` instead.

``` shell
pip install 'git+https://gitlab.iap.kit.edu/AirShowerPhysics/corsika.git#subdirectory=python'
```

## For developers

If you are developing the CORSIKA 8 framework or, more generally, want to install the code directly from a cloned repo.
**Note:** if you are NOT using a virtual environment, you may want to use the `pip install --user ...` instead.

``` shell
cd ./path/to/corsika/python
pip install -e .[test,examples]
```

## Examples

Using this library you can directly read the output from corsika simulations by using the `Library` class.
For example, after running a corsika shower `./bin/corsika -E 1e4 -f shower_output`, the output will be made in a directory called `shower_output` and can be read using:

``` python
from corsika8.io import Library

# Load the library
lib = Library("shower_output")

# Print out meta-data about the shower
print("Library configuration:")
print(lib.config)

# Print out the names of the available objects to load
# These reflect the names of the sub-dirs in the output dir
print("Library sub-directories:")
print(lib.names)

# Get the longitudinal profile as a pandas data-frame
profile = lib.get("profile").astype("pandas")

```

For more advanced examples see the `corsika/python/examples` directory.
Examples can be run like this:

```shell
python examples/shower_profile.py --input-dir <path-to-C8-output>
```
