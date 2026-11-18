# CORSIKA 8 Framework for Particle Cascades in Astroparticle Physics

The purpose of CORSIKA 8 is to simulate any particle cascades in
astroparticle physics or astrophysical context. Considerable emphasis has been
put on modularity, flexibility, completeness, validation and
correctness. To boost computational efficiency, different techniques
are provided, like thinning or cascade equations. The aim is for
CORSIKA 8 to remain the most comprehensive framework for simulating
particle cascades with stochastic and continuous processes.

The software makes extensive use of static design patterns and
compiler optimization. Thus, the most fundamental configuration
decisions of the user must be performed at compile time. At run time,
model parameters can still be changed.

CORSIKA 8 is by default released under the BSD 3-Clause License. See [license
file](https://gitlab.iap.kit.edu/AirShowerPhysics/corsika/blob/main/LICENSE)
which is part of every release and the source code.

If you use, or want to refer to, CORSIKA 8 please cite ["CORSIKA 8: A general framework for particle cascade simulations", Astropart. Phys. 182 (2026) 103265](https://doi.org/10.1016/j.astropartphys.2026.103265) 
For simulations including radio emission, please also cite 
["Simulating radio emission from particle cascades with CORSIKA 8", Astropart. Phys. 166 (2025)
103072](https://doi.org/10.1016/j.astropartphys.2024.103072).

We kindly ask (and require) any relevant improvement or addition to be offered or
contributed to the main CORSIKA 8 repository for the benefit of the
whole community.

CORSIKA 8 makes use of various third-party code, in particular interaction
models. Please check the [using and collaborating
agreement](https://gitlab.iap.kit.edu/AirShowerPhysics/corsika/blob/main/usage_and_collaboration_agreement.md)
for further information on this topic.

If you plan to contribute to CORSIKA 8, please check the guidelines outlined here:
[coding
guidelines](https://gitlab.iap.kit.edu/AirShowerPhysics/corsika/blob/main/documentation/code/contributing.rst). Code
that fails the review by the CORSIKA 8 author group must be improved
before it can be merged in the official code base. After your code has
been accepted and merged, you become a contributor of the CORSIKA 8
project (code author).

IMPORTANT: Before you contribute, you need to read and agree to the conditions set out in the
[using and collaborating
agreement](https://gitlab.iap.kit.edu/AirShowerPhysics/corsika/blob/main/usage_and_collaboration_agreement.md).
The agreement can be discussed, and eventually improved if necessary.


## Get in contact
  * Join our Matrix space using this [invite link](https://matrix.to/#/#mh7x43sm-o8dk.swi7GmnAS7cgsk9wncSjze:kit.edu). Select and install your preferred Matrix client and follow the instructions.
  * Connect to https://gitlab.iap.kit.edu, register yourself and join the "Air Shower Physics" group. Write to us on Matrix (in the User Questions channel), or directly contact one of the [steering committee members](https://gitlab.iap.kit.edu/AirShowerPhysics/corsika/-/wikis/Steering-Committee) in case there are problems with that.
  * Connect to corsika-devel@lists.kit.edu (self-register at
    https://www.lists.kit.edu/sympa/subscribe/corsika-devel) to get in
    touch with the project.

## Using CORSIKA 8 Container Releases

CORSIKA 8 is distributed not only as source code but also as a pre-built container image for use with Apptainer.
This container release provides a fully configured and portable runtime environment that includes all necessary dependencies and is ideal for users who want to run CORSIKA 8 without building from source.
For development or contributing to the framework, please refer to the installation instructions below.

You can download the latest container release on [zenodo](https://zenodo.org/records/18195386) or from the [GitLab releases page](https://gitlab.iap.kit.edu/AirShowerPhysics/corsika/-/releases). To use the container, ensure that [Apptainer] is installed on your system. To execute CORSIKA 8 inside the container, use the `apptainer run` command, specifying the `.sif` file. On most systems, you can also execute the container directly. In both cases, you will execute the binary `/opt/corsika8/install/bin/c8_air_shower`.
An example is shown below:

```shell
apptainer run /path/to/corsika8.sif -E 1e3 -p 2212 -f output/path/test
```

This is equivalent to using the container as an executable
```shell
/path/to/corsika8.sif -E 1e3 -p 2212 -f output/path/test
```

If you run the container outside your home directory, you might get failures which result from problems with file access. This often means that your current working directory is not available in the container. You can request additional mounts when running using, for example:
```shell
apptainer run --bind $(pwd):$(pwd) /path/to/corsika8.sif -E 1e3 -p 2212 -f output/path/test
```
Depending on the directory structure on your installation, you might need other or additional [bind mount options].

You can explore the contents of the container using an interactive shell
```shell
apptainer shell /path/to/corsika8.sif
```
or running individual commands directly, for example
```shell
apptainer exec /path/to/corsika8.sif ls /opt/corsika8/install/bin
```

[Apptainer]: https://apptainer.org
[bind mount options]: https://apptainer.org/docs/user/main/bind_paths_and_mounts.html#user-defined-bind-paths

## Installation

CORSIKA 8 is tested regularly at least on `gcc11.0.0` and `clang-14.0.0`.

### Prerequisites

You will also need:

- Python 3 (supported versions are Python >= 3.6), with pip
- cmake > 3.4
- git
- g++, gfortran, binutils, make
- rsync
- tar
- optional: FLUKA (see below)


On a bare Ubuntu machine, just add:
``` shell
sudo apt-get install python3 python3-pip cmake g++ gfortran git doxygen graphviz rsync tar
```

### Creating a virtual environment and Conan

It is recommended that you install CORSIKA 8 and its dependencies within a python3 virtual environment.
To do so, you can run the following.
``` shell
# Create the environment using your native python3 binary
python3 -m venv /path/to/new/virtual/environment/corsika-8
# Load the environment (should be run each time you open a new terminal)
source /path/to/new/virtual/environment/corsika-8/bin/activate

```

You will need to load the environment each time that you open a new terminal.

CORSIKA 8 uses the [conan](https://conan.io/) package manager to
manage our dependencies. Currently, version 2.0.0 or higher is required.
**Note**: if you are NOT using a virtual environment, you may want to use the `pip install --user` flag.

``` shell
pip install conan particle==0.25.1 numpy
```

### Enabling FLUKA support

For legal reasons, we do not distribute or bundle FLUKA together with the CORSIKA 8 source core.
As FLUKA is the standard low-energy hadronic interaction model for CORSIKA 8, you have to download
it separately from (http://www.fluka.eu/), which requires registering there as FLUKA user.
The following should be done *before* compiling CORSIKA 8:

 1. Note your system's version of gfortran (`gfortran --version`) and glibc (`ldd --version`)
 2. Download the FLUKA __binary__, ensuring that it matches the versions you found above
 3. Download the FLUKA __data file__ (will be named something similar to __fluka20xy.z-data.tar.gz__).
 4. Un-tar the files that you downloaded using `tar -xf <filename>`
 5. Set environment variables `export FLUFOR=gfortran` and `export FLUPRO=<path to where you unzipped the files>`. Note that the `FLUPRO` directory should contain __libflukahp.a__. Both of these variables will have to be set every time you open a new terminal.
 6. Go to the `FLUPRO` directory and run `make`. This will compile an exe, __flukahp__, in your current directory.
 7. Follow the normal steps to compile CORSIKA 8 (see below).

When you later install CORSIKA 8, you should see a message during the __cmake__ step indicating the FLUKA was correctly found.

``` shell
libflukahp.a found in directory <some location here> via FLUPRO environment variable
FLUKA support is enabled.
```

### Compiling CORSIKA 8

Once Conan is installed and FLUKA is provided, follow these steps to download and install CORSIKA 8:

``` shell
cd ./top/directory/for/corsika/installation
git clone --recursive git@gitlab.iap.kit.edu:AirShowerPhysics/corsika.git
# Or for https: git clone --recursive https://gitlab.iap.kit.edu/AirShowerPhysics/corsika.git
mkdir corsika-build
cd corsika-build
../corsika/conan-install.sh --source-directory ../corsika --release-with-debug
# conan-install.sh takes required options from command line to install dependencies for 'Debug', 'Release' and 'RelWithDebInfo' builds.
../corsika/corsika-cmake.sh -c "-DCMAKE_BUILD_TYPE="RelWithDebInfo" -DWITH_FLUKA=ON -DCMAKE_INSTALL_PREFIX=../corsika-install"
make -j4  # The number should match the number of available cores on your machine
make install
```

## Running Unit Tests

To run the unit tests, do the following.

```shell
cd ./corsika-build
ctest -j4  # The number should match the number of available cores on your machine
```

## Running applications and examples

### Standard applications

Applications for standard use-cases are located in the `applications` directory.
These are example scripts that can be used directly or slightly modified for your use case.
See [applications/README.md] for more.
The applications are compiled automatically after running `make` and will appear in your `corsika-build/bin` directory.
After running `make install` the binaries will also be copied into your `corsika-install/bin` directory.


For example, from inside your `corsika-install/bin` directory, run
```shell
c8_air_shower --pdg 2212 -E 1e5 -f my_shower
```
This will run a vertical 100 TeV proton shower and will create and put the output into `./my_shower`.


### Building the examples

Unlike the applications, the examples must be compiled as a second step.
From your top corsika directory, (the one that includes `corsika-build` and `corsika-install`) run
```shell
export CONAN_DEPENDENCIES=$PWD/corsika-install/lib/cmake/dependencies
cmake -DCMAKE_TOOLCHAIN_FILE=${CONAN_DEPENDENCIES}/conan_toolchain.cmake -DCMAKE_PREFIX_PATH=${CONAN_DEPENDENCIES} -DCMAKE_POLICY_DEFAULT_CMP0091=NEW -DCMAKE_BUILD_TYPE=RelWithDebInfo -Dcorsika_DIR=$PWD/corsika-build -DWITH_FLUKA=ON -S $PWD/corsika/examples -B $PWD/corsika-build-examples
cd corsika-build-examples
make -j4 # The number should match the number of available cores on your machine
```

You can run the examples by using the binaries in `corsika-build-examples/bin/`.
For example:
```shell
corsika-build-examples/bin/known_particles
```
This will print out all of the particles that are known by CORSIKA.


### Generating documentation

The ``docs`` target generates both the Doxygen API reference and the Sphinx user
manual. It requires Doxygen, Graphviz, and Sphinx. If you work with the docker
corsika/devel containers, these are already included.
Otherwise, e.g. on Ubuntu machines, do:
```shell
sudo apt-get install doxygen graphviz python3-sphinx
```
For the complete Sphinx dependency set, install
``documentation/rtfd-requirements.txt`` in your Python environment.

Switch to the `corsika-build` directory and do
```shell
make docs
make install
```
Open the installed user manual with Firefox:
```shell
firefox ../corsika-install/share/corsika/doc/manual/index.html
```
The Doxygen API reference is available at
`../corsika-install/share/corsika/doc/html/index.html`.
