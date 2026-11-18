Dependencies
============

CORSIKA 8 keeps the number of external dependencies as small as practical. The
maintained installation commands are documented in the repository README; this
page gives an overview of the main packages and why they are used.

The only dependency that must be available before bootstrapping the build is
``python3``. Conan is installed through ``pip`` and then used to obtain the
remaining packages on supported systems.

Core Packages
-------------

``eigen3``
   Header-only package for linear algebra and geometry. Avoid implicit Eigen
   include-directory dependencies such as ``Eigen3_INCLUDE``; use the explicit
   CMake target ``CORSIKAeigen3`` instead.

``boost``
   Header-only Boost components are used where needed, including
   ``boost-histogram``. Avoid implicit include-directory dependencies such as
   ``Boost_INCLUDE``; use the explicit CMake target ``CORSIKAboost`` instead.

``phys/Units``
   Header-only physics-units support. CORSIKA 8 provides its own version as
   part of the framework.

``catch2``
   Header-only unit-test framework used by the test suite.

Physics and Optional Packages
-----------------------------

``pythia8``
   Event and decay generator. CORSIKA 8 supports a system-level installation
   and provides a tarball-backed version through ``modules``.

``conex``
   Cascade-equation solver for air showers. CORSIKA 8 supports a system-level
   installation and a git-submodule-backed version through ``modules``.

``proposal``
   Muon and electromagnetic-physics package provided through ``modules``.

Development Checks
------------------

When changing build logic, make dependencies explicit in CMake targets instead
of relying on global include paths. Useful checks include:

* ``grep -Ri "Eigen3_INCLUDE" .``
* ``grep -R "Boost_INCLUDE" .``
* ``grep -R "#include <boost/" .``
