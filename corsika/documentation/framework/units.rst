Physics Units
=============

CORSIKA 8 represents dimensionful quantities with strongly typed physical
units.  Length, time, energy, momentum, density, and related dimensions are
part of the C++ type, not only a convention documented beside a variable.  The
compiler can therefore reject invalid expressions such as adding a length to a
time or assigning an inverse length where a frequency is expected.

Why units matter
----------------

Particle-cascade simulations combine interaction models, geometry, detector
descriptions, and tabulated environmental data.  These components commonly
use different natural scales.  A typed unit system makes conversions explicit
at their interface and prevents a silent mismatch from propagating through a
simulation.

With compiler optimization enabled, dimensional typing does not imply that a
quantity carries run-time metadata: the stored value remains a numerical
quantity in the framework's internal unit system.  Use the unit-aware aliases
and literals from ``corsika/framework/core/PhysicalUnits.hpp`` rather than
passing bare ``double`` values through framework interfaces.

Practical rules
---------------

* Keep dimensional values in their typed form for as long as possible.
* Convert at external boundaries only: parsing a text table, calling an
	external generator, or writing a user-selected output representation.
* Include units in output metadata when serialising a raw numeric value.
* Prefer framework geometry types for positions, vectors, and trajectories;
	they preserve both physical dimensions and coordinate-system information.

The unit system is an important part of the framework frontier described in
:doc:`architecture`: adapters for external packages are responsible for making
their unit conventions explicit before data enters the CORSIKA 8 core.
