Applications and Examples
=========================

Applications and examples provide the most direct starting points for running
or extending CORSIKA 8 after it has been installed.

Applications
------------

The `applications directory <https://gitlab.iap.kit.edu/AirShowerPhysics/corsika/tree/main/applications>`_
contains physics-complete programs intended for standard astroparticle-physics
simulations. In particular,
`c8_air_shower <https://gitlab.iap.kit.edu/AirShowerPhysics/corsika/blob/main/applications/c8_air_shower.cpp>`_
simulates air showers in a curved atmosphere. Build and invocation instructions
are maintained in the repository README.

Examples
--------

The `examples directory <https://gitlab.iap.kit.edu/AirShowerPhysics/corsika/tree/main/examples>`_
contains smaller programs for learning and extension work:

* ``cascade_examples`` configure full cascade simulations for particular use
  cases.
* ``framework_examples`` demonstrate use of the framework's internal APIs.
* ``physics_examples`` exercise individual particle interactions and cascades.

For help creating a new process, see :doc:`../framework/howto_create_module`.
Simulation output can be inspected with :doc:`python_library`.