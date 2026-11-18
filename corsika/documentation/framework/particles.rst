Particle Storage and Stack
==========================

Particles in memory are stored in a Stack. The Stack handles the memory management and access to particle data
properties. The stack can be extended with additional fields that then become part of the particle object.

The standard Stack object in CORSIKA 8 is the NuclearStackExtension, further extended with internal
fields like `geometry node` and `weight` or `history`. But the important part is the `NuclearStackExtension`.

There are several ways to create particles:

* Novel particles on the stack:

  * ``stack::addParticle(particle_data_type const&)``
  * ``stack::addParticle(nuclear_particle_data_type const&)``
  * ``stack::addParticle(particle_data_momentum_type const&)``
  * ``stack::addParticle(nuclear_particle_data_momentum_type const&)``
* Secondary particles:

  * ``stack::addSecondary(parent const&, particle_data_type const&)``
  * ``stack::addSecondary(parent const&, nuclear_particle_data_type const&)``
  * ``stack::addSecondary(parent const&, particle_data_momentum_type const&)``
  * ``stack::addSecondary(parent const&, nuclear_particle_data_momentum_type const&)``
  * Or directly:

    * ``particle::addSecondary(particle_data_type const&)``
    * ``particle::addSecondary(nuclear_particle_data_type const&)``
    * ``particle::addSecondary(particle_data_momentum_type const&)``
    * ``particle::addSecondary(nuclear_particle_data_momentum_type const&)``

The content of these method parameters are:

* ``particle_data_type``: PID, kinetic energy, direction, position, and time.
* ``particle_data_momentum_type``: PID, momentum, position, and time.
* ``nuclear_particle_data_type``: PID, kinetic energy, direction, position,
  time, mass number, and atomic number.
* ``nuclear_particle_data_momentum_type``: PID, momentum, position, time,
  mass number, and atomic number.
