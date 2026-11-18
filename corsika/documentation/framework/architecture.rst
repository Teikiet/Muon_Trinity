Framework Architecture
======================

CORSIKA 8 is a particle-cascade framework: it transports particles through a
configured environment while applying physics processes that either modify a
particle or produce its secondaries.  This separation is intentional.  A
simulation can combine different media, tracking implementations, particle
stack extensions, and physics modules without changing the central cascade
algorithm.

This architecture follows the design described in the `CORSIKA 8 framework
paper <https://doi.org/10.1007/s41781-018-0013-0>`_.  The paper is useful
background for the design goals—modularity, reproducibility, and efficient
use of compute resources—whereas this guide documents the interfaces exposed
by the current code.

Core simulation loop
--------------------

The central object is :dox:cpp:class:`corsika::Cascade`.  It owns no physics
configuration itself; instead, it is constructed with references to an
environment, a tracking algorithm, a process sequence, an output object, and
a particle stack.  At a high level, the execution flow is:

#. Seed the stack with one or more primary particles and assign their current
   geometry nodes.
#. Select a particle from the stack.
#. Build its trajectory with the configured tracking algorithm.
#. Determine the largest permitted step from the active continuous processes
   and geometry boundaries.
#. Transport the particle and execute continuous processes along the accepted
   step.
#. Sample and execute the applicable discrete event, such as a decay or an
   interaction.  This can add secondary particles to the stack.
#. Run processes that consume or annotate secondaries, write output, and
   continue until the stack is empty.

In pseudocode, the controlling idea is:

.. code-block:: text

   while the stack is not empty:
       particle = select a stack entry
       trajectory = tracking.getTrack(particle)
       step = minimum limit from processes and geometry
       propagate particle over step
       apply continuous processes
       sample and apply a discrete process, if one occurs
       add and process secondaries

The exact selection order and the available stack fields are configuration
choices.  Do not assume that an individual particle is a persistent,
independently allocated C++ object: a particle handle refers to data managed
by the stack.  See :doc:`particle_stack` for the data model.

Processes and step limits
-------------------------

The :dox:cpp:class:`corsika::process::ProcessSequence` combines compatible
process objects at compile time.  Processes may participate in several stages
of a simulation step:

* **Continuous processes** constrain the maximum track length and update a
  particle along a trajectory.  Examples include energy loss, observation
  surfaces, and quantities accumulated over the path.
* **Discrete processes** describe probabilistic events that can remove or
  transform a particle and create secondaries.  Interactions and decays are
  handled in this category.
* **Secondaries processes** receive produced particles, for example to record,
  filter, or thin them before their later transport.
* **Stack and boundary processes** provide work that is associated with the
  whole stack or with a transition between geometry volumes.
* **Cascade-equation processes** can operate on a stack through a numerical
  cascade-equation representation when supplied by the configured modules.

The step length is a shared accuracy and performance control.  Each active
continuous process can impose a maximum, and the cascade uses the smallest
allowed value.  A geometry boundary also terminates a step: a process must not
silently apply material properties from one volume to a trajectory segment in
another volume.  Consequently, process authors should return physically
meaningful limits and should not rely on another module to repair an overly
long step.

Compile-time composition and run-time data
------------------------------------------

The framework uses templates and static process composition in performance-
critical paths.  Fundamental choices—such as the stack type, tracking class,
environment model, and process list—are made by the C++ application and
therefore require recompilation when changed.  This permits the compiler to
optimize the selected combination and avoids dynamic dispatch in the inner
transport loop.

The values supplied to these components, such as primary energy, atmosphere
tables, observation geometry, and model parameters, are ordinary run-time
data.  Keep this distinction in mind when designing an application: choose
the *kinds* of components in C++, then configure their *values* through the
application's supported inputs.

Main building blocks
--------------------

* :doc:`environment` describes the volume hierarchy, coordinate system, and
  material models.
* :doc:`tracking` turns the current particle state into a trajectory and
  identifies crossings of volume boundaries.
* :doc:`particle_stack` stores particles and optional per-particle extensions.
* :doc:`modules` documents the framework process interfaces and available
  module API reference.
* :doc:`output` describes the structured output library written by modules.

When extending CORSIKA 8, prefer adding a narrowly scoped process or model to
one of these building blocks over embedding physics-specific control flow in
an application loop.  This keeps the simulation composable and makes it
possible for other applications to reuse the extension.