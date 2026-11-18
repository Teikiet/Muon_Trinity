Geometry and Environment
========================

An environment describes the space in which a cascade is transported.  It
combines a global coordinate system with a hierarchy of geometry volumes whose
media provide the material information needed by physics processes.  This lets
one simulation describe, for example, an atmosphere together with a denser
observation medium instead of assuming that every particle remains in air.

The high-level container is :dox:cpp:class:`corsika::media::Environment`.  It
owns a :dox:cpp:class:`corsika::media::Universe`, an unbounded root of the
volume hierarchy, and exposes the environment coordinate system.  Child
volumes are represented by
:dox:cpp:class:`corsika::media::VolumeTreeNode` objects.  Use
``Environment::createNode`` to create nodes with the environment's configured
medium interface.

Volumes, media, and boundaries
------------------------------

Each particle has a current geometry node.  Tracking uses that node and the
particle trajectory to find the next relevant boundary.  A step is ended at
that boundary so that the segment is evaluated with one well-defined medium.
After a crossing, boundary processing updates the particle's node before
subsequent physics is evaluated.

Material models provide quantities such as mass density and nuclear
composition.  Physics modules use these properties to calculate effects such
as interaction probabilities and energy loss.  Properties that depend on
position should be implemented by an inhomogeneous model rather than by
changing a particle's material data directly.  The media API reference is
available from :doc:`media`.

The supplied media include homogeneous and interpolated density models, as
well as helpers for layered spherical atmospheres.  These are building blocks,
not a limitation of the framework: applications may construct a different
volume tree or provide another compatible medium model.

Coordinate systems
------------------

Geometry objects carry their coordinate-system context.  The environment
coordinate system is the natural global frame for volume definitions and
particle positions.  Formulate geometry operations in terms of the framework
point, vector, and trajectory types rather than converting components by hand;
this preserves the coordinate-system and unit checks implemented by the
geometry layer.

For air-shower applications, a global Earth-centred frame and a shower-aligned
frame serve different purposes.  The former is convenient for large-scale
environment geometry, while the latter is useful for longitudinal quantities
and detector coordinates.  Transform explicitly at the interface where a
shower-specific quantity is required; do not redefine the environment's global
frame during a run.

See :doc:`tracking` for how trajectories traverse this hierarchy and
:doc:`architecture` for the role of boundaries in a cascade step.
