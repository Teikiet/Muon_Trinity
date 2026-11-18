Tracking
========

Tracking constructs the geometrical trajectory followed by a particle and
locates its intersections with the environment.  The cascade uses that
trajectory to limit a transport step, apply continuous physics processes, and
identify the next volume boundary.  Tracking therefore describes *where* a
particle goes; processes describe *what happens* while it travels.

Tracking has two linked responsibilities:

#. Propose a trajectory from the particle's position and momentum.
#. Intersect that trajectory with the volume hierarchy and report the next
	relevant boundary.

The currently configured tracking type must provide ``getTrack(particle)``;
:dox:cpp:class:`corsika::Cascade` obtains the particle's path through this
interface.  The returned trajectory is then shared with the process sequence.

Choosing a tracking implementation
----------------------------------

``corsika::tracking_line::Tracking`` provides straight-line tracking.  It is
appropriate for neutral particles or a configuration where magnetic deflection
is not modelled.  It intersects straight tracks with the supported geometry
objects and volume nodes.

Leap-frog tracking implementations are available for configurations requiring
stepwise numerical propagation, including curved trajectories in fields.  The
tracking algorithm and trajectory representation must be selected together;
the setup headers make this relationship explicit.  A tracking class should be
chosen before building the cascade because it is a template parameter in the
performance-critical simulation path.

Step construction
-----------------

A cascade step is not chosen by tracking alone.  The process sequence gathers
the maximum lengths requested by continuous processes and the geometry limits
the trajectory at the next crossing.  The smallest resulting limit is used.
Continuous modules can therefore control accuracy for quantities accumulated
along the path—such as energy loss, field deflection, or an optical
calculation—without duplicating boundary logic.

After transport, continuous processes receive the accepted trajectory segment.
Only then can the cascade sample a discrete event such as an interaction or
decay.  This order makes the path and the medium along it unambiguous.  See
:doc:`architecture` for the complete step sequence and :doc:`environment` for
the volume and medium model.

Implementation guidance
-----------------------

When implementing a new tracker, keep intersections numerically robust near
volume boundaries and make the returned track consistent with the particle
state.  Do not bypass the environment's node information or move a particle
across several media in one unreported step.  Physics effects that require
integration along a track normally belong in a continuous process; the tracker
should remain focused on trajectory construction and geometric intersections.


