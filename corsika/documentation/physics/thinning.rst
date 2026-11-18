Thinning
========

Thinning reduces the number of tracked particles in extensive air-shower
simulations while retaining a statistically representative cascade.

In CORSIKA 8, ``EMThinning`` is a secondaries process. It acts after a physics
process has created secondaries and before those particles are transported. The
current implementation applies only to electromagnetic ``1 -> 2`` splittings:
an electromagnetic primary must produce two electromagnetic secondaries. It is
therefore not a general thinning algorithm for every particle or interaction
type in a cascade.

Below its configured energy threshold, the process retains a statistically
selected subset of secondaries and assigns weights to preserve expectation
values. ``maxWeight`` limits those weights. When a Hillas selection would
exceed that limit, the implementation uses statistical thinning, retaining the
secondaries independently with probabilities chosen to respect the maximum
weight. Thinning reduces CPU time and output volume at the cost of statistical
fluctuations; weighted observables must be used consistently in downstream
analysis.

Standard air-shower application
-------------------------------

``c8_air_shower`` exposes the following controls:

``--emthin``
	Fraction of the primary energy below which electromagnetic thinning starts.
	The default is ``1e-6``.

``--max-weight``
	Maximum electromagnetic-particle weight. The default ``0`` selects one half
	of Kobal's optimum weight.

``--multithin``
	Keeps thinned particles with zero weight, which is useful when a workflow
	needs to retain their history or bookkeeping information.

For the CORSIKA 8 radio-validation study, the compared air showers used an
electromagnetic thinning threshold of ``1e-6`` of the primary energy and a
maximum weight of ``50``. These values support that particular comparison; they
are not a substitute for convergence studies of the energy, geometry,
observables, and detector response of a new simulation.

Further reading
---------------

* `J. A. M. H. J. Hillas, original ICRC 1981 thinning reference <https://gitlab.iap.kit.edu/AirShowerPhysics/corsika/-/wikis/uploads/a178f66e1b5c54b1e3d23b6866f63ebd/_JK_Hillas-Symposium.pdf>`_.
* J. M. Alameddine et al., *Simulating radio emission from particle cascades
  with CORSIKA 8*, Astroparticle Physics **166** (2025) 103072,
  `doi:10.1016/j.astropartphys.2024.103072
  <https://doi.org/10.1016/j.astropartphys.2024.103072>`_.