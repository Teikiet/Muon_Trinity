Interaction Models
==================

An interaction model is a discrete process in the CORSIKA 8 cascade. The model
generates secondaries when a particle interacts, while the framework continues
to handle stack management, tracking, boundaries, and subsequent processes.
Applications choose their process combination explicitly; model availability is
therefore configuration-dependent rather than a single global CORSIKA setting.

Standard air-shower application
-------------------------------

``c8_air_shower`` selects a high-energy hadronic model with
``--hadronModel``. Its supported choices are:

* ``SIBYLL-2.3d`` (the default)
* ``QGSJet-II.04``
* ``QGSJet-III``
* ``EPOS-LHC``
* ``EPOS-LHC-R``
* ``Pythia8``

The application switches between its high- and low-energy hadronic processes at
``--hadronModelTransitionEnergy``, whose default is :math:`10^{1.9}` GeV
(approximately 79.4 GeV). Below this energy it uses FLUKA when CORSIKA 8 was
built with ``WITH_FLUKA``; otherwise it uses UrQMD. FLUKA is optional because
it is separately distributed and must be installed before building CORSIKA 8.
See the project README for the required FLUKA setup.

The standard application also combines several specialized processes with the
hadronic models:

* PROPOSAL handles electromagnetic interactions and continuous transport for
  supported charged particles, using the configured low- and high-energy
  hadronic models for its photo-hadronic component.
* SOPHIA handles hadronic photon interactions in the resonance region.
* Pythia8 supplies ordinary decays and primary neutrino interactions; TAUOLA is
  selected for tau decays.

Charm tracking is disabled by default. When it is requested with
``--track-charm``, the standard application supports it only with SIBYLL,
Pythia8, EPOS-LHC-R, or QGSJet-III.

Choosing a model
----------------

The named models have different physics assumptions, particle coverage, and
validity ranges. Select and report both the high-energy model and the
transition energy for a production simulation. A comparison across models is a
physics-systematics study, not a replacement for validating one configuration.
The available model adapters are documented in the API reference; consult each
model's own literature and licensing terms before publication or redistribution.

Further reading
---------------

* The CORSIKA 8 framework and module interfaces are described by
  :doc:`../framework/architecture` and :doc:`../framework/modules`.
* `CORSIKA 8 installation and FLUKA setup <https://gitlab.iap.kit.edu/AirShowerPhysics/corsika/-/blob/main/README.md>`_.
* J. M. Alameddine et al., *Simulating radio emission from particle cascades
  with CORSIKA 8*, Astroparticle Physics **166** (2025) 103072,
  `doi:10.1016/j.astropartphys.2024.103072
  <https://doi.org/10.1016/j.astropartphys.2024.103072>`_.