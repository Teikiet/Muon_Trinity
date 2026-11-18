Cherenkov Light
================

Cherenkov light is electromagnetic radiation emitted by a charged particle moving
through a transparent medium faster than light *can propagate in that medium*. It
never means that the particle exceeds the vacuum speed of light :math:`c`. Light
has the phase velocity :math:`c/n` in a material with refractive index :math:`n`,
which is slightly below :math:`c` in air.

For an extensive air shower, relativistic electrons and positrons are the main
sources. Their light can be observed by imaging atmospheric Cherenkov telescopes.
This page introduces the physics, the steps an optical-photon simulation needs to
perform, and the refractive-index facilities that are currently available in
CORSIKA 8.

.. important::

   CORSIKA 8 currently provides atmospheric refractive-index models, but it does
   **not** yet provide a Cherenkov optical-photon emission, transport, or detector
   module. In particular, the legacy CORSIKA-7 Cherenkov input cards in
   ``modules/data/CHERENKOV`` are data files, not CORSIKA-8 configuration.
   The emission and transport workflow below describes the physical method a
   future optical module must implement; it is not an available CORSIKA-8 API.

When is light emitted?
----------------------

Let :math:`v` be the particle speed and define :math:`\beta=v/c`. Cherenkov light
is possible only if

.. math::

   \beta n(\lambda, \mathbf{x}) > 1.

Here :math:`n` can depend on both wavelength :math:`\lambda` and position
:math:`\mathbf{x}`. The condition is a threshold: below it the Cherenkov yield is
zero. For a particle above threshold, constructive interference produces a cone of
light around its instantaneous flight direction. The opening angle is

.. math::

   \cos\theta_\mathrm{C} = \frac{1}{\beta n},
   \qquad
   \sin^2\theta_\mathrm{C} = 1 - \frac{1}{\beta^2 n^2}.

The apparent analogy is a sonic boom: the particle leaves a coherent wavefront
behind it. A useful distinction is that Cherenkov light is directional and prompt,
whereas air fluorescence is approximately isotropic and follows energy deposition.

In the atmosphere, :math:`n-1` decreases as the air becomes thinner with altitude.
Consequently, the threshold, the cone angle, and the light yield vary along a
shower trajectory. For an ultrarelativistic particle in air, the cone is narrow;
this is why a correct local particle direction and local refractive index matter.

How many photons?
-----------------

For a particle with charge :math:`ze`, the Frank--Tamm formula gives the spectral
yield in the usual approximation:

.. math::

   \frac{\mathrm{d}^2N}{\mathrm{d}x\,\mathrm{d}\lambda}
   = \frac{2\pi\alpha z^2}{\lambda^2}
     \left(1 - \frac{1}{\beta^2 n^2(\lambda)}\right),

where :math:`\alpha` is the fine-structure constant. The expression is used only
where the threshold condition is satisfied. Its :math:`1/\lambda^2` dependence
means that more short-wavelength photons are produced, which contributes to the
familiar blue appearance.

For a constant refractive index over a selected wavelength interval
:math:`[\lambda_\min, \lambda_\max]`, integrating the formula gives

.. math::

   \frac{\mathrm{d}N}{\mathrm{d}x}
   = 2\pi\alpha z^2
     \left(1 - \frac{1}{\beta^2 n^2}\right)
     \left(\frac{1}{\lambda_\min} - \frac{1}{\lambda_\max}\right).

A real simulation should select a finite band compatible with atmospheric
transmission and detector sensitivity. It must also account for wavelength
dependence (dispersion) when that precision is required; a wavelength-independent
:math:`n` is a controlled approximation, not a universal property of air.

Physical emission workflow
--------------------------

The following method, derived from the photon-generation approach described in
the accompanying optical-photon work, is the practical translation of the equations
above:

#. **Obtain local track data.** At a sufficiently short part of a charged
   particle's path, use its charge, speed, position, and *instantaneous* direction
   together with the local refractive index.
#. **Test the threshold and calculate the yield.** If :math:`\beta n \leq 1`, emit
   no Cherenkov photons. Otherwise evaluate the Frank--Tamm yield in the chosen
   wavelength band and sample or weight the corresponding photon population.
#. **Sample the emission cone.** Draw an azimuth :math:`\phi` uniformly in
   :math:`[0,2\pi)`. In coordinates whose :math:`z` axis is the particle direction,
   a photon direction is

   .. math::

      \hat{\mathbf{k}} =
      (\sin\theta_\mathrm{C}\cos\phi,
       \sin\theta_\mathrm{C}\sin\phi,
       \cos\theta_\mathrm{C}).

   Rotate this direction into the simulation coordinate system.
#. **Sample the emission point.** Place the photon at a point of
   the charged particle track segment, this can be randomly sampled or tailored to the refractive index.
#. **Transport and detect.** Propagate photons through the atmosphere, applying
   the selected absorption, scattering, refraction, and detector-response model.
   These choices are experiment-specific and must be recorded with the simulation
   configuration.

Why tracks must be subdivided
------------------------------

Particle tracking commonly represents a step by its start and end point. That can
be sufficient for the cascade particles themselves, but it can be too coarse for
Cherenkov light: magnetic deflection and multiple scattering change the local
particle direction. Since the cone axis is the local direction, generating every
photon from one long straight segment biases photon origins and arrival directions.

An optical implementation should therefore substep only for photon emission. The
substeps approximate the curved or scattered path while preserving the original
particle-transport step where possible. The required substep size depends on the
angular and timing resolution demanded by the observation geometry. Reducing every
cascade step indiscriminately can greatly increase runtime, so adaptive substepping
is preferable when a reliable path interpolation is available.

Available refractivity support
------------------------------

CORSIKA 8 environments can expose a position-dependent refractive index through
``media::IRefractiveIndexModel`` and its
``getRefractiveIndex(Point const&)`` method. This interface accepts position only;
it currently has no wavelength argument. It can support local threshold and
angle calculations in a future nondispersive optical implementation, but it cannot
by itself model dispersion.

The current models are:

* ``media::UniformRefractiveIndex``: a constant :math:`n`, useful for controlled
  studies and tests.
* ``media::GladstoneDaleRefractiveIndex``: a density-scaled atmospheric model,

  .. math::

     n(\mathbf{x}) = 1 + (n_\mathrm{ref}-1)
       \frac{\rho(\mathbf{x})}{\rho_\mathrm{ref}}.

  This is the appropriate starting point when an atmospheric density model is
  already part of the environment.
* ``media::ExponentialRefractiveIndex`` and
  ``media::LinearTabulatedRefractiveIndex``: alternative spatial profiles for
  cases where their respective parameterisation or table is suitable.

See :doc:`../framework/media` for the medium-property API. The
``applications/c8_air_shower.cpp`` application demonstrates composition of an
environment with ``GladstoneDaleRefractiveIndex``; it does not enable optical
Cherenkov photon generation.


Further reading
---------------

* I. Frank and I. Tamm, *Coherent visible radiation of fast electrons passing
  through matter*, Doklady Akademii Nauk SSSR **14** (1937), 109--114.
* D. Baack, *Optical Photon Emission in Extended Airshowers -- Hybrid computing
  in the context of CORSIKA 8*, dissertation (2023). Its discussions of Cherenkov
  radiation, photon generation, and emission-track substepping motivated this
  overview.
