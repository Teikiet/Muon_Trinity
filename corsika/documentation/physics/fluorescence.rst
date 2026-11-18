Fluorescence Light
==================

Air fluorescence is light emitted when atmospheric molecules excited by the
energy deposited by a charged-particle cascade return to lower-energy states.
Unlike Cherenkov light, which is emitted on a cone around the particle direction,
fluorescence is approximately isotropic. Its yield follows local energy
deposition, so a fluorescence simulation must account for the atmosphere at the
emission point as well as the selected wavelength band and detector response.

.. important::

   CORSIKA 8 does not currently provide a fluorescence photon-emission,
   transport, or detector module. This page documents the physical requirements
   of a future implementation; it does not describe an available CORSIKA 8 API.

Requirements for an optical implementation
------------------------------------------

An implementation must obtain the energy deposited along a sufficiently local
particle-track segment, convert that deposit to a photon yield using an
atmospheric fluorescence-yield model, sample isotropic emission directions, and
transport the photons to a detector. The yield and transport depend on details
such as pressure, temperature, humidity, wavelength-dependent transmission, and
instrument response. These inputs must therefore be part of the recorded
simulation configuration.

The optical-photon work that motivated the CORSIKA 8 overview describes
fluorescence generation together with Cherenkov-photon generation and
substepping. CORSIKA 8 currently exposes atmospheric refractive-index support,
but that support alone is insufficient for fluorescence emission or optical
transport. See :doc:`cherenkov` for the corresponding Cherenkov status and
shared optical-transport considerations.

Further reading
---------------

* D. Baack, *Optical Photon Emission in Extended Airshowers -- Hybrid computing
  in the context of CORSIKA 8*, dissertation (2023).