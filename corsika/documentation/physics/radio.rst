Radio Emission
==============

CORSIKA 8 can calculate the time-domain electric field emitted by a particle
cascade at a configurable collection of observer positions. The radio module is
a continuous cascade process: each accepted particle track contributes to every
observer, and the observer stores the resulting electric-field trace. This makes
it possible to run radio-emission calculations alongside the shower simulation,
including more than one radio formalism for the same shower.

The implementation and its air-shower validation are described by Alameddine
*et al.*, `Simulating radio emission from particle cascades with CORSIKA 8
<https://doi.org/10.1016/j.astropartphys.2024.103072>`_.

Physics and formalisms
----------------------

Radio emission from an air shower is dominated by two mechanisms:

* **Geomagnetic emission** arises from the opposite deflection of electrons and
  positrons in the geomagnetic field. Its polarization follows the transverse
  current, approximately along :math:`\mathbf{v} \times \mathbf{B}`.
* **Charge-excess emission** results from the net negative charge in the shower
  front. Its polarization is radial with respect to the shower axis.

The module provides two microscopic, time-domain formalisms:

* ``CoREAS`` implements the endpoint calculation used by the CoREAS extension
  of CORSIKA 7. Near the Cherenkov angle, where the endpoint expression becomes
  numerically unstable, it uses a ZHS-like limiting calculation.
* ``ZHS`` evaluates the radiation from straight particle-track segments using
  the ZHS algorithm.

Both formalisms operate on the same CORSIKA 8 tracks and can therefore be
instantiated together for a direct comparison without shower-to-shower
fluctuations. Their agreement is a useful internal validation of a simulation
configuration.

Configuration model
-------------------

A ``RadioProcess<TObservers, TRadioImpl, TPropagator>`` combines three
user-configured parts:

* An **observer collection** specifies the locations that receive the signal.
  ``TimeDomainObserver`` stores the sampled electric-field components in a
  configurable time window.
* The **formalism** is either ``CoREAS`` or ``ZHS``.
* A **propagator** supplies the signal path and propagation time from an
  emission point to an observer. This separation permits a common formalism to
  be used with different medium models and propagation strategies.

The ``c8_air_shower`` application contains a complete example. It creates
separate observer collections for ``CoREAS`` and ``ZHS``, builds a tabulated
flat-atmosphere propagator, and registers both radio processes with the output
manager. Register a radio process before starting the output lifecycle so that
its observer traces are written with the rest of the shower output.

.. code-block:: cpp

   auto propagator = make_tabulated_flat_atmosphere_radio_propagator(
       environment, injectionPosition, surface, step);

   RadioProcess<decltype(detector),
                CoREAS<decltype(detector), decltype(propagator)>,
                decltype(propagator)>
       coreas(detector, propagator);
   output.add("CoREAS", coreas);

The output is organized per radio process and includes observer traces. The
Python ``corsika8.io.outputs.RadioProcess`` reader loads the observer data from
``observers.parquet`` into time and ``Ex``, ``Ey``, and ``Ez`` columns.

Track-step accuracy
-------------------

Radio calculations resolve the phase relationship of radiation from consecutive
particle-track segments. Consequently, track length can affect the predicted
fluence and radiation energy, especially near Cherenkov-like time compression.
The radio-module validation found that reducing the maximum magnetic deflection
per leap-frog tracking step from the default :math:`0.2\,\mathrm{rad}` to
:math:`0.001\,\mathrm{rad}` brought the ``CoREAS`` and ``ZHS`` predictions for
identical air showers into close agreement.

For precision air-shower radio studies, configure curved leap-frog tracking
with ``maxDeflection = 0.001`` as a starting point, then check convergence for
the chosen energy, geometry, medium, observer layout, and frequency band. The
smaller value increases the number of particle-track segments and hence the
runtime. It is an accuracy setting for radio emission, not a replacement for a
configuration-specific convergence study.

Further reading
---------------

* J. M. Alameddine et al., *Simulating radio emission from particle cascades
  with CORSIKA 8*, Astroparticle Physics **166** (2025) 103072,
  `doi:10.1016/j.astropartphys.2024.103072
  <https://doi.org/10.1016/j.astropartphys.2024.103072>`_.
* C. W. James, H. Falcke, T. Huege, and M. Ludwig, *General description of
  electromagnetic radiation processes based on instantaneous charge
  acceleration in endpoints*, Physical Review E **84** (2011) 056602.
* J. Alvarez-Muniz, A. Romero-Wolf, and E. Zas, *Cherenkov radio pulses from
  electromagnetic showers in the time-domain*, Physical Review D **81** (2010)
  123009.