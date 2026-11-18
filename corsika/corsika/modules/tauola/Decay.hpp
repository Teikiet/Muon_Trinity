/*
 * (c) Copyright 2023 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/modules/tauola/TauolaInterfaceParticle.hpp>
#include <corsika/framework/process/DecayProcess.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/modules/Random.hpp>

namespace corsika::tauola {

  /**
   * The various tau polarizations that we support.
   */
  enum class Helicity { LeftHanded, Unpolarized, RightHanded };

  auto inline format_as(Helicity code) { return fmt::underlying(code); }

  /**
   * An interface to TAUOLA to simulate tau lepton decays.
   *
   * This class interfaces CORSIKA8 with the TAUOLA decay generator
   * for fast tau- and tau+ decays. The generator can be configured
   * by various arguments to the constructor, namely:
   *
   *     helicity:  Since CORSIKA8 does not track polarization, this
   *                allows you to explicitly set the helicity of all
   *                taus decayed using TAUOLA. `helicity` is multiplied
   *                with the normalized tau momentum direction to
   *                calculate the polarization vector. The default setting is Unpolarized
   * tau radiative: This turns on/off radiative QED corrections for leptonic tau decays.
   * This defaults to off. spinncorr: This turns the spin-correlation effects on/off.
   *
   * `radiative` and `spincorr` can only be set once during the
   * first initialization of TAUOLA. Additional attempts to change
   * these parameters by constructing a new TAUOLA() instance will have
   * no effect since the TAUOLA kernel has already been initialized.
   *
   */
  class Decay : public DecayProcess<Decay> {

    const Helicity helicity_; ///< The helicity that is used for tau decays (+1, 0, -1)

    /**
     * Decay a single TauolaInterfaceParticle.
     *
     * This handles creating the polarization vector
     * and performing the decay of the particle. The particle
     * is modified in place and the daughters are added
     * to its particle recorded.
     *
     * @param particle    A TAUOLA particle instance.
     */
    template <typename InterfaceParticle>
    void decay(InterfaceParticle*);

  public:
    /**
     * Create a new TAUOLA Process.
     *
     * @param helicity     The helicity (LeftHanded, Unpolarized, RightHanded) to use for
     * decays
     * @param radiative    Whether to enable QCD radiative corrections.
     * @param spincorr     Whether to use spin-correlations.
     */
    Decay(const Helicity helicity = Helicity::Unpolarized, const bool radiative = false,
          const bool spincorr = true);

    ~Decay();

    /**
     * Check if TAUOLA should handle this decay.
     *
     * @param code    The PDG ID of the particle.
     *
     * @returns True if TAUOLA can handle this decay.
     */
    bool isDecayHandled(const Code);

    /**
     * Decay a tau lepton using TAUOLA.
     *
     * @param particle    The particle to decay.
     */
    template <typename TSecondaryView>
    void doDecay(TSecondaryView&);

    /**
     * Get the lifetime of a particle.
     *
     * @param particle    The particle to calculate the lifetime for.
     *
     */
    template <typename TParticle>
    TimeType getLifetime(TParticle const&);

    /**
     * Print a summary of the decays calculated with TAUOLA.
     */
    void PrintSummary() const;

    std::shared_ptr<spdlog::logger> logger_ = get_logger("corsika_tauola_decay");

    int count_ = 0; ///< The number of taus decayed with TAUOLA.

  private:
    /**
     * Set up RNG members to interface cleanly with TAUOLA
     */
    static corsika::rng_function_type rngFcn_;
    static double WrappedRNG();
  };

  // Define static member outside the class
  corsika::rng_function_type Decay::rngFcn_ = nullptr;

} // namespace corsika::tauola

#include <corsika/detail/modules/tauola/Decay.inl>
