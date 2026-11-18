/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/Vector.hpp>
#include <corsika/framework/process/ContinuousProcess.hpp>
#include <corsika/framework/core/Step.hpp>

#include <corsika/modules/writers/WriterOff.hpp>

#include <map>

namespace corsika {

  /**
   *   PDG2018, passage of particles through matter
   *
   * Note, that \f$I_{\mathrm{eff}}\f$ of composite media a determined from \f$ \ln I =
   * \sum_i a_i \ln(I_i) \f$ where \f$ a_i \f$ is the fraction of the electron population
   * (\f$\sim Z_i\f$) of the \f$i\f$-th element. This can also be used for shell
   * corrections or density effects.
   *
   * The \f$I_{\mathrm{eff}}\f$ of compounds is not better than a few percent, if not
   * measured explicitly.
   *
   * For shell correction, see Sec 6 of https://www.nap.edu/read/20066/chapter/8#115
   *
   */

  template <typename TOutput = WriterOff>
  class BetheBlochPDG : public ContinuousProcess<BetheBlochPDG<TOutput>>, public TOutput {

    using MeVgcm2 = decltype(1e6 * electronvolt / gram * square(1e-2 * meter));

  public:
    template <typename... TOutputArgs>
    BetheBlochPDG(TOutputArgs&&... args);

    /**
     * Interface function of ContinuousProcess.
     *
     * @param particle The particle to process in its current state
     * @param track The trajectory in space of this particle, on which doContinuous
     *should act
     * @param limitFlag flag to identify, if BetheBlochPDG::getMaxStepLength is the
     *        globally limiting factor (or not)
     clang-format-on **/
    template <typename TParticle>
    ProcessReturn doContinuous(Step<TParticle>& step, bool const limitFlag);

    template <typename TParticle, typename TTrajectory>
    LengthType getMaxStepLength(TParticle const&, TTrajectory const&) const;

    template <typename TParticle>
    static HEPEnergyType getBetheBloch(TParticle const&, const GrammageType);

    template <typename TParticle>
    static HEPEnergyType getRadiationLosses(TParticle const&, const GrammageType);

    template <typename TParticle>
    static HEPEnergyType getTotalEnergyLoss(TParticle const&, const GrammageType);

    YAML::Node getConfig() const override;
  };

} // namespace corsika

#include <corsika/detail/modules/energy_loss/BetheBlochPDG.inl>
