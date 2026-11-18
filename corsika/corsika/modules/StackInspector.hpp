/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/process/StackProcess.hpp>

#include <chrono>

namespace corsika {

  /**
   * StackProcess that will act each @f$n_{step}@f$ steps to perform diagnostics on the
   * full stack.
   *
   * The StackInspector can dump the entrie stack content for debugging, or also just
   * determine the total energy remaining on the stack. From the decrease of energy on the
   * stack an ETA for the completion of the simulation is determined.
   *
   * @tparam TStack Is the type of the particle stack.
   */

  template <typename TStack>
  class StackInspector : public StackProcess<StackInspector<TStack>> {

    typedef typename TStack::particle_type Particle;

    using StackProcess<StackInspector<TStack>>::getStep;

  public:
    StackInspector(int const nStep, bool const reportStack, HEPEnergyType const vE0);
    ~StackInspector();

    void doStack(TStack const&);

    /**
     * To set a new E0, for example when a new shower event is started.
     */
    void setE0(HEPEnergyType const E0) { E0_ = E0; }

  private:
    bool ReportStack_;
    int PrintoutCounter_ = 0;
    const int MaxNumberOfPrintouts_ = 10;
    HEPEnergyType E0_;
    const HEPEnergyType dE_threshold_ = 1_eV;
    std::chrono::system_clock::time_point StartTime_;
    HEPEnergyType energyPostInit_;
    std::chrono::system_clock::time_point timePostInit_;
  };

} // namespace corsika

#include <corsika/detail/modules/StackInspector.inl>
