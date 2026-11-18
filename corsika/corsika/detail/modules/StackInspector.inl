/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/modules/StackInspector.hpp>

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>

namespace corsika {

  template <typename TStack>
  inline StackInspector<TStack>::StackInspector(int const vNStep, bool const vReportStack,
                                                HEPEnergyType const vE0)
      : StackProcess<StackInspector<TStack>>(vNStep)
      , ReportStack_(vReportStack)
      , E0_(vE0)
      , StartTime_(std::chrono::system_clock::now())
      , energyPostInit_(HEPEnergyType::zero())
      , timePostInit_(std::chrono::system_clock::now()) {}

  template <typename TStack>
  inline StackInspector<TStack>::~StackInspector() {}

  template <typename TStack>
  inline void StackInspector<TStack>::doStack(TStack const& vS) {

    [[maybe_unused]] int i = 0;
    HEPEnergyType Etot = 0_GeV;

    for (auto const& iterP : vS) {
      HEPEnergyType E = iterP.getEnergy();
      Etot += E;
      if (ReportStack_) {
        auto pos = iterP.getPosition().getCoordinates();
        CORSIKA_LOG_INFO(
            "StackInspector: i= {:5d}"
            ", id= {:^15}"
            " E={:15e} GeV, "
            " pos= {}"
            " node = {}",
            (i++), iterP.getPID(), (E / 1_GeV), pos, fmt::ptr(iterP.getNode()));
      }
    }

    if ((E0_ - Etot) < dE_threshold_) return;

    // limit number of printouts
    if (PrintoutCounter_ < MaxNumberOfPrintouts_) {
      std::chrono::system_clock::time_point const now = std::chrono::system_clock::now();
      std::chrono::duration<double> const elapsed_seconds = now - StartTime_; // seconds

      // Select reference times and energies using either the true start
      // or the delayed start to avoid counting overhead in ETA
      bool const usePostVals = (energyPostInit_ != HEPEnergyType::zero());
      auto const dE = (usePostVals ? energyPostInit_ : E0_) - Etot;
      std::chrono::duration<double> const usedSeconds = now - timePostInit_;
      double const dT = usedSeconds.count();

      double const progress = (E0_ - Etot) / E0_;
      // for printout
      std::time_t const now_time = std::chrono::system_clock::to_time_t(now);

      double const ETA_seconds = (1.0 - progress) * dT * (E0_ / dE);
      std::chrono::system_clock::time_point const ETA =
          now + std::chrono::seconds((int)ETA_seconds);

      std::time_t const ETA_time = std::chrono::system_clock::to_time_t(ETA);

      int const yday0 = std::localtime(&now_time)->tm_yday;
      int const yday1 = std::localtime(&ETA_time)->tm_yday;
      int const dyday = yday1 - yday0;

      std::stringstream ETA_string;
      ETA_string << std::put_time(std::localtime(&ETA_time), "%T %Z");

      CORSIKA_LOG_DEBUG(
          "StackInspector: "
          " time={}"
          ", running={:.1f} seconds"
          " ( {:.1f}%)"
          ", nStep={}"
          ", stackSize={}"
          ", Estack={:.1f} GeV"
          ", ETA={}{}",
          std::put_time(std::localtime(&now_time), "%T %Z"), elapsed_seconds.count(),
          (progress * 100), getStep(), vS.getSize(), Etot / 1_GeV,
          (dyday == 0 ? "" : fmt::format("+{}d ", dyday)), ETA_string.str());

      PrintoutCounter_++;

      if (PrintoutCounter_ == MaxNumberOfPrintouts_) {
        CORSIKA_LOG_DEBUG("StackInspector reached allowed maximum of {} lines printout",
                          MaxNumberOfPrintouts_);
      }

      // Change reference time once the shower has begin (avoid counting overhead time)
      if (progress > 0.02 && energyPostInit_ == HEPEnergyType::zero()) {
        energyPostInit_ = Etot;
        timePostInit_ = std::chrono::system_clock::now();
      }
    }
  }

} // namespace corsika
