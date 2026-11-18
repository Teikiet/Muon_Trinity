/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/geometry/FourVector.hpp>
#include <corsika/modules/OnShellCheck.hpp>
#include <corsika/framework/core/Logging.hpp>

namespace corsika {

  inline OnShellCheck::OnShellCheck(double const vMassTolerance,
                                    double const vEnergyTolerance, bool const vError)
      : mass_tolerance_(vMassTolerance)
      , energy_tolerance_(vEnergyTolerance)
      , throw_error_(vError) {
    CORSIKA_LOGGER_DEBUG(logger_, "mass tolerance is set to {:3.2f}%",
                         mass_tolerance_ * 100);
    CORSIKA_LOGGER_DEBUG(logger_, "energy tolerance is set to {:3.2f}%",
                         energy_tolerance_ * 100);
  }

  inline OnShellCheck::~OnShellCheck() {
    logger_->info(
        " summary \n"
        " particles shifted: {} \n"
        " average energy shift (%): {} \n"
        " max. energy shift (%): {} ",
        int(count_), (count_ ? average_shift_ / count_ * 100 : 0), max_shift_ * 100.);
  }

  template <typename TView>
  inline void OnShellCheck::doSecondaries(TView& vS) {
    for (auto& p : vS) {
      auto const pid = p.getPID();
      if (is_nucleus(pid)) continue;
      auto const e_original = p.getEnergy();
      auto const p_original = p.getMomentum();
      auto const Plab = FourVector(e_original, p_original);
      auto const m_kinetic = Plab.getNorm();
      auto const m_corsika = get_mass(pid);
      auto const m_err_abs = abs(m_kinetic - m_corsika);
      if (m_err_abs >= mass_tolerance_ * m_corsika) {
        const HEPEnergyType e_shifted =
            sqrt(p_original.getSquaredNorm() + m_corsika * m_corsika);
        auto const e_shift_relative = (e_shifted / e_original - 1);
        count_ = count_ + 1;
        average_shift_ += abs(e_shift_relative);
        if (abs(e_shift_relative) > max_shift_) max_shift_ = abs(e_shift_relative);
        CORSIKA_LOGGER_TRACE(
            logger_,
            "shift particle mass for {} \n"
            "{:>45} {:7.5f} \n"
            "{:>45} {:7.5f} \n"
            "{:>45} {:7.5f} \n"
            "{:>45} {:7.5f} \n",
            pid, "corsika mass (GeV):", m_corsika / 1_GeV,
            "kinetic mass (GeV): ", m_kinetic / 1_GeV,
            "m_kin-m_cor (GeV): ", m_err_abs / 1_GeV,
            "mass tolerance (GeV): ", (m_corsika * mass_tolerance_) / 1_GeV);
        /*
          For now we warn if the necessary shift is larger than 1%.
          we could promote this to an error.
        */
        if (abs(e_shift_relative) > energy_tolerance_) {
          logger_->warn("warning! shifted particle energy by {} %",
                        e_shift_relative * 100);
          if (throw_error_) {
            throw std::runtime_error(
                "OnShellCheck: error! shifted energy by large amount!");
          }
        }

        // reset energy
        p.setEnergy(e_shifted);
      } else {
        CORSIKA_LOGGER_DEBUG(logger_, "particle mass for {} OK", pid);
      }
    }
  } // namespace corsika

} // namespace corsika
