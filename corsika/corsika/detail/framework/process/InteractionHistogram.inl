/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <fstream>
#include <string>

#include <corsika/framework/process/InteractionHistogram.hpp>
#include <corsika/detail/framework/process/InteractionHistogram.hpp> // for detail namespace

namespace corsika {

  inline InteractionHistogram::InteractionHistogram()
      : inthist_cms_{detail::hist_factory(num_bins_cms, lower_edge_cms, upper_edge_cms)}
      , inthist_lab_{detail::hist_factory(num_bins_lab, lower_edge_lab, upper_edge_lab)} {
  }

  inline void InteractionHistogram::reset() {
    inthist_cms_ = detail::hist_factory(num_bins_cms, lower_edge_cms, upper_edge_cms);
    inthist_lab_ = detail::hist_factory(num_bins_lab, lower_edge_lab, upper_edge_lab);
  }

  inline void InteractionHistogram::fill(Code const projectile_id,
                                         HEPEnergyType const lab_energy,
                                         HEPEnergyType const mass_target) {
    auto constexpr inv_eV = 1 / 1_eV;
    auto const projectile_mass = get_mass(projectile_id);
    auto const sqrtS = sqrt(projectile_mass * projectile_mass +
                            mass_target * mass_target + 2 * lab_energy * mass_target);

    CORSIKA_LOG_DEBUG("pM={}, tM={}, pid={}, Elab={}, sqrtS={}, pdg={} a={} z={}",
                      projectile_mass / 1_GeV, mass_target / 1_GeV, projectile_id,
                      lab_energy / 1_GeV, sqrtS / 1_GeV, get_PDG(projectile_id),
                      get_nucleus_A(projectile_id), get_nucleus_Z(projectile_id));

    inthist_cms_(static_cast<int>(get_PDG(projectile_id)), sqrtS * inv_eV);
    inthist_lab_(static_cast<int>(get_PDG(projectile_id)), lab_energy * inv_eV);
  }

  inline InteractionHistogram& InteractionHistogram::operator+=(
      InteractionHistogram const& other) {
    inthist_lab_ += other.inthist_lab_;
    inthist_cms_ += other.inthist_cms_;
    return *this;
  }

  inline InteractionHistogram InteractionHistogram::operator+(
      InteractionHistogram other) const {
    other.inthist_lab_ += inthist_lab_;
    other.inthist_cms_ += inthist_cms_;
    return other;
  }
} // namespace corsika
