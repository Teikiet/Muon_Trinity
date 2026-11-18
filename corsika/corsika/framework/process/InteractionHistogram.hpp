/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/utility/SaveBoostHistogram.hpp>

#include <boost/histogram.hpp>

#include <fstream>
#include <functional>
#include <map>
#include <utility>

#include <corsika/detail/framework/process/InteractionHistogram.hpp> // for detail namespace

namespace corsika {

  /** @ingroup Processes
      @{

       Class that creates and stores histograms of collisions
       @f$dN/dE_{lab}@f$, @f$dN/d\sqrt{s}@f$ which is used by class
       InteractionCounter

       Histograms are of type boost::histogram

  */

  class InteractionHistogram {
    static double constexpr lower_edge_cms = 1e3, upper_edge_cms = 1e17; // eV sqrt s
    static double constexpr lower_edge_lab = 1e3, upper_edge_lab = 1e21; // eV lab
    static unsigned int constexpr num_bins_lab = 18 * 10, num_bins_cms = 14 * 10;

    /**
       hist_type is a boost::histogram with two axes
        - a growing PDG id axis
        - a fixed logarithmic energy axis as configured by the user
     */
    using hist_type =
        decltype(detail::hist_factory(num_bins_lab, lower_edge_lab, upper_edge_lab));

    hist_type inthist_cms_, inthist_lab_;

  public:
    InteractionHistogram();

    /**
       fill both CMS and lab histograms at the same time
       @param projectile_id corsika::Code of particle
       @param lab_energy Energy in lab. frame
       @param mass_target Mass of target particle
       @param A if projectile_id is Nucleus : Mass of nucleus
       @param Z if projectile_id is Nucleus : Charge of nucleus
    */
    void fill(Code const projectile_id, HEPEnergyType const lab_energy,
              HEPEnergyType const mass_target);

    //! return histogram in c.m.s. frame
    hist_type const& CMSHist() const { return inthist_cms_; }
    /// return histogram in laboratory frame
    hist_type const& labHist() const { return inthist_lab_; }

    /**
     * Clear all bins and restore axes to their initial state (as after construction).
     * Use at shower boundaries for per-shower histograms.
     */
    void reset();

    InteractionHistogram& operator+=(InteractionHistogram const& other);
    InteractionHistogram operator+(InteractionHistogram other) const;
  };

  /** @} */

} // namespace corsika

#include <corsika/detail/framework/process/InteractionHistogram.inl> // for implementation
