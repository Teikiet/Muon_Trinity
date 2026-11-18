/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/random/RNGManager.hpp>
#include <corsika/framework/utility/CorsikaData.hpp>
#include <corsika/framework/geometry/FourVector.hpp>

#include <boost/filesystem/path.hpp>
#include <boost/multi_array.hpp>

#include <array>
#include <utility>
#include <string>

namespace corsika::urqmd {

  class UrQMD {
  public:
    /**
     * The UrQMD interaction model.
     *
     * @param path Location of UrQMD XS data file.
     * @param retryFlag Internal UrQMD flag for retrying interaction in case of empty
     * event, 0 means retry.
     */
    UrQMD(boost::filesystem::path const path = corsika_data("UrQMD/UrQMD-1.3.1-xs.dat"),
          int const retryFlag = 0);

    bool isValid(Code const projectileId, Code const targetId) const;

    CrossSectionType getTabulatedCrossSection(Code const, Code const,
                                              HEPEnergyType const) const;

    CrossSectionType getCrossSection(Code const projectileId, Code const targetId,
                                     FourMomentum const& projP4,
                                     FourMomentum const& targP4) const;

    template <typename TView>
    void doInteraction(TView&, Code const projectile, Code const targetId,
                       FourMomentum const& projP4, FourMomentum const& targP4);

  private:
    void readXSFile(boost::filesystem::path);

    // data members
    default_prng_type& RNG_ = RNGManager<>::getInstance().getRandomStream("urqmd");
    std::uniform_int_distribution<int> booleanDist_{0, 1};
    int iflb_; //! // flag for retrying interaction in case of empty event, 0 means retry
    boost::multi_array<CrossSectionType, 3> xs_interp_support_table_;
  };

} // namespace corsika::urqmd

#include <corsika/detail/modules/urqmd/UrQMD.inl>
