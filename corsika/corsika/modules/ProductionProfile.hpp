/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/media/ShowerAxis.hpp>
#include <corsika/framework/process/SecondariesProcess.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/modules/writers/WriterOff.hpp>

#include <array>
#include <fstream>
#include <limits>
#include <string>

namespace corsika {

  /**
   * \class ProductionProfile
   *
   * \todo test missing
   *
   * is a SecondariesProcess, which is constructed from an environment::ShowerAxis
   * object, and a dX in units of g/cm2  (GrammageType).
   *
   */

  template <typename TOutput = WriterOff>
  class ProductionProfile : public SecondariesProcess<ProductionProfile<TOutput>>,
                            public TOutput {

  public:
    template <typename... TArgs>
    ProductionProfile(TArgs&&... args);

    /**
     * Count particles which are secondaries from discrete processes.
     *
     * @tparam TStackView
     */
    template <typename TStackView>
    void doSecondaries(TStackView&);

    long getCount() { return count_; };

    YAML::Node getConfig() const;
    long count_ = 0;
  };

} // namespace corsika

#include <corsika/detail/modules/ProductionProfile.inl>
