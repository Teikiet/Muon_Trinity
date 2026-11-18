/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/Logging.hpp>
#include <corsika/framework/stack/Stack.hpp>

#include <tuple>
#include <utility>
#include <vector>

namespace corsika::weights {

  /**
   * Describe "particle weights" on a Stack.
   *
   * Corresponding defintion of a stack-readout object, the iteractor
   * dereference operator will deliver access to these function
   * defintion of a stack-readout object, the iteractor dereference
   * operator will deliver access to these function.
   */

  /**
   * \tparam TParentStack The stack to be exteneded here with weight data
   */
  template <typename TParentStack>
  struct WeightDataInterface : public TParentStack {

    typedef TParentStack super_type;

  public:
    // default version for particle-creation from input data
    void setParticleData(std::tuple<double> const v);
    void setParticleData(WeightDataInterface const& parent, std::tuple<double> const);
    void setParticleData();
    void setParticleData(WeightDataInterface const& parent);

    std::string asString() const;

    void setWeight(double const v);

    double getWeight() const;
  };

  // definition of stack-data object to store geometry information

  /**
   * @class WeightData
   *
   * definition of stack-data object to store geometry information.
   */
  class WeightData {

  public:
    typedef std::vector<double> weight_vector_type;

    WeightData() = default;
    WeightData(WeightData const&) = default;
    WeightData(WeightData&&) = default;
    WeightData& operator=(WeightData const&) = default;
    WeightData& operator=(WeightData&&) = default;

    // these functions are needed for the Stack interface
    void clear();

    unsigned int getSize() const;

    unsigned int getCapacity() const;

    void copy(int const i1, int const i2);

    void swap(int const i1, int const i2);

    // custom data access function
    void setWeight(int const i, double const v);

    double getWeight(int const i) const;

    // these functions are also needed by the Stack interface
    void incrementSize();

    void decrementSize();

    // custom private data section
  private:
    weight_vector_type weight_vector_;
  };

  template <typename TParentStack>
  struct MakeWeightDataInterface {
    typedef WeightDataInterface<TParentStack> type;
  };

} // namespace corsika::weights

#include <corsika/detail/stack/WeightStackExtension.inl>
