/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/stack/CombinedStack.hpp>
#include <corsika/stack/GeometryNodeStackExtension.hpp>
#include <corsika/stack/VectorStack.hpp>
#include <corsika/stack/WeightStackExtension.hpp>
#include <corsika/stack/history/HistorySecondaryProducer.hpp>
#include <corsika/stack/history/HistoryStackExtension.hpp>
#include <corsika/media/Environment.hpp>
#include <corsika/media/interfaces/IMagneticFieldModel.hpp>
#include <corsika/media/interfaces/IMediumModel.hpp>
#include <corsika/media/interfaces/IMediumPropertyModel.hpp>

namespace corsika {

  // maybe use a similar copy of this file with defined templates for tests?
  using DummyEnvironmentInterface =
      media::IMediumPropertyModel<media::IMagneticFieldModel<media::IMediumModel>>;
  using DummyEnvironment = media::Environment<DummyEnvironmentInterface>;

  namespace test::detail {

    // ------------------------------------------
    // add geometry node tracking data to stack:

    // the GeometryNode stack needs to know the type of geometry-nodes from the
    // environment:
    template <typename TStackIter>
    using SetupGeometryDataInterface =
        typename node::MakeGeometryDataInterface<TStackIter, DummyEnvironment>::type;

    // combine particle data stack with geometry information for tracking
    template <typename TStackIter>
    using StackWithGeometryInterface =
        CombinedParticleInterface<VectorStack::pi_type, SetupGeometryDataInterface,
                                  TStackIter>;

    using StackWithGeometry =
        CombinedStack<typename VectorStack::stack_data_type,
                      node::GeometryData<DummyEnvironment>, StackWithGeometryInterface,
                      DefaultSecondaryProducer>;

    template <typename TStackIter>
    using SetupWeightDataInterface =
        typename weights::MakeWeightDataInterface<TStackIter>::type;

    template <typename TStackIter>
    using StackWithWeightInterface =
        CombinedParticleInterface<StackWithGeometry::pi_type, SetupWeightDataInterface,
                                  TStackIter>;

    using StackWithWeight =
        CombinedStack<typename StackWithGeometry::stack_data_type, weights::WeightData,
                      StackWithWeightInterface, DefaultSecondaryProducer>;

    // ------------------------------------------
    // Add [optional] history data to stack, too:

    // combine dummy stack with geometry information for tracking
    template <typename TStackIter>
    using StackWithHistoryInterface =
        CombinedParticleInterface<StackWithWeight::pi_type,
                                  history::HistoryEventDataInterface, TStackIter>;

    using StackWithHistory =
        CombinedStack<typename StackWithWeight::stack_data_type,
                      history::HistoryEventData, StackWithHistoryInterface,
                      history::HistorySecondaryProducer>;

  } // namespace test::detail

} // namespace corsika
