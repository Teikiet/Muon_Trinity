/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/media/Environment.hpp>

namespace corsika {
  namespace media {

    template <typename IEnvironmentModel>
    inline media::Environment<IEnvironmentModel>::Environment()
        : coordinateSystem_{get_root_CoordinateSystem()}
        , universe_(std::make_unique<BaseNodeType>(
              std::make_unique<Universe>(coordinateSystem_))) {}

    template <typename IEnvironmentModel>
    inline typename media::Environment<IEnvironmentModel>::BaseNodeType::VTNUPtr&
    media::Environment<IEnvironmentModel>::getUniverse() {
      return universe_;
    }

    template <typename IEnvironmentModel>
    inline typename media::Environment<IEnvironmentModel>::BaseNodeType::VTNUPtr const&
    media::Environment<IEnvironmentModel>::getUniverse() const {
      return universe_;
    }

    template <typename IEnvironmentModel>
    inline const CoordinateSystem&
    media::Environment<IEnvironmentModel>::getCoordinateSystem() const {
      return coordinateSystem_;
    }

    // factory method for creation of VolumeTreeNodes
    template <typename IEnvironmentModel>
    template <class TVolumeType, typename... TVolumeArgs>
    std::unique_ptr<media::VolumeTreeNode<IEnvironmentModel> > inline media::Environment<
        IEnvironmentModel>::createNode(TVolumeArgs&&... args) {
      static_assert(std::is_base_of_v<IVolume, TVolumeType>,
                    "unusable type provided, needs to be derived from "
                    "\"Volume\"");

      return std::make_unique<BaseNodeType>(
          std::make_unique<TVolumeType>(std::forward<TVolumeArgs>(args)...));
    }

    template <typename IEnvironmentModel>
    std::set<Code> const get_all_elements_in_universe(
        media::Environment<IEnvironmentModel> const& env) {

      auto const& universe = *(env.getUniverse());
      auto const allElementsInUniverse = std::invoke([&]() {
        std::set<Code> allElementsInUniverse;
        auto collectElements = [&](auto& vtn) {
          if (vtn.hasModelProperties()) {
            auto const& comp =
                vtn.getModelProperties().getNuclearComposition().getComponents();
            for (auto const c : comp) allElementsInUniverse.insert(c);
          }
        };
        universe.walk(collectElements);
        return allElementsInUniverse;
      });
      return allElementsInUniverse;
    }

  } // namespace media
} // namespace corsika
