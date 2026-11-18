/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/media/interfaces/IMediumPropertyModel.hpp>

namespace corsika {
  namespace media {

    template <typename T>
    template <typename... Args>
    inline MediumPropertyModel<T>::MediumPropertyModel(Medium const medium,
                                                       Args&&... args)
        : T(std::forward<Args>(args)...)
        , medium_(medium) {}

    template <typename T>
    inline Medium MediumPropertyModel<T>::getMedium() const {
      return medium_;
    }

    template <typename T>
    inline void MediumPropertyModel<T>::setMedium(Medium const medium) {
      medium_ = medium;
    }

  } // namespace media
} // namespace corsika
