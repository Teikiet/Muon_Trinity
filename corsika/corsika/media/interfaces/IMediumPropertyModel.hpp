/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/media/medium/MediumProperties.hpp>

#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

namespace corsika {
  namespace media {

    /**
     * An interface for type of media, needed e.g. to determine energy losses.
     *
     * This is the base interface for media types.
     */
    template <typename TModel>
    class IMediumPropertyModel : public TModel {

    public:
      /**
       * Evaluate the medium type at a given location.
       *
       * @param  point    The location to evaluate at.
       * @returns    The media type
       */
      virtual Medium getMedium() const = 0;

      /**
       * A virtual default destructor.
       */
      virtual ~IMediumPropertyModel() = default;

    }; // END: class IMediumTypeModel

  } // namespace media
} // namespace corsika
