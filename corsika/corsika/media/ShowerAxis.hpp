/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/media/Environment.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/Vector.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

#include <cstdlib>
#include <fstream>
#include <functional>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <vector>

#include <iostream>

#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace corsika {
  namespace media {
    /**
     * \class ShowerAxis
     *
     * The environment::media::ShowerAxis is created from a Point and
     * a Vector and inside an Environment. It internally uses
     * a table with steps=10000 (default) rows for interpolation.
     *
     * The shower axis can convert location in the shower into a
     * projected grammage along the shower axis.
     *
     **/

    ///\todo documentation needs update ...
    class ShowerAxis {
    public:
      template <typename TEnvModel>
      ShowerAxis(Point const& pStart, Point const& pEnd,
                 media::Environment<TEnvModel> const& env, bool const doThrow = false,
                 int const steps = 10'000);

      template <typename TEnvModel>
      ShowerAxis(Point const& pStart, Vector<length_d> const& length,
                 media::Environment<TEnvModel> const& env, bool const doThrow = false,
                 int const steps = 10'000);

      LengthType getSteplength() const;

      GrammageType getMaximumX() const;

      GrammageType getMinimumX() const;

      /**
       * Returns the grammage along the shower axis of the projection of a point p
       * onto the shower axis.
       * Will return either getMinimumX() or getMaximumX() in case the projection is
       * outside the shower axis.
       *
       *  @param p Point to project onto the shower axis.
       *  @retval Grammage along shower axis for projection of point p.
       */
      GrammageType getProjectedX(Point const& p) const;

      /**
       * Returns the grammage along the shower axis for a given length along the shower
       * axis. Will return either getMinimumX() or getMaximumX() in case the length is
       * outside the shower axis.
       *
       *  @param l Length along shower axis.
       *  @retval Grammage along shower axis for length l.
       */
      GrammageType getX(LengthType l) const;

      DirectionVector const& getDirection() const;

      Point const& getStart() const;

    private:
      Point const pointStart_;
      Vector<length_d> const length_;
      bool throw_ = false;
      LengthType const max_length_, steplength_;
      DirectionVector const axis_normalized_;
      std::vector<GrammageType> X_;
    };

  } // namespace media
} // namespace corsika

#include <corsika/detail/media/ShowerAxis.inl>
