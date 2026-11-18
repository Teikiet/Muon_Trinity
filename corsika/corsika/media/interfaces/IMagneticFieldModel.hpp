/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/Vector.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

namespace corsika {
  namespace media {
    /**
     * An interface for magnetic field models.
     *
     * This is the base interface for magnetic field model mixins.
     *
     */
    template <typename Model>
    class IMagneticFieldModel : public Model {

      // a type-alias for a magnetic field vector
      using MagneticFieldVector = Vector<magnetic_flux_density_d>;

    public:
      /**
       * Evaluate the magnetic field at a given location.
       *
       * @param  point    The location to evaluate the field at.
       * @returns    The magnetic field vector at that point.
       */
      virtual auto getMagneticField(Point const&) const -> MagneticFieldVector = 0;

      /**
       * A virtual default destructor.
       */
      virtual ~IMagneticFieldModel() = default; // LCOV_EXCL_LINE

    }; // END: class MagneticField

  } // namespace media
} // namespace corsika
