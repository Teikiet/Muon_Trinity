/*
 * (c) Copyright 2021 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/stack/SecondaryView.hpp>

namespace corsika {

  /**
   * Class to handle the generation of new secondaries. Used as default mix-in for
   * SecondaryView.
   */
  template <class T1, template <class> class T2>
  class DefaultSecondaryProducer {

    using view_type = SecondaryView<T1, T2, DefaultSecondaryProducer>;

  public:
    static bool constexpr has_event{false};

    /**
     * Method is called after a new secondary has been created on the
     * SecondaryView. Extra logic can be introduced here.
     *
     * The input Particle is the new secondary that was produced and
     * is of course a reference into the SecondaryView itself.
     */
    template <typename Particle>
    void new_secondary(Particle&&) const {
      CORSIKA_LOG_TRACE("DefaultSecondaryProducer::new_secondary(Particle&&)");
    }

    /**
     * Method is called when a new SecondaryView is being created
     * created. Extra logic can be introduced here.
     *
     * The input Particle is a reference object into the original
     * parent stack! It is not a reference into the SecondaryView
     * itself.
     */
    template <typename Particle>
    DefaultSecondaryProducer(Particle const&) {

      CORSIKA_LOG_TRACE("DefaultSecondaryProducer::DefaultSecondaryProducer(Particle&)");
    }
  };

} // namespace corsika