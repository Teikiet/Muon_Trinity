/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/stack/SecondaryView.hpp>
#include <corsika/stack/history/Event.hpp>

#include <corsika/framework/core/Logging.hpp>

#include <boost/type_index.hpp>

#include <memory>
#include <type_traits>
#include <utility>

namespace corsika::history {

  //! mix-in class for SecondaryView that fills secondaries into an \class Event
  template <class T1, template <class> class T2>
  class HistorySecondaryProducer {
  public:
    EventPtr event_;
    static bool constexpr has_event{true};

  public:
    template <typename Particle>
    HistorySecondaryProducer(Particle const& p)
        : event_{std::make_shared<Event>()} {
      CORSIKA_LOG_TRACE("HistorySecondaryProducer::HistorySecondaryProducer");
      event_->setProjectileIndex(p.getIndex());
      event_->setParentEvent(p.getEvent());
    }

    /**
     * Method is called after a new Secondary has been created on the
     * SecondaryView. Extra logic can be introduced here.
     *
     * The input Particle is the new secondary that was produced and
     * is of course a reference into the SecondaryView itself.
     */
    template <typename Particle>
    auto new_secondary(Particle& sec) {
      CORSIKA_LOG_TRACE("HistorySecondaryProducer::new_secondary(sec)");

      // store particles at production time in Event here
      auto const sec_index =
          event_->addSecondary(sec.getEnergy(), sec.getMomentum(), sec.getPID());
      sec.setParentEventIndex(sec_index);
      sec.setEvent(event_);
    }
  };

} // namespace corsika::history
