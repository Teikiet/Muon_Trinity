/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/Logging.hpp>
#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/stack/history/EventType.hpp>
#include <corsika/stack/history/SecondaryParticle.hpp>

#include <memory>
#include <optional>
#include <vector>

namespace corsika::history {

  class Event;
  using EventPtr = std::shared_ptr<history::Event>;

  class Event {

    size_t projectile_index_ = 0; //!< index of projectile on stack
    std::vector<SecondaryParticle> secondaries_;
    EventPtr parent_event_;

    EventType type_ = EventType::Uninitialized;

    std::optional<Code> targetCode_;

  public:
    Event() = default;

    void setParentEvent(EventPtr const& evt) { parent_event_ = evt; }

    bool hasParentEvent() const { return bool(parent_event_); }
    EventPtr& parentEvent() { return parent_event_; }
    EventPtr const& parentEvent() const { return parent_event_; }

    void setProjectileIndex(size_t i) { projectile_index_ = i; }
    size_t projectileIndex() const { return projectile_index_; }

    size_t addSecondary(HEPEnergyType energy, Vector<hepmomentum_d> const& momentum,
                        Code pid) {
      secondaries_.emplace_back(energy, momentum, pid);
      return secondaries_.size() - 1;
    }

    std::vector<SecondaryParticle> const& secondaries() const { return secondaries_; }

    void setTargetCode(const Code t) { targetCode_ = t; }

    std::string as_string() const {
      return fmt::format("hasParent={}, projIndex={}, Nsec={}", hasParentEvent(),
                         projectile_index_, secondaries_.size());
    }

    EventType eventType() const { return type_; }

    void setEventType(EventType t) { type_ = t; }
  };

} // namespace corsika::history
