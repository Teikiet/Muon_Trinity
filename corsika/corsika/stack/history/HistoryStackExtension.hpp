/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/Logging.hpp>
#include <corsika/framework/stack/Stack.hpp>

#include <memory>
#include <utility>
#include <vector>

namespace corsika::history {

  /**
   * @class HistoryData
   *
   * definition of stack-data object to store history information this
   * is vector with shared_ptr<TEvent>, where TEvent is a free
   * template parameter for customization.
   */
  template <typename TEvent>
  class HistoryData {
    using EventPtr =
        std::shared_ptr<TEvent>; //!< Pointer to the event where this particle was created
    using ParentEventIndex = int; //!< index to TEvent::secondaries_
    using DataType = std::pair<EventPtr, ParentEventIndex>;

  public:
    // these functions are needed for the Stack interface
    void clear() { historyData_.clear(); }
    unsigned int getSize() const { return historyData_.size(); }
    unsigned int getCapacity() const { return historyData_.size(); }
    void copy(const int i1, const int i2) { historyData_[i2] = historyData_[i1]; }
    void swap(const int i1, const int i2) {
      std::swap(historyData_[i1], historyData_[i2]);
    }

    // custom data access function
    void setEvent(const int i, EventPtr v) { historyData_[i].first = std::move(v); }
    const EventPtr& getEvent(const int i) const { return historyData_[i].first; }
    EventPtr& getEvent(const int i) { return historyData_[i].first; }

    void setParentEventIndex(const int i, ParentEventIndex v) {
      historyData_[i].second = std::move(v);
    }
    ParentEventIndex getParentEventIndex(const int i) const {
      return historyData_[i].second;
    }

    // these functions are also needed by the Stack interface
    void incrementSize() { historyData_.push_back(DataType{}); }
    void decrementSize() {
      if (historyData_.size() > 0) { historyData_.pop_back(); }
    }

    // custom private data section
  private:
    std::vector<DataType> historyData_;
  };

  /**
   * @class HistoryDataInterface
   *
   * corresponding definition of a stack-readout object, the iteractor
   * dereference operator will deliver access to these function.
   */
  template <typename T, typename TEvent>
  class HistoryDataInterface : public T {
  protected:
    using T::getStack;
    using T::getStackData;

  public:
    using T::getIndex;

  public:
    // create a new particle from scratch
    void setParticleData() {
      CORSIKA_LOG_TRACE("HistoyDatatInterface::setParticleData()");
      getStackData().setParentEventIndex(getIndex(), -1);
    }

    // create a new particle as secondary of a parent
    void setParticleData(HistoryDataInterface& /*parent*/) {
      CORSIKA_LOG_TRACE("HistoyDatatInterface::setParticleData(parnt)");
      setParticleData();
    }

    void setEvent(const std::shared_ptr<TEvent>& v) {
      getStackData().setEvent(getIndex(), v);
    }

    void setParentEventIndex(int index) {
      getStackData().setParentEventIndex(getIndex(), index);
    }

    std::shared_ptr<TEvent> getEvent() const {
      return getStackData().getEvent(getIndex());
    }

    int getParentEventIndex() const {
      return getStackData().getParentEventIndex(getIndex());
    }

    std::string asString() const {
      return fmt::format("i_parent={}, [evt: {}]", getParentEventIndex(),
                         (bool(getEvent()) ? getEvent()->as_string() : "n/a"));
    }
  };

  template <typename T, typename TEvent>
  struct MakeHistoryDataInterface {
    typedef HistoryDataInterface<T, TEvent> type;
  };

} // namespace corsika::history

// for user-friendlyness we create the HistoryDataInterface type
// with the histoy::Event data content right here:

#include <corsika/stack/history/Event.hpp>

namespace corsika::history {

  template <typename TStackIter>
  using HistoryEventDataInterface =
      typename history::MakeHistoryDataInterface<TStackIter, history::Event>::type;

  using HistoryEventData = history::HistoryData<history::Event>;

} // namespace corsika::history
