/*
 * (c) Copyright 2021 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

namespace corsika {

  /**
   * @ingroup Processes
   * To index individual processes (continuous processes) inside a
   * ProcessSequence.
   */

  class ContinuousProcessIndex {
  public:
    ContinuousProcessIndex()
        : id_(nullptr) {} // default
    ContinuousProcessIndex(void const* id)
        : id_(id) {}
    void setIndex(void const* id) { id_ = id; }
    void const* getIndex() const { return id_; }
    bool operator==(ContinuousProcessIndex const v) const { return id_ == v.id_; }
    bool operator!=(ContinuousProcessIndex const v) const { return !(*this == v); }

  private:
    void const* id_;
  };

} // namespace corsika
