/*
 * (c) Copyright 2021 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/process/ContinuousProcessIndex.hpp>

namespace corsika {

  /**
   * @ingroup Processes
   * To store step length in LengthType and unique index in ProcessSequence of shortest
   * step ContinuousProcess.
   */

  class ContinuousProcessStepLength {
  public:
    ContinuousProcessStepLength(LengthType const step)
        : step_(step)
        , id_(0) {}
    ContinuousProcessStepLength(LengthType const step, ContinuousProcessIndex const id)
        : step_(step)
        , id_(id) {}
    void setLength(LengthType const step) { step_ = step; }
    LengthType getLength() const { return step_; }
    void setIndex(ContinuousProcessIndex const id) { id_ = id; }
    ContinuousProcessIndex getIndex() const { return id_; }
    operator ContinuousProcessIndex() const { return id_; }
    operator LengthType() const { return step_; }
    bool operator<(ContinuousProcessStepLength const& r) const { return step_ < r.step_; }

  private:
    LengthType step_;
    ContinuousProcessIndex id_;
  };

} // namespace corsika
