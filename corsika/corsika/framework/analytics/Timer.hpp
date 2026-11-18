/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <chrono>
#include <utility>

#include <type_traits>

namespace corsika {

  template <typename TClock = std::chrono::high_resolution_clock,
            typename TDuration = std::chrono::microseconds>
  class Timer {
  protected:
    /// Default clock used for time measurement
    using clock_type = TClock;

    /// Internal resolution of the time measurement
    using duration_type = TDuration;

    /// Startpoint of time measurement
    typename clock_type::time_point start_;

    /// Measured runtime of the function
    duration_type timeDiff_;

  public:
    Timer()
        : timeDiff_(0){};

    /// Start the timer
    void startTimer() { start_ = clock_type::now(); }

    /// Stop the timer
    void stopTimer() {
      timeDiff_ = std::chrono::duration_cast<duration_type>(clock_type::now() - start_);
    }

    /**
     * Returns the runtime of the last call to the wrapped function
     * @return Returns the measured runtime of the wrapped function/functor in the unit
     * given by TDuration
     **/
    duration_type getTime() const { return timeDiff_; }
  };

  std::false_type is_timer_impl(...);
  template <typename T, typename U>
  std::true_type is_timer_impl(Timer<T, U> const volatile&);

  template <typename T>
  constexpr bool is_timer_v =
      std::is_same_v<decltype(is_timer_impl(std::declval<T&>())), std::true_type>;

} // namespace corsika
