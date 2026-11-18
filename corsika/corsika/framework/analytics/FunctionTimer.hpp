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

#include <corsika/framework/analytics/Timer.hpp>

namespace corsika {

  /// Wraps and measures the runtime of a single function type object
  /**
   *
   * @tparam TFunc funtion pointer that should be wrapped
   * @tparam TClock type of the clock that should be used for measurements, default is
   * high_resolution_clock
   * @tparam TDuration type of std::duration to measure the elapsed time, default is
   * microseconds
   */
  template <typename TFunc, class TTimer = Timer<std::chrono::high_resolution_clock,
                                                 std::chrono::microseconds>>
  class FunctionTimer : public TTimer {
    static_assert(
        is_timer_v<TTimer>,
        "TTimer is not a timer!"); // Better
                                   // https://en.cppreference.com/w/cpp/language/constraints
                                   // but not available in C++17
  public:
    /** Constructs the wrapper with the given functionpointer
     *  @param f  Function or functor whose runtime should be measured
     **/
    FunctionTimer(TFunc f);

    /** Functor for calling the wrapped function
     *  This functor calls the wrapped function and measures the elapsed time between call
     *and return. The return value needs to be copy constructible.
     *  @tparam TArgs Parameter types that are forwarded to the function. The use of
     *correct types is the responsibility of the user, no checks are done.
     *  @param args Arguments are forwarded to the wrapped function. This method does not
     *support overloaded function resolution.
     *  @return The return value of the wrapped function is temporarily copied and then
     *returned by value
     **/
    template <typename... TArgs>
    auto operator()(TArgs&&... args) -> std::invoke_result_t<TFunc, TArgs...>;

  private:
    TFunc function_;
  };

} // namespace corsika

#include <corsika/detail/framework/analytics/FunctionTimer.inl>