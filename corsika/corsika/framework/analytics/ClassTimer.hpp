/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

// Another possibility:
// https://en.wikibooks.org/wiki/More_C%2B%2B_Idioms/Execute-Around_Pointer
// for a more global approach
//
// In this case here only a single function is measured via member function pointer.

#pragma once

#include <chrono>
#include <utility>

#include <corsika/framework/analytics/Timer.hpp>

namespace corsika {

  template <class TClass, typename TTimer>
  class ClassTimerImpl : public TTimer {
    static_assert(
        is_timer_v<TTimer>,
        "TTimer is not a timer!"); // Better
                                   // https://en.cppreference.com/w/cpp/language/constraints
                                   // but not available in C++17

  protected:
    /// Reference to the class object on which the function should be called
    TClass& obj_;

  public:
    ClassTimerImpl(TClass& obj)
        : obj_(obj){};
  };

  /// Measure the runtime of a single class function
  /**
   * @tparam TClassFunc Type of the member function pointer that should be wrapped
   * @tparam TFunc Actual function of the type defined in TClass
   */
  template <typename TClassFunc, TClassFunc TFunc,
            typename TTimer =
                Timer<std::chrono::high_resolution_clock, std::chrono::microseconds>>
  class ClassTimer;

  /// Measure the runtime of a single class function
  /** Specialisation to capture exact information about the composition of the member
   * function pointer used.
   *
   *  This class wrapes a single function and allowes the measureing of its runtime if it
   * called via the "call(...)" function
   *
   * @tparam TClass Class of the function that should be wrapped
   * @tparam TRet   Return value of the wrapped function
   * @tparam TArgs  Arguments passed to the wrapped function
   * @tparam TFuncPtr Actual function of the type defined by TRet
   * TClass::TFuncPtr(TArgs...)
   */
  template <typename TClass, typename TRet, typename... TArgs,
            TRet (TClass::*TFuncPtr)(TArgs...), typename TTimer>
  class ClassTimer<TRet (TClass::*)(TArgs...), TFuncPtr, TTimer>
      : public ClassTimerImpl<TClass, TTimer> {
  private:
  public:
    ClassTimer(TClass& obj);

    /// Executes the wrapped function
    /** This function executes and measure the runtime of the wrapped function with the
     * highest precision available (high_resolution_clock).
     *
     * @param args Arguments are perfect forwarded to the wrapped function.
     * @return Returns the return value of the wrapped function. This value get copied
     * during the process and therefore must be copy constructible!
     */
    TRet call(TArgs... args);
  };

  /// Specialisation for member functions without return value
  template <typename TClass, typename... TArgs, void (TClass::*TFuncPtr)(TArgs...),
            typename TTimer>
  class ClassTimer<void (TClass::*)(TArgs...), TFuncPtr, TTimer>
      : public ClassTimerImpl<TClass, TTimer> {
  public:
    ClassTimer(TClass& obj);

    void call(TArgs... args);
  };

  /// Specialisation for const member functions
  template <typename TClass, typename TRet, typename... TArgs,
            TRet (TClass::*TFuncPtr)(TArgs...) const, typename TTimer>
  class ClassTimer<TRet (TClass::*)(TArgs...) const, TFuncPtr, TTimer>
      : public ClassTimerImpl<TClass, TTimer> {
  public:
    ClassTimer(TClass& obj);

    TRet call(TArgs... args);
  };

  /// Specialisation for const member functions without return value
  template <typename TClass, typename... TArgs, void (TClass::*TFuncPtr)(TArgs...) const,
            typename TTimer>
  class ClassTimer<void (TClass::*)(TArgs...) const, TFuncPtr, TTimer>
      : public ClassTimerImpl<TClass, TTimer> {
  public:
    ClassTimer(TClass& obj);

    void call(TArgs... args);
  };

} // namespace corsika

#include <corsika/detail/framework/analytics/ClassTimer.inl>