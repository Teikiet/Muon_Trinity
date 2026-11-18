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

namespace corsika {

  // Common

  template <typename TClass, typename TRet, typename... TArgs,
            TRet (TClass::*TFuncPtr)(TArgs...), typename TTimer>
  inline ClassTimer<TRet (TClass::*)(TArgs...), TFuncPtr, TTimer>::ClassTimer(TClass& obj)
      : ClassTimerImpl<TClass, TTimer>(obj) {}

  template <typename TClass, typename TRet, typename... TArgs,
            TRet (TClass::*TFuncPtr)(TArgs...), typename TTimer>
  inline TRet ClassTimer<TRet (TClass::*)(TArgs...), TFuncPtr, TTimer>::call(
      TArgs... args) {
    this->start_ = ClassTimerImpl<TClass, TTimer>::clock_type::now();
    auto tmp = (this->obj_.*TFuncPtr)(std::forward<TArgs>(args)...);
    this->timeDiff_ = std::chrono::duration_cast<
        typename ClassTimerImpl<TClass, TTimer>::duration_type>(
        ClassTimerImpl<TClass, TTimer>::clock_type::now() - this->start_);
    return tmp;
  }

  // Specialisation 1

  template <typename TClass, typename... TArgs, void (TClass::*TFuncPtr)(TArgs...),
            typename TTimer>
  inline ClassTimer<void (TClass::*)(TArgs...), TFuncPtr, TTimer>::ClassTimer(TClass& obj)
      : ClassTimerImpl<TClass, TTimer>(obj) {}

  template <typename TClass, typename... TArgs, void (TClass::*TFuncPtr)(TArgs...),
            typename TTimer>
  inline void ClassTimer<void (TClass::*)(TArgs...), TFuncPtr, TTimer>::call(
      TArgs... args) {
    this->start_ = ClassTimerImpl<TClass, TTimer>::clock_type::now();
    (this->obj_.*TFuncPtr)(std::forward<TArgs>(args)...);
    this->timeDiff_ = std::chrono::duration_cast<
        typename ClassTimerImpl<TClass, TTimer>::duration_type>(
        ClassTimerImpl<TClass, TTimer>::clock_type::now() - this->start_);
    return;
  }

  /// Specialisation 2

  template <typename TClass, typename TRet, typename... TArgs,
            TRet (TClass::*TFuncPtr)(TArgs...) const, typename TTimer>
  inline ClassTimer<TRet (TClass::*)(TArgs...) const, TFuncPtr, TTimer>::ClassTimer(
      TClass& obj)
      : ClassTimerImpl<TClass, TTimer>(obj) {}

  template <typename TClass, typename TRet, typename... TArgs,
            TRet (TClass::*TFuncPtr)(TArgs...) const, typename TTimer>
  inline TRet ClassTimer<TRet (TClass::*)(TArgs...) const, TFuncPtr, TTimer>::call(
      TArgs... args) {
    this->start_ = ClassTimerImpl<TClass, TTimer>::clock_type::now();
    auto tmp = (this->obj_.*TFuncPtr)(std::forward<TArgs>(args)...);
    this->timeDiff_ = std::chrono::duration_cast<
        typename ClassTimerImpl<TClass, TTimer>::duration_type>(
        ClassTimerImpl<TClass, TTimer>::clock_type::now() - this->start_);
    return tmp;
  }

  /// Specialisation 3
  template <typename TClass, typename... TArgs, void (TClass::*TFuncPtr)(TArgs...) const,
            typename TTimer>
  inline ClassTimer<void (TClass::*)(TArgs...) const, TFuncPtr, TTimer>::ClassTimer(
      TClass& obj)
      : ClassTimerImpl<TClass, TTimer>(obj) {}

  template <typename TClass, typename... TArgs, void (TClass::*TFuncPtr)(TArgs...) const,
            typename TTimer>
  inline void ClassTimer<void (TClass::*)(TArgs...) const, TFuncPtr, TTimer>::call(
      TArgs... args) {
    this->start_ = ClassTimerImpl<TClass, TTimer>::clock_type::now();
    (this->obj_.*TFuncPtr)(std::forward<TArgs>(args)...);
    this->timeDiff_ = std::chrono::duration_cast<
        typename ClassTimerImpl<TClass, TTimer>::duration_type>(
        ClassTimerImpl<TClass, TTimer>::clock_type::now() - this->start_);
    return;
  }

} // namespace corsika
