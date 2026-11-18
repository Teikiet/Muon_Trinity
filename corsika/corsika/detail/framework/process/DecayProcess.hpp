/*
 * (c) Copyright 2021 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/process/ProcessTraits.hpp>
#include <corsika/framework/utility/HasMethodSignature.hpp>

namespace corsika {

  /**
     traits test for DecayProcess::doDecay method
  */

  template <class TProcess, typename TReturn, typename... TArgs>
  struct has_method_doDecay : public detail::has_method_signature<TReturn, TArgs...> {

    //! type of process to be studied
    typedef std::decay_t<TProcess> process_type;

    ///! method signature
    using detail::has_method_signature<TReturn, TArgs...>::testSignature;

    //! the default value
    template <class T>
    static std::false_type test(...);

    //! signature of templated method
    template <class T>
    static decltype(testSignature(&T::template doDecay<TArgs...>)) test(std::nullptr_t);

    //! signature of non-templated method
    template <class T>
    static decltype(testSignature(&T::doDecay)) test(std::nullptr_t);

  public:
    /**
        @name traits results
        @{
    */
    using type = decltype(test<process_type>(nullptr));
    static const bool value = type::value;
    //! @}
  };

  //! @file DecayProcess.hpp
  //! value traits type
  template <class TProcess, typename TReturn, typename... TArgs>
  bool constexpr has_method_doDecay_v =
      has_method_doDecay<TProcess, TReturn, TArgs...>::value;

  /**
     traits test for DecayProcess::getLifetime method
  */

  template <class TProcess, typename TReturn, typename... TArgs>
  struct has_method_getLifetime : public detail::has_method_signature<TReturn, TArgs...> {

    //! the moethod signature
    using detail::has_method_signature<TReturn, TArgs...>::testSignature;

    //! the default value
    template <class T>
    static std::false_type test(...);

    //! signature of templated method
    template <class T>
    static decltype(testSignature(&T::template getLifetime<TArgs...>)) test(
        std::nullptr_t);

    //! signature of non-templated method
    template <class T>
    static decltype(testSignature(&T::getLifetime)) test(std::nullptr_t);

  public:
    /**
        @name traits results
        @{
    */
    using type = decltype(test<std::decay_t<TProcess>>(nullptr));
    static const bool value = type::value;
    //! @}
  };

  //! value traits type
  template <class TProcess, typename TReturn, typename... TArgs>
  bool constexpr has_method_getLifetime_v =
      has_method_getLifetime<TProcess, TReturn, TArgs...>::value;

} // namespace corsika
