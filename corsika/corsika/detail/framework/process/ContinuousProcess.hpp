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
     traits test for ContinuousProcess::doContinuous method
  */
  template <class TProcess, typename TReturn, typename TArg1, typename TArg2>
  struct has_method_doContinuous
      : public detail::has_method_signature<TReturn, TArg1, TArg2, bool> {

    //! method signature
    using detail::has_method_signature<TReturn, TArg1, TArg2, bool>::testSignature;

    //! the default value
    template <class T>
    static std::false_type test(...);

    //! templated parameter option
    template <class T>
    static decltype(testSignature(&T::template doContinuous<TArg1, TArg2>)) test(
        std::nullptr_t);

    //! non templated parameter option
    template <class T>
    static decltype(testSignature(&T::doContinuous)) test(std::nullptr_t);

  public:
    /**
        @name traits results
        @{
    */
    using type = decltype(test<std::decay_t<TProcess>>(nullptr));
    static const bool value = type::value;
    //! @}
  };

  //! @file ContinuousProcess.hpp
  //! value traits type
  template <class TProcess, typename TReturn, typename TArg1, typename TArg2>
  bool constexpr has_method_doContinuous_v =
      has_method_doContinuous<TProcess, TReturn, TArg1, TArg2>::value;

  /**
     traits test for ContinuousProcess::getMaxStepLength method
  */

  template <class TProcess, typename TReturn, typename... TArgs>
  struct has_method_getMaxStepLength
      : public detail::has_method_signature<TReturn, TArgs...> {

    //! method signature
    using detail::has_method_signature<TReturn, TArgs...>::testSignature;

    //! the default value
    template <class T>
    static std::false_type test(...);

    //! templated option
    template <class T>
    static decltype(testSignature(&T::template getMaxStepLength<TArgs...>)) test(
        std::nullptr_t);

    //! non templated option
    template <class T>
    static decltype(testSignature(&T::getMaxStepLength)) test(std::nullptr_t);

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
  bool constexpr has_method_getMaxStepLength_v =
      has_method_getMaxStepLength<TProcess, TReturn, TArgs...>::value;

} // namespace corsika
