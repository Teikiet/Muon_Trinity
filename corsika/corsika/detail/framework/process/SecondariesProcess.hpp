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
   * traits test for SecondariesProcess::doSecondaries method.
   */
  template <class TProcess, typename TReturn, typename... TArg>
  struct has_method_doSecondaries
      : public detail::has_method_signature<TReturn, TArg...> {

    //! method signature
    using detail::has_method_signature<TReturn, TArg...>::testSignature;

    //! the default value
    template <class T>
    static std::false_type test(...);

    //! templated parameter option
    template <class T>
    static decltype(testSignature(&T::template doSecondaries<TArg...>)) test(
        std::nullptr_t);

    //! non templated parameter option
    template <class T>
    static decltype(testSignature(&T::doSecondaries)) test(std::nullptr_t);

  public:
    /**
     *  @name traits results
     *  @{
     */
    using type = decltype(test<std::decay_t<TProcess>>(nullptr));
    static const bool value = type::value;
    //! @}
  };

  /**
   * @file SecondariesProcess.hpp
   *
   * @brief value traits type.
   */
  template <class TProcess, typename TReturn, typename... TArg>
  bool constexpr has_method_doSecondaries_v =
      has_method_doSecondaries<TProcess, TReturn, TArg...>::value;

} // namespace corsika
