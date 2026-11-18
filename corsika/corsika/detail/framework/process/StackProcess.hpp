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
     traits test for StackProcess::doStack method
  */

  template <class TProcess, typename TReturn, typename... TArgs>
  struct has_method_doStack : public detail::has_method_signature<TReturn, TArgs...> {

    //! method signature
    using detail::has_method_signature<TReturn, TArgs...>::testSignature;

    //! the default value
    template <class T>
    static std::false_type test(...);

    //! templated parameter option
    template <class T>
    static decltype(testSignature(&T::template doStack<TArgs...>)) test(std::nullptr_t);

    //! templated-template parameter option
    template <template <typename> typename T>
    static decltype(testSignature(&T<TArgs...>::template doStack)) test(std::nullptr_t);

    //! non-templated parameter option
    template <class T>
    static decltype(testSignature(&T::doStack)) test(std::nullptr_t);

  public:
    /**
        @name traits results
        @{
    */
    using type = decltype(test<std::decay_t<TProcess>>(nullptr));
    static const bool value = type::value;
    //! @}
  };

  //! @file StackProcess.hpp
  //! value traits type
  template <class TProcess, typename TReturn, typename... TArgs>
  bool constexpr has_method_doStack_v =
      has_method_doStack<TProcess, TReturn, TArgs...>::value;

} // namespace corsika
