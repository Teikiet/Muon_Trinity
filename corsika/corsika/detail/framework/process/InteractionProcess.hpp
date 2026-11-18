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
   *  traits test for InteractionProcess::doInteraction method.
   */

  template <class TProcess, typename TReturn, typename TTemplate, typename... TArgs>
  struct has_method_doInteract : public detail::has_method_signature<TReturn, TArgs...> {

    ///! method signature
    using detail::has_method_signature<TReturn, TArgs...>::testSignature;

    //! the default value
    template <class T>
    static std::false_type test(...);

    //! signature of templated method
    template <class T>
    static decltype(testSignature(&T::template doInteraction<TTemplate>)) test(
        std::nullptr_t);

    //! signature of non-templated method
    template <class T>
    static decltype(testSignature(&T::doInteraction)) test(std::nullptr_t);

  public:
    /**
     *  @name traits results
     *  @{
     */
    using type = decltype(test<std::decay_t<TProcess>>(nullptr));
    static const bool value = type::value;
    //! @}
  };

  //! value traits type
  template <class TProcess, typename TReturn, typename TTemplate, typename... TArgs>
  bool constexpr has_method_doInteract_v =
      has_method_doInteract<TProcess, TReturn, TTemplate, TArgs...>::value;

  /**
   *  traits test for TEMPLATED InteractionProcess::getCrossSection method (PROPOSAL).
   */

  template <class TProcess, typename TReturn, typename TTemplate, typename... TArgs>
  struct has_method_getCrossSectionTemplate
      : public detail::has_method_signature<TReturn, TArgs...> {

    ///! method signature
    using detail::has_method_signature<TReturn, TArgs...>::testSignature;

    //! the default value
    template <class T>
    static std::false_type test(...);

    //! templated parameter option
    template <class T>
    static decltype(testSignature(&T::template getCrossSection<TTemplate>)) test(
        std::nullptr_t);

    //! non templated parameter option
    template <class T>
    static decltype(testSignature(&T::getCrossSection)) test(std::nullptr_t);

  public:
    /**
     *  @name traits results
     *  @{
     */
    using type = decltype(test<std::decay_t<TProcess>>(nullptr));
    static const bool value = type::value;
    //! @}
  };

  //! value traits type shortcut
  template <class TProcess, typename TReturn, typename TTemplate, typename... TArgs>
  bool constexpr has_method_getCrossSectionTemplate_v =
      has_method_getCrossSectionTemplate<TProcess, TReturn, TTemplate, TArgs...>::value;

  /**
   *  traits test for InteractionProcess::getCrossSection method.
   */

  template <class TProcess, typename TReturn, typename... TArgs>
  struct has_method_getCrossSection
      : public detail::has_method_signature<TReturn, TArgs...> {

    ///! method signature
    using detail::has_method_signature<TReturn, TArgs...>::testSignature;

    //! the default value
    template <class T>
    static std::false_type test(...);

    //! templated parameter option
    template <class T>
    static decltype(testSignature(&T::template getCrossSection<TArgs...>)) test(
        std::nullptr_t);

    //! non templated parameter option
    template <class T>
    static decltype(testSignature(&T::getCrossSection)) test(std::nullptr_t);

  public:
    /**
     *  @name traits results
     * @{
     */
    using type = decltype(test<std::decay_t<TProcess>>(nullptr));
    static const bool value = type::value;
    //! @}
  };

  //! value traits type shortcut
  template <class TProcess, typename TReturn, typename... TArgs>
  bool constexpr has_method_getCrossSection_v =
      has_method_getCrossSection<TProcess, TReturn, TArgs...>::value;

} // namespace corsika
