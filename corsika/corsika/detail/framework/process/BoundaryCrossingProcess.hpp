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
     traits test for BoundaryCrossingProcess::doBoundaryCrossing method
  */
  template <class TProcess, typename TReturn, typename TParticle>
  struct has_method_doBoundaryCrossing
      : public detail::has_method_signature<TReturn, TParticle&,
                                            typename TParticle::node_type const&,
                                            typename TParticle::node_type const&> {

    ///! method signature
    using detail::has_method_signature<
        TReturn, TParticle&, typename TParticle::node_type const&,
        typename TParticle::node_type const&>::testSignature;

    //! the default value
    template <class T>
    static std::false_type test(...);

    //! templated parameter option
    template <class T>
    static decltype(testSignature(&T::template doBoundaryCrossing<TParticle>)) test(
        std::nullptr_t);

    //! non templated parameter option
    template <class T>
    static decltype(testSignature(&T::doBoundaryCrossing)) test(std::nullptr_t);

  public:
    /**
        @name traits results
        @{
    */
    using type = decltype(test<std::decay_t<TProcess>>(nullptr));
    static const bool value = type::value;
    //! @}
  };

  //! @file BoundaryCrossingProcess.hpp
  //! value traits type

  template <class TProcess, typename TReturn, typename TParticle>
  bool constexpr has_method_doBoundaryCrossing_v =
      has_method_doBoundaryCrossing<TProcess, TReturn, TParticle>::value;

} // namespace corsika
