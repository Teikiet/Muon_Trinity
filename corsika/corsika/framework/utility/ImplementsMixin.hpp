/*
 * (c) Copyright 2021 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */
#pragma once

/**
 * @file ImplementsMixin.hpp
 */

#include <type_traits>

namespace corsika {

  namespace detail {

    /**
       Helper traits class (partial) for static compile time checking.

       This method checks whether a given class implements a particular
       mixin - this is particularly useful in the media heirarchy that utilizes
       mixins for the interface definition.
    */
    template <template <typename> typename Mixin, typename T>
    struct implements_mixin {

      template <typename U>
      static std::true_type test(Mixin<U>&);
      static std::false_type test(...);

    public:
      using type = decltype(test(std::declval<T&>()));
      static constexpr bool value = type::value;
    };

    template <template <typename> typename Mixin, typename T>
    bool constexpr implements_mixin_v = implements_mixin<Mixin, T>::value;

  } // namespace detail
} // namespace corsika
