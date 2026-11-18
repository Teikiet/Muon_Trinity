#pragma once

/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

/**
 * @file HasMethodSignature.hpp
 */

#include <type_traits>

namespace corsika {

  namespace detail {

    /**
       Helper traits class (partial) for static compile time checking.

       Note, this is a poor replacement for C++20 concepts... they are
       eagerly awaited!

       It defines the default body of a generic test function returning
       std::false_type.

       In addition it defines the pattern for class-method matching with a
       return type TReturn and function arguments TArgs... . Right now
       both method signatures, "const" and "not const", are matched.
    */
    template <typename TReturn, typename... TArgs>
    struct has_method_signature {

      // the non-const version
      template <class T>
      static std::true_type testSignature(TReturn (T::*)(TArgs...));

      // the const version
      template <class T>
      static std::true_type testSignature(TReturn (T::*)(TArgs...) const);
    };

  } // namespace detail
} // namespace corsika
