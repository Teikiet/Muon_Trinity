/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

namespace corsika {
  namespace media {

    namespace detail {

      struct NoExtraModelInner {};

      template <typename M>
      struct NoExtraModel {};

      template <template <typename> typename M>
      struct has_extra_models : std::true_type {};

      template <>
      struct has_extra_models<NoExtraModel> : std::false_type {};

    } // namespace detail
  }   // namespace media
} // namespace corsika
