/*
 * (c) Copyright 2023 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <FLUKA.hpp>

#include <rng_impl.hpp>

#include <cctype>
#include <string>
#include <string_view>
#include <functional>

IMPLEMENT_RNG(fluka)

namespace fluka {
  extern "C" {
  double flrndm_() { return ::draw_std_rnd(); }
  void flrnlp_(double* array, int const* N) { rng_ptr(array, *N); }

  /**
   * The following (function) pointers make sure the corresponding objects in libflukahp.a
   * won't get dropped from the final file during linking.
   */

  [[maybe_unused]] extern auto* const hepevt_ptr = &hepevt_;
  [[maybe_unused]] extern auto* const stpxyc_ptr = &stpxyz_;
  [[maybe_unused]] extern auto* const evtxyz_ptr = &evtxyz_;
  [[maybe_unused]] extern auto* const sgmxyz_ptr = &sgmxyz_;
  [[maybe_unused]] extern auto* const ndmhep_ptr = &ndmhep_;

#ifdef FLUKA_HAS_FLKAVR
  //! support for querying (INFN-) FLUKA version was added in 2024.1.0
  extern struct {
    int mjflvr, mnflvr, mrflvr;
  } flkavr_;

  extern struct {
    char chflvr;
  } flkavc_;
#endif
  }

  char const* get_version() {
#ifdef FLUKA_HAS_FLKAVR
    std::string_view const dot = ".";
    static std::string const str =
        std::to_string(flkavr_.mjflvr) + std::string{dot} +
        std::to_string(flkavr_.mnflvr) + std::string{dot} +
        std::to_string(flkavr_.mrflvr) +
        (std::isalnum(flkavc_.chflvr) ? std::string{flkavc_.chflvr} : std::string{});

    return str.c_str();
#elif defined FLUKA_EXTRACTED_VERSION
#define XSTR(x) STR(x)
#define STR(x) #x
    return XSTR(FLUKA_EXTRACTED_VERSION);
#else
    return "undefined";
#endif
  }

} // namespace fluka
