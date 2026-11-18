/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <cfenv>

/*
 * Same declaration of function as provided in GLIBC
 * Repetition allowed in the case where cfenv defines the functions already, no clash.
 */
extern "C" {

int feenableexcept(int excepts) noexcept;
int fedisableexcept(int excepts) noexcept;
}

#ifdef CORSIKA_HAS_FEENABLEEXCEPT
// Nothing to do, OS privides the functions
#else
#ifdef CORSIKA_OS_MAC
#include <corsika/detail/framework/utility/CorsikaFenvOSX.inl>
#else
#include <corsika/detail/framework/utility/CorsikaFenvFallback.inl>
#endif
#endif
