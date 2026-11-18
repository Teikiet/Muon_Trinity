/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

/**
 * \file ProcessReturn.hpp
 **/

namespace corsika {

  /**
     @ingroup
     @{

     since in a process sequence many status updates can accumulate
     for a single particle, this enum should define only bit-flags
     that can be accumulated easily with "|="
   */

  enum class ProcessReturn : int {
    Ok = (1 << 0),
    ParticleAbsorbed = (1 << 2),
    Interacted = (1 << 3),
    Decayed = (1 << 4),
  };

  inline ProcessReturn operator|(ProcessReturn a, ProcessReturn b) {
    return static_cast<ProcessReturn>(static_cast<int>(a) | static_cast<int>(b));
  }

  inline ProcessReturn& operator|=(ProcessReturn& a, const ProcessReturn b) {
    return a = a | b;
  }

  inline ProcessReturn operator&(const ProcessReturn a, const ProcessReturn b) {
    return static_cast<ProcessReturn>(static_cast<int>(a) & static_cast<int>(b));
  }

  inline bool operator==(const ProcessReturn a, const ProcessReturn b) {
    return (static_cast<int>(a) & static_cast<int>(b)) != 0;
  }

  inline bool isOk(const ProcessReturn a) {
    return static_cast<int>(a & ProcessReturn::Ok);
  }

  inline bool isAbsorbed(const ProcessReturn a) {
    return static_cast<int>(a & ProcessReturn::ParticleAbsorbed);
  }

  inline bool isDecayed(const ProcessReturn a) {
    return static_cast<int>(a & ProcessReturn::Decayed);
  }

  inline bool isInteracted(const ProcessReturn a) {
    return static_cast<int>(a & ProcessReturn::Interacted);
  }

  //! @}

} // namespace corsika
