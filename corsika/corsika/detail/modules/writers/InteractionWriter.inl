/*
 * (c) Copyright 2024 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/Logging.hpp>

namespace corsika {

  template <typename TTracking, typename TOutput>
  inline InteractionWriter<TTracking, TOutput>::InteractionWriter(
      media::ShowerAxis const& axis, ObservationPlane<TTracking, TOutput> const& obsPlane)
      : obsPlane_(obsPlane)
      , showerAxis_(axis)
      , interactionCounter_(0)
      , showerId_(0) {}

  template <typename TTracking, typename TOutput>
  template <typename TStackView>
  inline void InteractionWriter<TTracking, TOutput>::doSecondaries(TStackView& vS) {
    // Dumps the stack to the output parquet stream
    // The primary and secondaries are all written along with the slant depth

    if (interactionCounter_++) { return; } // Only run on the first interaction

    auto primary = vS.getProjectile();
    auto dX = showerAxis_.getProjectedX(primary.getPosition());
    CORSIKA_LOG_INFO("First interaction at dX {}", dX);
    CORSIKA_LOG_INFO("Primary: {}, E_kin {}", primary.getPID(),
                     primary.getKineticEnergy());

    // get the location of the primary w.r.t. observation plane
    Vector const displacement = primary.getPosition() - obsPlane_.getPlane().getCenter();

    auto const x = displacement.dot(obsPlane_.getXAxis());
    auto const y = displacement.dot(obsPlane_.getYAxis());
    auto const z = displacement.dot(obsPlane_.getPlane().getNormal());
    auto const nx = primary.getDirection().dot(obsPlane_.getXAxis());
    auto const ny = primary.getDirection().dot(obsPlane_.getYAxis());
    auto const nz = primary.getDirection().dot(obsPlane_.getPlane().getNormal());

    auto const px = primary.getMomentum().dot(obsPlane_.getXAxis());
    auto const py = primary.getMomentum().dot(obsPlane_.getYAxis());
    auto const pz = primary.getMomentum().dot(obsPlane_.getPlane().getNormal());

    summary_["shower_" + std::to_string(showerId_)]["pdg"] =
        static_cast<int>(get_PDG(primary.getPID()));
    summary_["shower_" + std::to_string(showerId_)]["name"] =
        static_cast<std::string>(get_name(primary.getPID()));
    summary_["shower_" + std::to_string(showerId_)]["total_energy"] =
        (primary.getKineticEnergy() + get_mass(primary.getPID())) / 1_GeV;
    summary_["shower_" + std::to_string(showerId_)]["kinetic_energy"] =
        primary.getKineticEnergy() / 1_GeV;
    summary_["shower_" + std::to_string(showerId_)]["x"] = x / 1_m;
    summary_["shower_" + std::to_string(showerId_)]["y"] = y / 1_m;
    summary_["shower_" + std::to_string(showerId_)]["z"] = z / 1_m;
    summary_["shower_" + std::to_string(showerId_)]["nx"] = static_cast<double>(nx);
    summary_["shower_" + std::to_string(showerId_)]["ny"] = static_cast<double>(ny);
    summary_["shower_" + std::to_string(showerId_)]["nz"] = static_cast<double>(nz);
    summary_["shower_" + std::to_string(showerId_)]["px"] =
        static_cast<double>(px / 1_GeV);
    summary_["shower_" + std::to_string(showerId_)]["py"] =
        static_cast<double>(py / 1_GeV);
    summary_["shower_" + std::to_string(showerId_)]["pz"] =
        static_cast<double>(pz / 1_GeV);
    summary_["shower_" + std::to_string(showerId_)]["time"] = primary.getTime() / 1_s;
    summary_["shower_" + std::to_string(showerId_)]["slant_depth"] =
        dX / (1_g / 1_cm / 1_cm);

    uint nSecondaries = 0;

    // Loop through secondaries
    auto particle = vS.begin();
    while (particle != vS.end()) {

      // Get the momentum of the secondary, write w.r.t. observation plane
      auto const p_2nd = particle.getMomentum();

      *(output_.getWriter())
          << showerId_ << static_cast<int>(get_PDG(particle.getPID()))
          << static_cast<float>(p_2nd.dot(obsPlane_.getXAxis()) / 1_GeV)
          << static_cast<float>(p_2nd.dot(obsPlane_.getYAxis()) / 1_GeV)
          << static_cast<float>(p_2nd.dot(obsPlane_.getPlane().getNormal()) / 1_GeV)
          << parquet::EndRow;

      CORSIKA_LOG_INFO(" 2ndary: {}, E_kin {}", particle.getPID(),
                       particle.getKineticEnergy());
      ++particle;
      ++nSecondaries;
    }

    summary_["shower_" + std::to_string(showerId_)]["n_secondaries"] = nSecondaries;
  }

  template <typename TTracking, typename TOutput>
  inline void InteractionWriter<TTracking, TOutput>::startOfLibrary(
      boost::filesystem::path const& directory) {
    output_.initStreamer((directory / ("interactions.parquet")).string());

    // enable compression with the default level
    output_.enableCompression();

    output_.addField("pdg", parquet::Repetition::REQUIRED, parquet::Type::INT32,
                     parquet::ConvertedType::INT_32);
    output_.addField("px", parquet::Repetition::REQUIRED, parquet::Type::FLOAT,
                     parquet::ConvertedType::NONE);
    output_.addField("py", parquet::Repetition::REQUIRED, parquet::Type::FLOAT,
                     parquet::ConvertedType::NONE);
    output_.addField("pz", parquet::Repetition::REQUIRED, parquet::Type::FLOAT,
                     parquet::ConvertedType::NONE);

    output_.buildStreamer();

    showerId_ = 0;
    interactionCounter_ = 0;
    summary_ = YAML::Node();
  }

  template <typename TTracking, typename TOutput>
  inline void InteractionWriter<TTracking, TOutput>::startOfShower(
      unsigned int const showerId) {
    showerId_ = showerId;
    interactionCounter_ = 0;
  }

  template <typename TTracking, typename TOutput>
  inline void InteractionWriter<TTracking, TOutput>::endOfShower(unsigned int const) {}

  template <typename TTracking, typename TOutput>
  inline void InteractionWriter<TTracking, TOutput>::endOfLibrary() {
    output_.closeStreamer();
  }

  template <typename TTracking, typename TOutput>
  inline YAML::Node InteractionWriter<TTracking, TOutput>::getConfig() const {
    YAML::Node node;
    node["type"] = "Interactions";
    node["units"]["energy"] = "GeV";
    node["units"]["length"] = "m";
    node["units"]["time"] = "ns";
    node["units"]["grammage"] = "g/cm^2";
    return node;
  }

  template <typename TTracking, typename TOutput>
  inline YAML::Node InteractionWriter<TTracking, TOutput>::getSummary() const {
    return summary_;
  }

} // namespace corsika
