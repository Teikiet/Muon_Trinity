/*
 * (c) Copyright 2025 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <boost/filesystem.hpp>

#include <catch2/catch_all.hpp>

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/Vector.hpp>
#include <corsika/framework/process/InteractionCounter.hpp>
#include <corsika/modules/writers/InteractionHistogramWriter.hpp>

#include <numeric>
#include <string>

#include <yaml-cpp/yaml.h>

using namespace corsika;

namespace {

  struct DummyProcess {
    CrossSectionType getCrossSection(Code const, Code const, FourMomentum const&,
                                     FourMomentum const&) {
      return 100_mb;
    }

    template <typename TParticle>
    void doInteraction(TParticle&, Code const, Code const, FourMomentum const&,
                       FourMomentum const&) {}
  };

  struct DummyOutput {};

  std::string resetOutputDir(std::string const& path) {
    if (boost::filesystem::exists(path)) { boost::filesystem::remove_all(path); }
    boost::filesystem::create_directory(path);
    return path;
  }

} // namespace

TEST_CASE("InteractionHistogramWriter", "process") {
  logging::set_level(logging::level::info);

  DummyProcess process;
  InteractionCounter counter(process);
  DummyOutput stackView;

  auto const fillNucleusOxygen = [&]() {
    unsigned short constexpr A = 14, Z = 7;
    Code const pid = get_nucleus_code(A, Z);
    counter.doInteraction(stackView, pid, Code::Oxygen,
                          {sqrt(static_pow<2>(105_TeV) + static_pow<2>(get_mass(pid))),
                           {105_TeV, 0_GeV, 0_GeV}},
                          {Oxygen::mass, {0_eV, 0_eV, 0_eV}});
  };

  SECTION("per_shower_reset_and_parquet") {
    std::string const outputDir = "./output_interaction_hist_test_per_shower";
    resetOutputDir(outputDir);

    InteractionHistogramWriter writer([&]() { return counter.getHistogram(); },
                                      [&]() { counter.resetHistogram(); });

    writer.startOfLibrary(outputDir);

    // startOfShower: reset before this shower
    writer.startOfShower(0);
    fillNucleusOxygen();

    auto const& lab0 = counter.getHistogram().labHist();
    CHECK(std::accumulate(lab0.cbegin(), lab0.cend(), 0) == 1);

    writer.endOfShower(0);

    // Second shower: reset must clear so we only ever see one interaction per write
    writer.startOfShower(1);
    CHECK(std::accumulate(counter.getHistogram().labHist().cbegin(),
                          counter.getHistogram().labHist().cend(), 0) == 0);
    fillNucleusOxygen();
    CHECK(std::accumulate(counter.getHistogram().labHist().cbegin(),
                          counter.getHistogram().labHist().cend(), 0) == 1);

    writer.endOfShower(1);
    writer.endOfLibrary();

    CHECK(boost::filesystem::exists(outputDir + "/lab_energy_0.parquet"));
    CHECK(boost::filesystem::exists(outputDir + "/cms_energy_0.parquet"));
    CHECK(boost::filesystem::exists(outputDir + "/lab_energy_1.parquet"));
    CHECK(boost::filesystem::exists(outputDir + "/cms_energy_1.parquet"));

    CHECK(boost::filesystem::exists(boost::filesystem::path(outputDir) / "config.yaml"));
    {
      auto const cfgPath = boost::filesystem::path(outputDir) / "config.yaml";
      YAML::Node const fileCfg = YAML::LoadFile(cfgPath.string());
      CHECK(fileCfg["type"].as<std::string>() == "InteractionHistogram");
      CHECK(fileCfg["name"].as<std::string>() ==
            boost::filesystem::path(outputDir).filename().string());
    }

    auto const cfg = writer.getConfig();
    CHECK(cfg["format"].as<std::string>() == "parquet");
    CHECK(cfg["per_shower_histogram"].as<bool>());

    auto const summary = writer.getSummary();
    CHECK(summary["total_rows"].as<unsigned long long>() > 0ULL);
    CHECK(summary["shower_0"]["rows"].as<unsigned long long>() > 0ULL);
    CHECK(summary["shower_1"]["rows"].as<unsigned long long>() > 0ULL);
    CHECK(summary["shower_0"]["projectile_pdgs_lab_energy"].IsSequence());
    CHECK(summary["shower_0"]["projectile_pdgs_cms_energy"].IsSequence());
    CHECK(summary["shower_0"]["integrated_interactions_per_pdg"].IsMap());
    CHECK(summary["shower_0"]["integrated_interactions_per_pdg"].size() > 0);

    if (boost::filesystem::exists(outputDir)) {
      boost::filesystem::remove_all(outputDir);
    }
  }

  SECTION("resetHistogram_clears_bins") {
    fillNucleusOxygen();
    CHECK(std::accumulate(counter.getHistogram().labHist().cbegin(),
                          counter.getHistogram().labHist().cend(), 0) == 1);

    counter.resetHistogram();

    CHECK(std::accumulate(counter.getHistogram().labHist().cbegin(),
                          counter.getHistogram().labHist().cend(), 0) == 0);
    CHECK(std::accumulate(counter.getHistogram().CMSHist().cbegin(),
                          counter.getHistogram().CMSHist().cend(), 0) == 0);
  }

  SECTION("no_reset_callback_keeps_cumulative_histogram") {
    std::string const outputDir = "./output_interaction_hist_test_cumulative";
    InteractionHistogramWriter writer([&]() { return counter.getHistogram(); });

    writer.startOfLibrary(resetOutputDir(outputDir));
    writer.startOfShower(0);
    fillNucleusOxygen();
    writer.endOfShower(0);

    writer.startOfShower(1);
    fillNucleusOxygen();
    writer.endOfShower(1);
    writer.endOfLibrary();

    auto const cfg = writer.getConfig();
    CHECK(!cfg["per_shower_histogram"].as<bool>());

    if (boost::filesystem::exists(outputDir)) {
      boost::filesystem::remove_all(outputDir);
    }
  }
}
