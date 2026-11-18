/*
 * (c) Copyright 2025 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/output/ParquetStreamer.hpp>

#include <boost/histogram/indexed.hpp>

#include <parquet/stream_writer.h>

#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace corsika {

  namespace detail_interaction_hist_parquet {

    template <typename Histogram>
    inline std::vector<int> pdgColumnOrder(Histogram const& h) {
      auto const& ax0 = h.axis(0);
      std::vector<int> v;
      auto const n = ax0.size();
      v.reserve(static_cast<std::size_t>(n));
      for (std::size_t i = 0; i < static_cast<std::size_t>(n); ++i) {
        v.push_back(ax0.value(static_cast<int>(i)));
      }
      return v;
    }

    template <typename Histogram>
    inline void fillTable2D(Histogram const& h, std::vector<std::vector<double>>& t,
                            std::size_t& n0, std::size_t& n1) {
      namespace bh = boost::histogram;
      auto const& ax0 = h.axis(0);
      auto const& ax1 = h.axis(1);
      n0 = static_cast<std::size_t>(ax0.size());
      n1 = static_cast<std::size_t>(ax1.size());
      t.assign(n0, std::vector<double>(n1, 0.));
      for (auto const& x : bh::indexed(h, bh::coverage::inner)) {
        t[static_cast<std::size_t>(x.index(0))][static_cast<std::size_t>(x.index(1))] =
            *x;
      }
    }

    template <typename Histogram>
    inline std::map<int, double> integratedPerPdg(
        Histogram const& h, std::vector<std::vector<double>> const& t, std::size_t n0,
        std::size_t n1) {
      auto const& ax0 = h.axis(0);
      std::map<int, double> m;
      for (std::size_t i = 0; i < n0; ++i) {
        double s = 0.;
        for (std::size_t j = 0; j < n1; ++j) { s += t[i][j]; }
        m[ax0.value(static_cast<int>(i))] = s;
      }
      return m;
    }

    template <typename Histogram>
    inline void writeParquetWide(boost::filesystem::path const& filepath,
                                 std::uint32_t const showerId, Histogram const& h,
                                 std::vector<std::vector<double>> const& t,
                                 std::size_t n0, std::size_t n1, bool labEnergy) {
      // parquet::StreamWriter closes to an invalid Parquet file (often 4 bytes: magic
      // only) if no row is ever written; do not open a file in that case.
      if (n1 == 0) { return; }

      ParquetStreamer ps;
      ps.initStreamer(filepath.string());
      ps.enableCompression();
      if (labEnergy) {
        ps.addField("energy_low_eV", parquet::Repetition::REQUIRED, parquet::Type::FLOAT,
                    parquet::ConvertedType::NONE);
        ps.addField("energy_high_eV", parquet::Repetition::REQUIRED, parquet::Type::FLOAT,
                    parquet::ConvertedType::NONE);
      } else {
        ps.addField("sqrt_s_low_eV", parquet::Repetition::REQUIRED, parquet::Type::FLOAT,
                    parquet::ConvertedType::NONE);
        ps.addField("sqrt_s_high_eV", parquet::Repetition::REQUIRED, parquet::Type::FLOAT,
                    parquet::ConvertedType::NONE);
      }
      auto const& ax0 = h.axis(0);
      for (std::size_t i = 0; i < n0; ++i) {
        int const pdg = ax0.value(static_cast<int>(i));
        std::string const col = "pdg_" + std::to_string(pdg);
        ps.addField(col, parquet::Repetition::REQUIRED, parquet::Type::FLOAT,
                    parquet::ConvertedType::NONE);
      }
      ps.buildStreamer();

      auto const& ax1 = h.axis(1);
      // Release local shared_ptr before closeStreamer() so writer_.reset() destroys the
      // StreamWriter before outfile_->Close() (same pattern as ParticleWriterParquet).
      {
        auto w = ps.getWriter();
        for (std::size_t j = 0; j < n1; ++j) {
          float const elo = ax1.bin(static_cast<int>(j)).lower();
          float const ehi = ax1.bin(static_cast<int>(j)).upper();
          *w << showerId << elo << ehi;
          for (std::size_t i = 0; i < n0; ++i) { *w << static_cast<float>(t[i][j]); }
          *w << parquet::EndRow;
        }
      }
      ps.closeStreamer();
    }

  } // namespace detail_interaction_hist_parquet

  inline void InteractionHistogramWriterParquet::startOfLibrary(
      boost::filesystem::path const& directory) {
    directory_ = directory;
    totalRows_ = 0;
    rowsLastShower_ = 0;
    last_lab_pdgs_.clear();
    last_cms_pdgs_.clear();
    last_integrated_per_pdg_.clear();
  }

  inline void InteractionHistogramWriterParquet::endOfShower(
      unsigned int const showerId, InteractionHistogram const& h) {
    using namespace detail_interaction_hist_parquet;
    auto const& labH = h.labHist();
    auto const& cmsH = h.CMSHist();

    std::vector<std::vector<double>> tLab;
    std::size_t n0l = 0;
    std::size_t n1l = 0;
    fillTable2D(labH, tLab, n0l, n1l);
    last_integrated_per_pdg_ = integratedPerPdg(labH, tLab, n0l, n1l);
    last_lab_pdgs_ = pdgColumnOrder(labH);

    std::vector<std::vector<double>> tCms;
    std::size_t n0c = 0;
    std::size_t n1c = 0;
    fillTable2D(cmsH, tCms, n0c, n1c);
    last_cms_pdgs_ = pdgColumnOrder(cmsH);

    std::string const sid = std::to_string(showerId);
    auto const uShower = static_cast<std::uint32_t>(showerId);
    writeParquetWide(directory_ / ("lab_energy_" + sid + ".parquet"), uShower, labH, tLab,
                     n0l, n1l, true);
    writeParquetWide(directory_ / ("cms_energy_" + sid + ".parquet"), uShower, cmsH, tCms,
                     n0c, n1c, false);

    rowsLastShower_ = n1l + n1c;
    totalRows_ += rowsLastShower_;
  }

  inline void InteractionHistogramWriterParquet::endOfLibrary() {}

} // namespace corsika
