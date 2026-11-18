/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/modules/Random.hpp>
#include <corsika/modules/urqmd/UrQMD.hpp>
#include <corsika/modules/urqmd/ParticleConversion.hpp>

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/core/EnergyMomentumOperations.hpp>
#include <corsika/framework/geometry/QuantityVector.hpp>
#include <corsika/framework/geometry/Vector.hpp>
#include <corsika/framework/utility/COMBoost.hpp>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/multi_array.hpp>

#include <algorithm>
#include <functional>
#include <iostream>
#include <fstream>
#include <sstream>

#include <urqmd.hpp>

namespace corsika::urqmd {

  inline UrQMD::UrQMD(boost::filesystem::path xs_file, int const retryFlag)
      : iflb_(retryFlag) {
    readXSFile(xs_file);
    corsika::connect_random_stream(RNG_, ::urqmd::set_rng_function);
    ::urqmd::iniurqmdc8_();
  }

  inline bool UrQMD::isValid(Code const projectileId, Code const targetId) const {

    if (!is_hadron(projectileId) || !corsika::urqmd::canInteract(projectileId)) {
      return false;
    }
    if (!is_nucleus(targetId)) { return false; }
    return true;
  }

  inline CrossSectionType UrQMD::getTabulatedCrossSection(
      Code const projectileId, Code const targetId, HEPEnergyType const labEnergy) const {

    // translated to C++ from CORSIKA 7 subroutine cxtot_u

    auto const kinEnergy = labEnergy - get_mass(projectileId);

    if (kinEnergy < HEPEnergyType::zero()) {
      CORSIKA_LOG_ERROR("Kinetic energy {} below zero", kinEnergy);
      throw std::runtime_error("Negative kinetic energy");
    }

    double const logKinEnergy = std::log10(kinEnergy * (1 / 1_GeV));
    double const ye = std::max(10 * logKinEnergy + 10.5, 1.);
    int const je = std::min(int(ye), int(xs_interp_support_table_.shape()[2] - 2));
    std::array<double, 3> w;
    w[2 - 1] = ye - je;
    w[3 - 1] = w[2 - 1] * (w[2 - 1] - 1.) * .5;
    w[1 - 1] = 1 - w[2 - 1] + w[3 - 1];
    w[2 - 1] = w[2 - 1] - 2 * w[3 - 1];

    int projectileIndex;
    switch (projectileId) {
      case Code::Proton:
        projectileIndex = 0;
        break;
      case Code::AntiProton:
        projectileIndex = 1;
        break;
      case Code::Neutron:
        projectileIndex = 2;
        break;
      case Code::AntiNeutron:
        projectileIndex = 3;
        break;
      case Code::PiPlus:
        projectileIndex = 4;
        break;
      case Code::PiMinus:
        projectileIndex = 5;
        break;
      case Code::KPlus:
        projectileIndex = 6;
        break;
      case Code::KMinus:
        projectileIndex = 7;
        break;
      case Code::K0Short:
      case Code::K0Long:
      /* since K0Short and K0Long are treated the same, we can also add K0 and K0Bar
       * to the list. This is a deviation from CORSIKA 7. */
      case Code::K0:
      case Code::K0Bar:
        projectileIndex = 8;
        break;
      default: { // LCOV_EXCL_START since this can never happen due to canInteract
        CORSIKA_LOG_WARN("UrQMD cross-section not tabulated for {}", projectileId);
        return CrossSectionType::zero();
        // LCOV_EXCL_STOP
      }
    }

    int targetIndex;
    switch (targetId) {
      case Code::Nitrogen:
        targetIndex = 0;
        break;
      case Code::Oxygen:
        targetIndex = 1;
        break;
      case Code::Argon:
        targetIndex = 2;
        break;
      default:
        std::stringstream ss;
        ss << "UrQMD cross-section not tabluated for target " << targetId;
        throw std::runtime_error(ss.str().data());
    }

    auto result = CrossSectionType::zero();
    for (int i = 0; i < 3; ++i) {
      result +=
          xs_interp_support_table_[projectileIndex][targetIndex][je + i - 1 - 1] * w[i];
    }

    CORSIKA_LOG_TRACE(
        "UrQMD::GetTabulatedCrossSection proj={}, targ={}, E={} GeV, sigma={}",
        get_name(projectileId), get_name(targetId), labEnergy / 1_GeV, result);

    return result;
  }

  inline CrossSectionType UrQMD::getCrossSection(Code const projectileId,
                                                 Code const targetId,
                                                 FourMomentum const& projectileP4,
                                                 FourMomentum const& targetP4) const {

    if (!isValid(projectileId, targetId)) {
      /*
       * unfortunately unavoidable at the moment until we have tools to get the actual
       * inealstic cross-section from UrQMD
       */
      return CrossSectionType::zero();
    }

    // define projectile, in lab frame
    auto const sqrtS2 = (projectileP4 + targetP4).getNormSqr();
    HEPEnergyType const Elab = (sqrtS2 - static_pow<2>(get_mass(projectileId)) -
                                static_pow<2>(get_mass(targetId))) /
                               (2 * get_mass(targetId));

    bool const tabulated = true;
    if (tabulated) { return getTabulatedCrossSection(projectileId, targetId, Elab); }

    // the following is a translation of ptsigtot() into C++
    if (!is_nucleus(projectileId) &&
        !is_nucleus(targetId)) { // both particles are "special"

      double sqrtS = sqrt(sqrtS2) / 1_GeV;

      // we must set some UrQMD globals first...
      auto const [ityp, iso3] = corsika::urqmd::convertToUrQMD(projectileId);
      ::urqmd::inputs_.spityp[0] = ityp;
      ::urqmd::inputs_.spiso3[0] = iso3;

      auto const [itypTar, iso3Tar] = corsika::urqmd::convertToUrQMD(targetId);
      ::urqmd::inputs_.spityp[1] = itypTar;
      ::urqmd::inputs_.spiso3[1] = iso3Tar;

      int one = 1;
      int two = 2;
      return ::urqmd::sigtot_(one, two, sqrtS) * 1_mb;
    }

    // at least one of them is a nucleus
    int const Ap = is_nucleus(projectileId) ? get_nucleus_A(projectileId) : 1;
    int const At = is_nucleus(targetId) ? get_nucleus_A(targetId) : 1;
    double const maxImpact = ::urqmd::nucrad_(Ap) + ::urqmd::nucrad_(At) +
                             2 * ::urqmd::options_.CTParam[30 - 1];
    return 10_mb * M_PI * static_pow<2>(maxImpact);
    // is a constant cross-section really reasonable?
  }

  template <typename TView>
  inline void UrQMD::doInteraction(TView& view, Code const projectileId,
                                   Code const targetId, FourMomentum const& projectileP4,
                                   FourMomentum const& targetP4) {

    // define projectile, in lab frame
    auto const sqrtS2 = (projectileP4 + targetP4).getNormSqr();
    HEPEnergyType const Elab = (sqrtS2 - static_pow<2>(get_mass(projectileId)) -
                                static_pow<2>(get_mass(targetId))) /
                               (2 * get_mass(targetId));

    if (!isValid(projectileId, targetId)) {
      throw std::runtime_error("invalid target,projectile,energy combination");
    }

    size_t const targetA = get_nucleus_A(targetId);
    size_t const targetZ = get_nucleus_Z(targetId);

    ::urqmd::inputs_.nevents = 1;
    ::urqmd::sys_.eos = 0; // could be configurable in principle
    ::urqmd::inputs_.outsteps = 1;
    ::urqmd::sys_.nsteps = 1;

    // initialization regarding projectile
    if (is_nucleus(projectileId)) {
      // is this everything?
      ::urqmd::inputs_.prspflg = 0;

      ::urqmd::sys_.Ap = get_nucleus_A(projectileId);
      ::urqmd::sys_.Zp = get_nucleus_Z(projectileId);
      ::urqmd::rsys_.ebeam = (Elab - get_mass(projectileId)) / 1_GeV / ::urqmd::sys_.Ap;

      ::urqmd::rsys_.bdist = ::urqmd::nucrad_(targetA) +
                             ::urqmd::nucrad_(::urqmd::sys_.Ap) +
                             2 * ::urqmd::options_.CTParam[30 - 1];

      int const id = 1;
      ::urqmd::cascinit_(::urqmd::sys_.Zp, ::urqmd::sys_.Ap, id);
    } else {
      ::urqmd::inputs_.prspflg = 1;
      ::urqmd::sys_.Ap =
          1; // even for non-baryons this has to be set, see vanilla UrQMD.f
      ::urqmd::rsys_.bdist = ::urqmd::nucrad_(targetA) + ::urqmd::nucrad_(1) +
                             2 * ::urqmd::options_.CTParam[30 - 1];
      ::urqmd::rsys_.ebeam = (Elab - get_mass(projectileId)) / 1_GeV;

      auto const [ityp, iso3] = corsika::urqmd::convertToUrQMD(
          (projectileId == Code::K0Long || projectileId == Code::K0Short)
              ? (booleanDist_(RNG_) ? Code::K0 : Code::K0Bar)
              : projectileId);
      // todo: conversion of K_long/short into strong eigenstates;
      ::urqmd::inputs_.spityp[0] = ityp;
      ::urqmd::inputs_.spiso3[0] = iso3;
    }

    // initialization regarding target
    if (is_nucleus(targetId)) {
      ::urqmd::sys_.Zt = targetZ;
      ::urqmd::sys_.At = targetA;
      ::urqmd::inputs_.trspflg = 0; // nucleus as target
      int const id = 2;
      ::urqmd::cascinit_(::urqmd::sys_.Zt, ::urqmd::sys_.At, id);
    } else {
      ::urqmd::inputs_.trspflg = 1; // special particle as target
      auto const [ityp, iso3] = corsika::urqmd::convertToUrQMD(targetId);
      ::urqmd::inputs_.spityp[1] = ityp;
      ::urqmd::inputs_.spiso3[1] = iso3;
    }

    int iflb =
        iflb_; // flag for retrying interaction in case of empty event, 0 means retry
    ::urqmd::urqmd_(iflb);

    // now retrieve secondaries from UrQMD
    COMBoost const boost(projectileP4, targetP4);

    for (int i = 0; i < ::urqmd::sys_.npart; ++i) {
      auto code = corsika::urqmd::convertFromUrQMD(::urqmd::isys_.ityp[i],
                                                   ::urqmd::isys_.iso3[i]);
      if (code == Code::K0 || code == Code::K0Bar) {
        code = booleanDist_(RNG_) ? Code::K0Short : Code::K0Long;
      }

      MomentumVector const momentum(QuantityVector<hepmomentum_d>(
          boost.getRotation() *
          QuantityVector<hepmomentum_d>{::urqmd::coor_.px[i] * 1_GeV,
                                        ::urqmd::coor_.py[i] * 1_GeV,
                                        ::urqmd::coor_.pz[i] * 1_GeV}
              .getEigenVector()));

      CORSIKA_LOG_DEBUG(" {} {} {} ", i, code, momentum.getComponents());

      HEPEnergyType const mass = get_mass(code);
      HEPEnergyType const Ekin = calculate_kinetic_energy(momentum.getNorm(), mass);
      if (Ekin <= 0_GeV) {
        if (Ekin < 0_GeV) {
          CORSIKA_LOG_WARN("Negative kinetic energy {} {}. Skipping.", code, Ekin);
        }
        view.addSecondary(std::make_tuple(code, 0_eV, DirectionVector{{0, 0, 0}}));
      } else {
        view.addSecondary(std::make_tuple(code, Ekin, momentum.normalized()));
      }
    }
    CORSIKA_LOG_DEBUG("UrQMD generated {} secondaries!", ::urqmd::sys_.npart);
  }

  inline void UrQMD::readXSFile(boost::filesystem::path const filename) {
    boost::filesystem::ifstream file(filename, std::ios::in);

    if (!file.is_open()) {
      // LCOV_EXCL_START since this is pointless to test
      throw std::runtime_error(filename.native() + " could not be opened.");
      // LCOV_EXCL_STOP
    }

    std::string line;

    std::getline(file, line);
    std::stringstream ss(line);

    char dummy; // this is '#'
    int nTargets, nProjectiles, nSupports;
    ss >> dummy >> nTargets >> nProjectiles >> nSupports;

    decltype(xs_interp_support_table_)::extent_gen extents;
    xs_interp_support_table_.resize(extents[nProjectiles][nTargets][nSupports]);

    for (int i = 0; i < nTargets; ++i) {
      for (int j = 0; j < nProjectiles; ++j) {
        for (int k = 0; k < nSupports; ++k) {
          std::getline(file, line);
          std::stringstream s(line);
          double energy, sigma;
          s >> energy >> sigma;
          xs_interp_support_table_[j][i][k] = sigma * 1_mb;
        }

        std::getline(file, line);
        std::getline(file, line);
      }
    }
    file.close();
  }

} // namespace corsika::urqmd
