/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <array>

#include <rng_decl.hpp>

/**
 * \file urqmd.hpp
 *
 * Interface file for the urqmd library.
 */

DECLARE_RNG(urqmd)

namespace urqmd {

  namespace constants {
    // from coms.f
    int constexpr nmax = 500;

    // from options.f
    int constexpr numcto = 400;
    int constexpr numctp = 400;

    // from inputs.f
    int constexpr aamax = 300;

  } // namespace constants

  template <typename T>
  using nmaxArray = std::array<T, constants::nmax>;
  using nmaxIntArray = nmaxArray<int>;
  using nmaxDoubleArray = nmaxArray<double>;

  extern "C" {

  void urqmd_(int&);
  void iniurqmdc8_();
  double ranf_(int&);
  void cascinit_(int const&, int const&, int const&);
  double nucrad_(int const&);
  int pdgid_(int&, int&);
  double sigtot_(int&, int&, double&);

  // defined in coms.f
  extern struct {
    int npart, nbar, nmes, ctag, nsteps, uid_cnt, ranseed, event;
    int Ap; // projectile mass number (in case of nucleus)
    int At; // target mass number (in case of nucleus)
    int Zp; // projectile charge number (in case of nucleus)
    int Zt; // target charge number (in case of nucleus)
    int eos, dectag, NHardRes, NSoftRes, NDecRes, NElColl, NBlColl;
  } sys_;

  extern struct {
    double time, acttime, bdist, bimp, bmin;
    double ebeam; // lab-frame energy of projectile
    double ecm;
  } rsys_;

  // defined in coms.f
  extern struct {
    nmaxIntArray spin, ncoll, charge, ityp, lstcoll, iso3, origin, strid, uid;
  } isys_;

  // defined in coor.f
  extern struct {
    nmaxDoubleArray r0, rx, ry, rz, p0, px, py, pz, fmass, rww, dectime;
  } coor_;

  // defined in inputs.f
  extern struct {
    int nevents;
    std::array<int, 2> spityp; // particle codes of: [0]: projectile, [1]: target
    int prspflg;               // projectile special flag
    int trspflg; // target special flag, set to 1 unless target is nucleus > H
    std::array<int, 2> spiso3; // particle codes of: [0]: projectile, [1]: target
    int outsteps, bflag, srtflag, efuncflag, nsrt, npb, firstev;
  } inputs_;

  // defined in inputs.f
  extern struct {
    double srtmin, srtmax, pbeam, betann, betatar, betapro, pbmin, pbmax;
  } input2_;

  // defined in options.f
  extern struct {
    std::array<double, constants::numcto> CTOption;
    std::array<double, constants::numctp> CTParam;
  } options_;

  extern struct {
    int fixedseed, bf13, bf14, bf15, bf16, bf17, bf18, bf19, bf20;
  } loptions_;

  // defined in urqmdInterface.F
  extern struct {
    std::array<double, 3> xs, bim;
  } cxs_u2_;
  }

} // namespace urqmd
