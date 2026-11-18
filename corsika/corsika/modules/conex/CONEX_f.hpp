/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

//----------------------------------------------
//  C++ interface for the CONEX fortran code
//----------------------------------------------
// wrapper

#include <conexConfig.h>
#include <conexHEModels.h>
#include <array>

namespace conex {
  // the CONEX fortran interface

  extern "C" {
  extern struct { std::array<double, 16> dptl; } cxoptl_;

  //! common block for atmosphere composition
  extern struct {
    std::array<double, 3> airz, aira, airw; //!< nuclear Z, A, composition fraction
    double airavz, airava;                  //!< average Z, A
    std::array<double, 3> airi;             //!< ionization potential, not used in cxroot
  } cxair_;

  void cegs4_(int&, int&);

  void initconex_(int&, int*, int&, int&,
#ifdef CONEX_EXTENSIONS
                  int&,
#endif
                  const char*, int);
  void conexrun_(int& ipart, double& energy, double& theta, double& phi,
                 double& injectionHeight, double& dimpact, int ioseed[3]);
  void conexcascade_();
  void hadroncascade_(int&, int&, int&, int&);
  void solvemomentequations_(int&);
  void show_(int& iqi, double& ei, double& xmi, double& ymi, double& zmi, double& dmi,
             double& xi, double& yi, double& zi, double& tmi, double& ui, double& vi,
             double& wi, int& iri, double& wti, int& latchi);

  int get_number_of_depth_bins_();

  void get_shower_data_(const int&, const int&, const int&, float&, float&, float&,
                        float&, float&);
  void get_shower_edep_(const int&, const int&, float&, float&);
  void get_shower_muon_(const int&, const int&, float&, float&);
  void get_shower_gamma_(const int&, const int&, float&);
  void get_shower_electron_(const int&, const int&, float&);
  void get_shower_positron_(const int&, const int&, float&);
  void get_shower_hadron_(const int&, const int&, float&);
  }

} // namespace conex
