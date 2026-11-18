/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/core/Logging.hpp>
#include <corsika/framework/geometry/CoordinateSystemRegistry.hpp>
#include <corsika/media/CORSIKA7Atmospheres.hpp>
#include <corsika/modules/conex/CONEXhybrid.hpp>
#include <corsika/modules/conex/CONEXrandom.hpp>
#include <corsika/modules/Random.hpp>
#include <corsika/modules/conex/CONEX_f.hpp>
#include <corsika/framework/random/RNGManager.hpp>
#include <corsika/framework/core/PhysicalConstants.hpp>
#include <corsika/framework/core/Logging.hpp>

#include <conexConfig.h>

#include <algorithm>
#include <fstream>
#include <numeric>
#include <utility>

namespace corsika {

  template <typename TOutputE, typename TOutputN>
  inline CONEXhybrid<TOutputE, TOutputN>::CONEXhybrid(
      Point const& center, media::ShowerAxis const& showerAxis, LengthType groundDist,
      LengthType injectionHeight, HEPEnergyType primaryEnergy, PDGCode primaryPDG,
      TOutputE& args1, TOutputN& args2)
      : SubWriter<TOutputE>(args1)
      , SubWriter<TOutputN>(args2)
      , center_{center}
      , showerAxis_{showerAxis}
      , groundDist_{groundDist}
      , injectionHeight_{injectionHeight}
      , primaryEnergy_{primaryEnergy}
      , primaryPDG_{primaryPDG}
      , showerCore_{showerAxis_.getStart() + showerAxis_.getDirection() * groundDist_}
      , conexObservationCS_{std::invoke([&]() {
        auto const& c8cs = get_root_CoordinateSystem();
        auto const translation = showerCore_ - center;
        auto const intermediateCS =
            make_translation(c8cs, translation.getComponents(c8cs));
        auto const transformCS = make_rotationToZ(intermediateCS, translation);

        CORSIKA_LOG_DEBUG("translation C8/CONEX obs: ", translation.getComponents());

        /*
        auto const transform = CoordinateSystem::getTransformation(
            intermediateCS2, c8cs); // either this way or vice versa... TODO: test this!
        */
        return transformCS;
      })}
      , x_sf_{std::invoke([&]() {
        Vector<length_d> const a{conexObservationCS_, {0._m, 0._m, 1._m}};
        auto b = a.cross(showerAxis_.getDirection());
        auto const lengthB = b.getNorm();
        if (lengthB < 1e-10_m) {
          b = Vector<length_d>{conexObservationCS_, {1_m, 0_m, 0_m}};
        }

        return b.normalized();
      })}
      , y_sf_{showerAxis_.getDirection().cross(x_sf_)} {

    CORSIKA_LOG_DEBUG("x_sf (conexObservationCS): {}",
                      x_sf_.getComponents(conexObservationCS_));
    CORSIKA_LOG_DEBUG("x_sf (C8): {}", x_sf_.getComponents(get_root_CoordinateSystem()));

    CORSIKA_LOG_DEBUG("y_sf (conexObservationCS): {}",
                      y_sf_.getComponents(conexObservationCS_));

    CORSIKA_LOG_DEBUG("y_sf (C8): {}", y_sf_.getComponents(get_root_CoordinateSystem()));

    CORSIKA_LOG_DEBUG("showerAxisDirection (conexObservationCS): {}",
                      showerAxis_.getDirection().getComponents(conexObservationCS_));
    CORSIKA_LOG_DEBUG(
        "showerAxisDirection (C8): {}",
        showerAxis_.getDirection().getComponents(get_root_CoordinateSystem()));

    CORSIKA_LOG_DEBUG("showerCore (conexObservationCS): {}",
                      showerCore_.getCoordinates(conexObservationCS_));
    CORSIKA_LOG_DEBUG("showerCore (C8): {}",
                      showerCore_.getCoordinates(get_root_CoordinateSystem()));

    auto const& components = ::corsika::standardAirComposition.getComponents();
    auto const& fractions = ::corsika::standardAirComposition.getFractions();
    if (::corsika::standardAirComposition.getSize() != 3) {
      throw std::runtime_error{"CONEXhybrid only usable with standard 3-component air"};
    }

    std::transform(components.cbegin(), components.cend(), ::conex::cxair_.aira.begin(),
                   get_nucleus_A);
    std::transform(components.cbegin(), components.cend(), ::conex::cxair_.airz.begin(),
                   get_nucleus_Z);
    std::copy(fractions.cbegin(), fractions.cend(), ::conex::cxair_.airw.begin());

    ::conex::cxair_.airava =
        std::inner_product(::conex::cxair_.airw.cbegin(), ::conex::cxair_.airw.cend(),
                           ::conex::cxair_.aira.cbegin(), 0.);
    ::conex::cxair_.airavz =
        std::inner_product(::conex::cxair_.airw.cbegin(), ::conex::cxair_.airw.cend(),
                           ::conex::cxair_.airz.cbegin(), 0.);

    // this is the CONEX default but actually unused there
    ::conex::cxair_.airi = {82.0e-09, 95.0e-09, 188.e-09};

    int randomSeeds[3] = {1234, 0, 0}; // SEEDS ARE NOT USED. All random numbers are
                                       // obtained from the CORSIKA 8 stream "conex"
    corsika::connect_random_stream("conex", ::conex::set_rng_function);
    int heModel = eSibyll23;

    int nShower = 1; // large to avoid final stats.
    int maxDetail = 0;
#ifdef CONEX_EXTENSIONS
    int particleListMode = 0;
#endif

    std::string configPath = CONEX_CONFIG_PATH;
    ::conex::initconex_(nShower, randomSeeds, heModel, maxDetail,
#ifdef CONEX_EXTENSIONS
                        particleListMode,
#endif
                        configPath.c_str(), configPath.size());
  }

  template <typename TOutputE, typename TOutputN>
  inline void CONEXhybrid<TOutputE, TOutputN>::initCascadeEquations() {

    // set phi, theta
    Vector<length_d> ez{conexObservationCS_, {0._m, 0._m, -1_m}};
    auto const c = showerAxis_.getDirection().dot(ez) / 1_m;
    double theta = std::acos(c) * 180 / M_PI;

    auto const showerAxisConex =
        showerAxis_.getDirection().getComponents(conexObservationCS_);
    double phi = std::atan2(-showerAxisConex.getY().magnitude(),
                            showerAxisConex.getX().magnitude()) *
                 180 / M_PI;

    CORSIKA_LOG_DEBUG(
        "theta (deg) = {}"
        "; phi (deg) = {}",
        theta, phi);

    int ipart = static_cast<int>(primaryPDG_);

    double dimpact = 0.; // valid only if shower core is fixed on the observation plane;
                         // for skimming showers an offset is needed like in CONEX

    // SEEDS ARE NOT USED. All random numbers are obtained from
    // the CORSIKA 8 stream "conex" and "epos"!
    std::array<int, 3> ioseed{1, 1, 1};

    double eprima = primaryEnergy_ / 1_GeV;
    double xminp = injectionHeight_ / 1_m;

    ::conex::conexrun_(ipart, eprima, theta, phi, xminp, dimpact, ioseed.data());
  }

  template <typename TOutputE, typename TOutputN>
  template <typename TStackView>
  inline void CONEXhybrid<TOutputE, TOutputN>::doSecondaries(TStackView& vS) {
    auto p = vS.begin();
    while (p != vS.end()) {
      Code const pid = p.getPID();
      if (addParticle(pid, p.getEnergy(), p.getMass(), p.getPosition(),
                      p.getMomentum().normalized(), p.getTime())) {
        p.erase();
      }
      ++p;
    }
  }

  template <typename TOutputE, typename TOutputN>
  inline bool CONEXhybrid<TOutputE, TOutputN>::addParticle(
      Code pid, HEPEnergyType energy, HEPEnergyType mass, Point const& position,
      DirectionVector const& direction, TimeType t, double weight) {

    auto const it = std::find_if(egs_em_codes_.cbegin(), egs_em_codes_.cend(),
                                 [=](auto const& p) { return pid == p.first; });
    if (it == egs_em_codes_.cend()) { return false; }

    // EM particle
    auto const egs_pid = it->second;
    CORSIKA_LOG_DEBUG("position conexObs: {}",
                      position.getCoordinates(conexObservationCS_));

    auto const coords = position.getCoordinates(conexObservationCS_) / 1_m;
    double const x = coords[0].magnitude();
    double const y = coords[1].magnitude();

    double const altitude = ((position - center_).getNorm() - conex::earthRadius) / 1_m;
    auto const d = position - showerCore_;

    // distance from core to particle projected along shower axis
    double const slantDistance = -d.dot(showerAxis_.getDirection()) / 1_m;

    // lateral coordinates in CONEX shower frame
    auto const dShowerPlane = d - d.getParallelProjectionOnto(showerAxis_.getDirection());
    double const lateralX = dShowerPlane.dot(x_sf_) / 1_m;
    double const lateralY = dShowerPlane.dot(y_sf_) / 1_m;

    double const slantX = showerAxis_.getProjectedX(position) * (1_cm * 1_cm / 1_g);

    double const time = (t * constants::c - groundDist_) / 1_m;

    // fill u,v,w momentum direction in EGS frame
    double const u = direction.dot(y_sf_).magnitude();
    double const v = direction.dot(x_sf_).magnitude();
    double const w = direction.dot(showerAxis_.getDirection()).magnitude();

    // generation, TO BE CHANGED WHEN WE HAVE THAT INFORMATION AVAILABLE
    int const latchin = 1;

    double const E = energy / 1_GeV;
    double const m = mass / 1_GeV;

    CORSIKA_LOG_DEBUG("CONEXhybrid: removing {} {:5e} GeV", egs_pid, energy);

    CORSIKA_LOG_DEBUG("#### parameters to cegs4_() ####");
    CORSIKA_LOG_DEBUG("egs_pid = {}", egs_pid);
    CORSIKA_LOG_DEBUG("E = {}", E);
    CORSIKA_LOG_DEBUG("m = {}", m);
    CORSIKA_LOG_DEBUG("x = {}", x);
    CORSIKA_LOG_DEBUG("y = {}", y);
    CORSIKA_LOG_DEBUG("altitude = {}", altitude);
    CORSIKA_LOG_DEBUG("slantDistance = {}", slantDistance);
    CORSIKA_LOG_DEBUG("lateralX = {}", lateralX);
    CORSIKA_LOG_DEBUG("lateralY = {}", lateralY);
    CORSIKA_LOG_DEBUG("slantX = {}", slantX);
    CORSIKA_LOG_DEBUG("time = {}", time);
    CORSIKA_LOG_DEBUG("u = {}", u);
    CORSIKA_LOG_DEBUG("v = {}", v);
    CORSIKA_LOG_DEBUG("w = {}", w);

    ::conex::cxoptl_.dptl[10 - 1] = egs_pid;
    ::conex::cxoptl_.dptl[4 - 1] = E;
    ::conex::cxoptl_.dptl[5 - 1] = m;
    ::conex::cxoptl_.dptl[6 - 1] = x;
    ::conex::cxoptl_.dptl[7 - 1] = y;
    ::conex::cxoptl_.dptl[8 - 1] = altitude;
    ::conex::cxoptl_.dptl[9 - 1] = time;
    ::conex::cxoptl_.dptl[11 - 1] = weight;
    ::conex::cxoptl_.dptl[12 - 1] = latchin;
    ::conex::cxoptl_.dptl[13 - 1] = slantX;
    ::conex::cxoptl_.dptl[14 - 1] = lateralX;
    ::conex::cxoptl_.dptl[15 - 1] = lateralY;
    ::conex::cxoptl_.dptl[16 - 1] = slantDistance;
    ::conex::cxoptl_.dptl[2 - 1] = u;
    ::conex::cxoptl_.dptl[1 - 1] = v;
    ::conex::cxoptl_.dptl[3 - 1] = w;

    int n = 1, i = 1;
    ::conex::cegs4_(n, i);

    return true;
  }

  template <typename TOutputE, typename TOutputN>
  template <typename TStack>
  inline void CONEXhybrid<TOutputE, TOutputN>::doCascadeEquations(TStack&) {

    ::conex::conexcascade_();

    int nX = ::conex::get_number_of_depth_bins_(); // make sure this works!

    int icut = 1;
    int icutg = 2;
    int icute = 3;
    int icutm = 2;
    int icuth = 3;
    int iSec = 0;

    const int maxX = nX;

    auto X = std::make_unique<float[]>(maxX);
    auto H = std::make_unique<float[]>(maxX);
    auto D = std::make_unique<float[]>(maxX);
    auto N = std::make_unique<float[]>(maxX);
    auto dEdX = std::make_unique<float[]>(maxX);
    auto Mu = std::make_unique<float[]>(maxX);
    auto dMu = std::make_unique<float[]>(maxX);
    auto Photon = std::make_unique<float[]>(maxX);
    auto Electrons = std::make_unique<float[]>(maxX);
    auto Hadrons = std::make_unique<float[]>(maxX);

    float EGround[3], fitpars[13];

    ::conex::get_shower_data_(icut, iSec, nX, X[0], N[0], fitpars[0], H[0], D[0]);
    ::conex::get_shower_edep_(icut, nX, dEdX[0], EGround[0]);
    ::conex::get_shower_muon_(icutm, nX, Mu[0], dMu[0]);
    ::conex::get_shower_gamma_(icutg, nX, Photon[0]);
    ::conex::get_shower_electron_(icute, nX, Electrons[0]);
    ::conex::get_shower_hadron_(icuth, nX, Hadrons[0]);

    // make sure CONEX binning is same to C8:
    GrammageType const dX = (X[1] - X[0]) * (1_g / square(1_cm));

    for (int i = 0; i < nX; ++i) {
      GrammageType const curX = X[i] * (1_g / square(1_cm));
      SubWriter<TOutputE>::write(curX, curX + dX,
                                 Code::Unknown, // this is sum of all dEdX
                                 dEdX[i] * 1_GeV / 1_g * square(1_cm) * dX);
      SubWriter<TOutputN>::write(curX, curX + dX, Code::Photon, Photon[i]);
      SubWriter<TOutputN>::write(curX, curX + dX, Code::Proton, Hadrons[i]);
      SubWriter<TOutputN>::write(curX, curX + dX, Code::Electron, Electrons[i]);
      SubWriter<TOutputN>::write(curX, curX + dX, Code::MuMinus, Mu[i]);
    }
  }

  template <typename TOutputE, typename TOutputN>
  inline YAML::Node CONEXhybrid<TOutputE, TOutputN>::getConfig() const {

    return YAML::Node();
  }

} // namespace corsika