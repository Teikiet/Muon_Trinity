/*
 * (c) Copyright 2021 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/media/interfaces/IRefractiveIndexModel.hpp>
#include <corsika/media/LayeredSphericalAtmosphereBuilder.hpp>
#include <corsika/framework/utility/ImplementsMixin.hpp>
#include <corsika/media/composition/NuclearComposition.hpp>

// for detail namespace, NoExtraModelInner, NoExtraModel and traits
#include <corsika/detail/media/LayeredSphericalAtmosphereBuilder.hpp>

namespace corsika {
  namespace media {

    /**
     * Atmosphere Ids following the CORSIKA 7 5-layered atmosphere models.
     *
     * Each model corresponds to a standard 5-layered atmosphere model. 4 Layers are
     * exponential, the outer layer is with constant density.
     *
     * All atmospheres are valid for heights (above Earth sea level) up to 112.8 km.
     */
    enum class AtmosphereId : uint8_t {
      LinsleyUSStd = 0,
      MiddleEuropeJan,
      MiddleEuropeFeb,
      MiddleEuropeMay,
      MiddleEuropeJun,
      MiddleEuropeAug,
      MiddleEuropeOct,
      MiddleEuropeDec,
      SouthPoleMar,
      SouthPoleJul,
      SouthPoleOct,
      SouthPoleDec,
      SouthPoleJan,
      SouthPoleAug,
      MalargueWinterI,
      MalargueWinterII,
      MalargueSpring,
      MalargueSummer,
      MalargueAutumn,
      USStdBK,
      LastAtmosphere
    };

    namespace {
      /**
       * Struct to hold parameters of one layer of a CORSIKA7 atmosphere.
       *
       * The definition of each layer is according to a BaseExponential:
       * @f[
       *   \varrho = offset/scaleHeight \cdot
       * \exp\left(-(height-altitude)/scaleHeight\right)
       * @f],
       * where @f$ altitude @f$ is the height where the atmosphere layer starts, @f$
       * offset/scaleHeight
       * @f$ is the density at this height.
       */
      struct AtmosphereLayerParameters {
        LengthType altitude;
        GrammageType offset;
        LengthType scaleHeight;
      };

      /**
       * All the 5 layers of a CORSIKA7 atmosphere.
       */
      typedef std::array<AtmosphereLayerParameters, 5> AtmosphereParameters;

      /**
       * Local, internal helper function to provide "Grammage" type.
       *
       * @param v
       * @return value v with g/cm2 units
       */
      auto constexpr grammage(double const v) { return v * 1_g / (1_cm * 1_cm); }

      std::array<
          AtmosphereParameters,
          static_cast<uint8_t>(
              media::AtmosphereId::LastAtmosphere)> constexpr atmosphereParameterList{
          {{{{4_km, grammage(1222.6562), 994186.38_cm},
             {10_km, grammage(1144.9069), 878153.55_cm},
             {40_km, grammage(1305.5948), 636143.04_cm},
             {100_km, grammage(540.1778), 772170.16_cm},
             {112.8_km, grammage(1), 1e9_cm}}},
           {{{4_km, grammage(1173.9861), 919546_cm},
             {10_km, grammage(1205.7625), 963267.92_cm},
             {40_km, grammage(1386.7807), 614315_cm},
             {100_km, grammage(555.8935), 739059.6_cm},
             {112.8_km, grammage(1), 1e9_cm}}},
           {{{4_km, grammage(1240.48), 933697_cm},
             {10_km, grammage(1117.85), 765229_cm},
             {40_km, grammage(1210.9), 636790_cm},
             {100_km, grammage(608.2128), 733793.8_cm},
             {112.8_km, grammage(1), 1e9_cm}}},
           {{{4_km, grammage(1285.2782), 1088310_cm},
             {10_km, grammage(1173.1616), 935485_cm},
             {40_km, grammage(1320.4561), 635137_cm},
             {100_km, grammage(680.6803), 727312.6_cm},
             {112.8_km, grammage(1), 1e9_cm}}},
           {{{4_km, grammage(1251.474), 1032310_cm},
             {10_km, grammage(1173.321), 925528_cm},
             {40_km, grammage(1307.826), 645330_cm},
             {100_km, grammage(763.1139), 720851.4_cm},
             {112.8_km, grammage(1), 1e9_cm}}},
           {{{4_km, grammage(1103.3362), 923077_cm},
             {10_km, grammage(1226.5761), 1109960_cm},
             {40_km, grammage(1382.6933), 630217_cm},
             {100_km, grammage(685.6073), 726901.3_cm},
             {112.8_km, grammage(1), 1e9_cm}}},
           {{{4_km, grammage(1262.7013), 1059360_cm},
             {10_km, grammage(1139.0249), 888814_cm},
             {10_km, grammage(1270.2886), 639902_cm},
             {40_km, grammage(681.4061), 727251.8_cm},
             {112.8_km, grammage(1), 1e9_cm}}},
           {{{4_km, grammage(1210.4), 970276_cm},
             {10_km, grammage(1103.8629), 820946_cm},
             {40_km, grammage(1215.3545), 639074_cm},
             {100_km, grammage(629.7611), 731776.5_cm},
             {112.8_km, grammage(1), 1e9_cm}}},
           {{{4_km, grammage(1130.74), 867358_cm},
             {10_km, grammage(1052.05), 741208_cm},
             {40_km, grammage(1137.21), 633846_cm},
             {100_km, grammage(442.512), 759850_cm},
             {112.8_km, grammage(1), 5.4303203e9_cm}}},
           {{{4_km, grammage(1183.70), 875221_cm},
             {10_km, grammage(1108.06), 753212_cm},
             {40_km, grammage(1424.02), 545846_cm},
             {100_km, grammage(207.595), 793043_cm},
             {112.8_km, grammage(1), 5.9787908e9_cm}}},
           {{{4_km, grammage(1177.19), 861745_cm},
             {10_km, grammage(1125.11), 765925_cm},
             {40_km, grammage(1304.77), 581351_cm},
             {100_km, grammage(433.823), 775155_cm},
             {112.8_km, grammage(1), 7.4095699e9_cm}}},
           {{{4_km, grammage(1139.99), 861913_cm},
             {10_km, grammage(1073.82), 744955_cm},
             {40_km, grammage(1052.96), 675928_cm},
             {100_km, grammage(492.503), 829627_cm},
             {112.8_km, grammage(1), 5.8587010e9_cm}}},
           {{{2.67_km, grammage(1133.10), 861730_cm},
             {5.33_km, grammage(1101.20), 826340_cm},
             {8_km, grammage(1085.00), 790950_cm},
             {100_km, grammage(1098.00), 682800_cm},
             {112.8_km, grammage(1), 26798156e9_cm}}},
           {{{6.67_km, grammage(1079.00), 764170_cm},
             {13.33_km, grammage(1071.90), 699910_cm},
             {20_km, grammage(1182.00), 635650_cm},
             {100_km, grammage(1647.10), 551010_cm},
             {112.8_km, grammage(1), 59.329575e9_cm}}},
           {{{8_km, grammage(1198.5972), 945766.30_cm},
             {18.1_km, grammage(1198.8796), 681780.12_cm},
             {34.5_km, grammage(1419.4152), 620224.52_cm},
             {100_km, grammage(730.6380), 728157.92_cm},
             {112.8_km, grammage(1), 1e9_cm}}},
           {{{8.3_km, grammage(1179.5010), 939228.66_cm},
             {12.9_km, grammage(1172.4883), 787969.34_cm},
             {34_km, grammage(1437.4911), 620008.53_cm},
             {100_km, grammage(761.3281), 724585.33_cm},
             {112.8_km, grammage(1), 1e9_cm}}},
           {{{5.9_km, grammage(1202.8804), 977139.52_cm},
             {12.0_km, grammage(1148.6275), 858087.01_cm},
             {34.5_km, grammage(1432.0312), 614451.60_cm},
             {100_km, grammage(696.42788), 730875.73_cm},
             {112.8_km, grammage(1), 1e9_cm}}},
           {{{9_km, grammage(1175.3347), 986169.72_cm},
             {14.6_km, grammage(1180.3694), 793171.45_cm},
             {33_km, grammage(1614.5404), 600120.95_cm},
             {100_km, grammage(755.56438), 725247.87_cm},
             {112.8_km, grammage(1), 1e9_cm}}},
           {{{8_km, grammage(1196.9290), 985241.1_cm},
             {13_km, grammage(1173.2537), 819245.00_cm},
             {33.5_km, grammage(1502.1837), 611220.86_cm},
             {100_km, grammage(750.89704), 725797.06_cm},
             {112.8_km, grammage(1), 1e9_cm}}},
           {{{7_km, grammage(1183.6071), 954248.34_cm},
             {11.4_km, grammage(1143.0425), 800005.34_cm},
             {37_km, grammage(1322.9748), 629568.93_cm},
             {100_km, grammage(655.67307), 737521.77_cm},
             {112.8_km, grammage(1), 1e9_cm}}}}};
    } // namespace

    /**
     * Function to create a CORSIKA 7 5-layer atmosphere.
     *
     * @tparam TEnvironmentInterface
     * @tparam TExtraEnv
     * @tparam TEnvironment
     * @tparam TArgs
     * @param env
     * @param atmId
     * @param center
     * @param args
     */
    template <typename TEnvironmentInterface,
              template <typename> typename TExtraEnv = detail::NoExtraModel,
              typename TEnvironment, typename... TArgs>
    void create_5layer_atmosphere(TEnvironment& env, AtmosphereId const atmId,
                                  Point const& center, TArgs... args);

    //! The standard/default air composition with fraction values based on CORSIKA 7
    //! Composition (N2,O2,Ar) = (78.084, 20.946, 0.934)
    //! Pfraction(Ar) = Ar/(2*N2 + 2*O2 + Ar) = 0.00469
    static inline NuclearComposition const standardAirComposition{
        {Code::Nitrogen, Code::Oxygen, Code::Argon}, {0.78479, .21052, 0.00469}};

  } // namespace media
} // namespace corsika

#include <corsika/detail/media/CORSIKA7Atmospheres.inl>
