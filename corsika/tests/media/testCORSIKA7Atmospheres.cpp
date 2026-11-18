/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/media/CORSIKA7Atmospheres.hpp>

#include <catch2/catch_all.hpp>

using namespace corsika;
using Catch::Approx;

TEST_CASE("CORSIKA7Atmospheres") {

  logging::set_level(logging::level::info);

  Point const gOrigin(0_m, 0_m, 0_m);

  media::Environment<media::IMediumModel> env;

  // build a Linsley US Standard atmosphere into `env`
  media::create_5layer_atmosphere<media::IMediumModel>(
      env, media::AtmosphereId::LinsleyUSStd, gOrigin);

  typedef
      typename media::Environment<media::IMediumModel>::BaseNodeType::VTN_type node_type;
  node_type const* universe = env.getUniverse().get();

  Point const p(constants::EarthRadius::Mean, 0_m, 0_m);
  CHECK(universe->getContainingNode(p)->getModelProperties().getMassDensity(p) / 1_g *
            static_pow<3>(1_cm) ==
        Approx(1222.6562 / 994186.38));
}
