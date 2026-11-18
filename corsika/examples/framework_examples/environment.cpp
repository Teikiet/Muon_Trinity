/*
 * (c) Copyright 2022 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/media/Environment.hpp>
#include <corsika/media/density_and_composition/HomogeneousMedium.hpp>
#include <corsika/media/interfaces/IMediumModel.hpp>
#include <corsika/media/medium/MediumProperties.hpp>
#include <corsika/media/medium/MediumPropertyModel.hpp>

using namespace corsika;

// Example to show how to construct an environment in CORSIKA

// MyMediumInterface is constructed via "mixin inheritance"
// Here, we construct our base class with IMediumModel and IMediumPropertyModel:
// the former allow us define a density profile and nuclear composite,
// the later is used to determine energy losses parameters.
using MyMediumInterface = media::IMediumPropertyModel<media::IMediumModel>;

using MyEnvType = media::Environment<MyMediumInterface>;

int main() {
  // create environment and universe, empty so far
  MyEnvType env;

  // create a geometry object: a sphere with radius of 1000 m
  CoordinateSystem const& rootCS = env.getCoordinateSystem();
  Point const center{rootCS, 0_m, 0_m, 0_m};
  LengthType const radius = 1000_m;
  auto sphere = std::make_unique<Sphere>(center, radius);

  // create node from geometry object
  auto node =
      std::make_unique<media::VolumeTreeNode<MyMediumInterface>>(std::move(sphere));

  // set medium properties to our node, say it is water
  media::NuclearComposition const nucl_comp{{Code::Hydrogen, Code::Oxygen}, {0.11, 0.89}};
  MassDensityType const density = 1_g / (1_cm * 1_cm * 1_cm);

  // create concrete implementation of MyMediumInterface by combining parts that
  // implement the corresponding interfaces. HomogeneousMedium implements IMediumModel,
  // MediumPropertyModel implements IMediumPropertyModel.
  auto const medium = std::make_shared<
      media::MediumPropertyModel<media::HomogeneousMedium<MyMediumInterface>>>(
      media::Medium::WaterLiquid, density, nucl_comp);
  node->setModelProperties(medium);

  // put our node into universe
  // note: this has to be done after setting node model properties, since
  // std::move will make our previous defined node, which is a unique pointer
  // un-referenceable in the context
  media::VolumeTreeNode<MyMediumInterface>* const universe = env.getUniverse().get();
  universe->addChild(std::move(node));

  // example to explore the media properties of the node
  for (auto h = 0_m; h < radius; h += 100_m) {
    Point const ptest{rootCS, 0_m, 0_m, h};
    MyMediumInterface const& media_prop =
        env.getUniverse()->getContainingNode(ptest)->getModelProperties();
    MassDensityType const rho = media_prop.getMassDensity(ptest);
    media::NuclearComposition const& nuc_comp = media_prop.getNuclearComposition();
    std::vector<Code> const& nuclei = nuc_comp.getComponents();
    std::vector<double> const& fractions = nuc_comp.getFractions();
    CORSIKA_LOG_INFO(
        "radius: {:.2f} m, density: {:.2f} g/cm^3, nuclei: ({:d}, {:d}), fractions: "
        "({:.2f}, "
        "{:.2f})",
        h / 1_m, rho / (1_g / static_pow<3>(1_cm)), get_PDG(nuclei.at(0)),
        get_PDG(nuclei.at(1)), fractions.at(0), fractions.at(1));
  }

  return 0;
}
