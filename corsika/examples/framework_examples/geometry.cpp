/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/geometry/Point.hpp>
#include <corsika/framework/geometry/CoordinateSystemRegistry.hpp>
#include <corsika/framework/geometry/Sphere.hpp>
#include <corsika/framework/geometry/Vector.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/core/Logging.hpp>

#include <cstdlib>
#include <typeinfo>

using namespace corsika;

int main() {
  logging::set_level(logging::level::info);

  CORSIKA_LOG_INFO("geometry_example");

  // define the root coordinate system
  CoordinateSystem const& root = get_root_CoordinateSystem();

  // another CS defined by a translation relative to the root CS
  CoordinateSystem const& cs2 = make_translation(root, {0_m, 0_m, 1_m});

  // rotations are possible, too; parameters are axis vector and angle
  CoordinateSystem const& cs3 =
      make_rotation(root, QuantityVector<length_d>{1_m, 0_m, 0_m}, 90 * degree_angle);

  // now let's define some geometrical objects:
  Point const p1(root, {0_m, 0_m, 0_m}); // the origin of the root CS
  Point const p2(cs2, {0_m, 0_m, 0_m});  // the origin of cs2

  Vector<length_d> const diff =
      p2 -
      p1; // the distance between the points, basically the translation vector given above
  auto const norm = diff.getSquaredNorm(); // squared length with the right dimension

  // print the components of the vector as given in the different CS
  CORSIKA_LOG_INFO(
      "p2-p1 components in root: {} \n"
      "p2-p1 components in cs2: {} \n"
      "p2-p1 components in cs3: {}\n"
      "p2-p1 norm^2: {} \n",
      diff.getComponents(root), diff.getComponents(cs2), diff.getComponents(cs3), norm);
  assert(norm == 1 * meter * meter);

  Sphere s(p1, 10_m); // define a sphere around a point with a radius
  CORSIKA_LOG_INFO("p1 inside s:{} ", s.contains(p2));
  assert(s.contains(p2) == 1);

  Sphere s2(p1, 3_um); // another sphere
  CORSIKA_LOG_INFO("p1 inside s2: {}", s2.contains(p2));
  assert(s2.contains(p2) == 0);

  // let's try parallel projections:
  auto const v1 = Vector<length_d>(root, {1_m, 1_m, 0_m});
  auto const v2 = Vector<length_d>(root, {1_m, 0_m, 0_m});

  auto const v3 = v1.getParallelProjectionOnto(v2);

  // cross product
  auto const cross =
      v1.cross(v2).normalized(); // normalized() returns dimensionless, normalized vectors

  // if a CS is not given as parameter for getComponents(), the components
  // in the "home" CS are returned
  CORSIKA_LOG_INFO(
      "v1: {} \n"
      "v2: {}\n "
      "parallel projection of v1 onto v2:  {} \n"
      "normalized cross product of v1 x v2 {} \n",
      v1.getComponents(), v2.getComponents(), v3.getComponents(), cross.getComponents());

  return EXIT_SUCCESS;
}
