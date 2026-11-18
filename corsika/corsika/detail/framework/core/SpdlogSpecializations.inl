/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <fmt/format.h>
#include <corsika/framework/core/ParticleProperties.hpp>
#include <boost/filesystem/path.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <Eigen/Dense>

//-----------------------------------
// STD
//-----------------------------------
namespace std {

  auto inline format_as(std::_Put_time<char> const& arg) {
    std::ostringstream os;
    os << arg;
    return os.str();
  }

} // namespace std

//-----------------------------------
// CORSIKA
//-----------------------------------
namespace corsika {

  // formatters for particle codes declared on ParticleProperties.hpp
  auto inline format_as(Code code) { return get_name(code); }
  auto inline format_as(PDGCode code) { return fmt::underlying(code); }

  template <typename Type>
  auto inline format_as(Type const& arg) {
    std::ostringstream os;
    os << arg;
    return os.str();
  }

} // namespace corsika

//-----------------------------------
// boost::filesystem
//-----------------------------------
namespace boost::filesystem {

  auto inline format_as(path const& fname) { return fname.string().c_str(); }

} // namespace boost::filesystem

//-----------------------------------
// phys::units
//-----------------------------------
namespace phys::units {

  template <typename Dimensions>
  auto inline format_as(phys::units::quantity<Dimensions> const& arg) {

    return io::to_string(arg);
  }

} // namespace phys::units

//----------------------------------
// Eigen
//----------------------------------
namespace Eigen {

  template <typename Scalar, int M, int N>
  auto inline format_as(Matrix<Scalar, M, N> const& arg) {

    std::ostringstream os;
    os << arg;
    return os.str();
  }

  template <typename Matrix1, typename Matrix2>
  auto inline format_as(Product<Matrix1, Matrix2> const& arg) {

    std::ostringstream os;
    os << arg;
    return os.str();
  }

} // namespace Eigen
