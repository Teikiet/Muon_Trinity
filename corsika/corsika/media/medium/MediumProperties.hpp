/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

namespace corsika {
  namespace media {
    /**
     * @defgroup MediaProperties Media Properties
     *
     * Medium types are useful most importantly for effective models
     * like energy losses. a particular medium (mixture of components)
     * may have specif properties not reflected by its mixture of
     * components.
     *
     * The data provided here is automatically parsed from the file
     * properties8.dat from NIST.
     *
     * The data of each known medium can be access via the global functions in namespace
     * corsika, or via a static class object with the following interface (here at the
     * example of the class HydrogenGas):
     *
     * @code
     * static constexpr Medium medium() { return media::Medium::HydrogenGas; }
     *
     * static std::string const getName() { return data_.getName(); }
     * static std::string const getPrettyName() { return data_.getPrettyName(); }
     * static double getWeight() { return data_.getWeight (); }
     * static int weight_significant_figure() { return data_.weight_significant_figure ();
     * } static int weight_error_last_digit() { return data_.weight_error_last_digit(); }
     * static double Z_over_A() { return data_.Z_over_A(); }
     * static double getSternheimerDensity() { return data_.getSternheimerDensity(); }
     * static double getCorrectedDensity() { return data_.getCorrectedDensity(); }
     * static StateOfMatter getStateOfMatter() { return data_.getStateOfMatter(); }
     * static MediumType getType() { return data_.getType(); }
     * static std::string const getSymbol() { return data_.getSymbol(); }
     *
     * static double getIeff() { return data_.getIeff(); }
     * static double getCbar() { return data_.getCbar(); }
     * static double getX0() { return data_.getX0(); }
     * static double getX1() { return data_.getX1(); }
     * static double getAA() { return data_.getAA(); }
     * static double getSK() { return data_.getSK(); }
     * static double getDlt0() { return data_.getDlt0(); }
     *
     * inline static const MediumData data_ { "hydrogen_gas", "hydrogen gas
     * (H%2#)", 1.008, 3, 7, 0.99212, 8.3748e-05, 8.3755e-05, StateOfMatter::DiatomicGas,
     * MediumType::Element, "H", 19.2, 9.5835, 1.8639, 3.2718, 0.14092, 5.7273, 0.0 };
     * @endcode
     *
     * The numeric data known to CORSIKA 8 (and obtained from NIST) can be browsed in
     * @ref MediaPropertiesClasses.
     *
     * @ingroup MediaProperties
     * @{
     */

    //! General type of medium
    enum class MediumType {
      Unknown,
      Element,
      RadioactiveElement,
      InorganicCompound,
      OrganicCompound,
      Polymer,
      Mixture,
      BiologicalDosimetry
    };

    //! Physical state of medium
    enum class StateOfMatter { Unknown, Solid, Liquid, Gas, DiatomicGas };

    /**
     * Enum for all known Media types.
     */
    enum class Medium : int16_t;

    /**
     * The internal integer type of a medium enum.
     */
    typedef std::underlying_type<Medium>::type MediumIntType;

    /**
     *
     * Simple object to group together the properties of a medium.
     */
    struct MediumData {
      std::string name_;
      std::string pretty_name_;
      double weight_;
      int weight_significant_figure_;
      int weight_error_last_digit_;
      double Z_over_A_;
      double sternheimer_density_;
      double corrected_density_;
      StateOfMatter state_;
      MediumType type_;
      std::string symbol_;
      double Ieff_;
      double Cbar_;
      double x0_;
      double x1_;
      double aa_;
      double sk_;
      double dlt0_;

      //! @name MediumDataInterface Interface methods
      //! Interface functions for MediumData
      //! @{
      std::string getName() const { return name_; }              /// returns name
      std::string getPrettyName() const { return pretty_name_; } /// returns pretty name
      double getWeight() const { return weight_; }               /// return weight
      const int& weight_significant_figure() const {
        return weight_significant_figure_;
      } /// return significnat figures of weight
      const int& weight_error_last_digit() const {
        return weight_error_last_digit_;
      }                                                    /// return error of weight
      const double& Z_over_A() const { return Z_over_A_; } /// Z_over_A_
      double getSternheimerDensity() const {
        return sternheimer_density_;
      } /// Sternheimer density
      double getCorrectedDensity() const {
        return corrected_density_;
      }                                                         /// corrected density
      StateOfMatter getStateOfMatter() const { return state_; } /// state
      MediumType getType() const { return type_; }              /// type
      std::string getSymbol() const { return symbol_; }         /// symbol
      double getIeff() const { return Ieff_; }                  /// Ieff
      double getCbar() const { return Cbar_; }                  /// Cbar
      double getX0() const { return x0_; }                      /// X0
      double getX1() const { return x1_; }                      /// X1
      double getAA() const { return aa_; }                      /// AA
      double getSK() const { return sk_; }                      /// Sk
      double getDlt0() const { return dlt0_; }                  /// Delta0
      //! @}
    };

    //! @}

  } // namespace media
} // namespace corsika

#include <corsika/media/GeneratedMediaProperties.inc>

namespace corsika {
  namespace media {

    /**
     * @file MediaProperties.hpp
     *
     * @ingroup MediaProperties
     * @{
     *
     * Returns MediumData object for medium identifed by enum Medium.
     */

    constexpr MediumData const& mediumData(Medium const m) {
      return corsika::media::detail::medium_data[static_cast<MediumIntType>(m)];
    }

    //! @}

  } // namespace media
} // namespace corsika
