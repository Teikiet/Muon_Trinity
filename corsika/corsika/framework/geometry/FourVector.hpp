/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/core/PhysicalGeometry.hpp>
#include <type_traits>

/**
 * @file FourVector.hpp
 * @author Ralf Ulrich
 * @brief General FourVector object.
 * @date 2021-10-16
 */

namespace corsika {

  /**
   * Description of physical four-vectors
   *
   * FourVector fully supports units, e.g. E in [GeV/c] and p in [GeV],
   * or also t in [s] and r in [m], etc.
   *
   * However, for HEP applications it is also possible to use E and p
   * both in [GeV].
   *
   * Thus, the input units of time-like and space-like coordinates
   * must either be idential (e.g. GeV) or scaled by "c" as in
   * [E/c]=[p].
   *
   * The FourVector can return its squared-norm \ref getNormSqr and its
   * norm \ref getNorm, whereas norm is sqrt(abs(norm-squared)). The
   * physical units are always calculated and returned properly.
   *
   * FourVector can also return if it is TimeLike, SpaceLike or PhotonLike.
   *
   * When a FourVector is initialized with a lvalue references,
   * e.g. as `FourVector<TimeType&, Vector<length_d>&>`, references
   * are also used as internal data types, which should lead to
   * complete disappearance of the FourVector class during
   * optimization.
   */

  template <typename TTimeType, typename TSpaceVecType>
  class FourVector {

  public:
    using space_vec_type = typename std::decay<TSpaceVecType>::type;
    using space_type = typename space_vec_type::quantity_type;
    using time_type = typename std::decay<TTimeType>::type;

    // check the types and the physical units here:
    static_assert(std::is_same<time_type, space_type>::value ||
                      std::is_same<time_type, decltype(std::declval<space_type>() /
                                                       meter * second)>::value,
                  "Units of time-like and space-like coordinates must either be idential "
                  "(e.g. GeV) or [E/c]=[p]");

    using norm_type = space_type;
    using norm_square_type =
        decltype(std::declval<norm_type>() * std::declval<norm_type>());

  public:
    // resource management
    FourVector() = default;
    FourVector(FourVector&&) = default;
    FourVector(FourVector const&) = default;
    FourVector& operator=(const FourVector&) = default;
    FourVector& operator=(FourVector&&) = default;
    ~FourVector() = default;

    FourVector(TTimeType const& eT, TSpaceVecType const& eS)
        : timeLike_(eT)
        , spaceLike_(eS) {}

    /**
     * @return timeLike_
     */
    TTimeType getTimeLikeComponent() const;

    /**
     * @return spaceLike_
     */
    TSpaceVecType& getSpaceLikeComponents();

    /**
     *
     * @return spaceLike_;
     */
    TSpaceVecType const& getSpaceLikeComponents() const;

    /**
     *
     * @return \f$p_0^2 - \vec{p}^2\f$
     */
    norm_square_type getNormSqr() const;

    /**
     *
     * @return \f$\sqrt(p_0^2 - \vec{p}^2)\f$
     */
    norm_type getNorm() const;

    /**
     * \todo FIXME: a better alternative would be to define an enumeration
     * enum { SpaceLike =-1, TimeLike, LightLike } V4R_Category;
     * and a method called  V4R_Category GetCategory() const;
     * RU: then you have to decide in the constructor which avoids "lazyness"
     **/
    ///\return if \f$|p_0|>|\vec{p}|\f$
    bool isTimelike() const;
    ///\return if \f$|p_0|<|\vec{p}|\f$
    bool isSpacelike() const;

    /**
     * \name Math operators (class members)
     * @{
     */
    FourVector& operator+=(FourVector const&);
    FourVector& operator-=(FourVector const&);
    FourVector& operator*=(double const);
    FourVector& operator/=(double const);

    /**
     * Scalar product of two FourVectors.
     *
     *  Note that the product between two 4-vectors assumes that you use
     *  the same "c" convention for both. Only the LHS vector is checked
     *  for this. You cannot mix different conventions due to
     *  unit-checking.
     */
    norm_type operator*(FourVector const& b) const;

    /** @} */

  protected:
    // the data members
    TTimeType timeLike_;
    TSpaceVecType spaceLike_;

    /**
     * \name Free math operators
     * We need to define them inline here since we do not want to
     * implement them as template functions. They are valid only for
     * the specific types as defined right here.
     *
     * Note, these are "free function" (even if they don't look as
     *  such). Thus, even if the input object uses internal references
     *  for storage, the free math operators, of course, must provide
     *  value-copies.
     * @{
     *
     */
    friend FourVector<time_type, space_vec_type> operator+(FourVector const& a,
                                                           FourVector const& b) {
      return FourVector<time_type, space_vec_type>(a.timeLike_ + b.timeLike_,
                                                   a.spaceLike_ + b.spaceLike_);
    }

    friend FourVector<time_type, space_vec_type> operator-(FourVector const& a,
                                                           FourVector const& b) {
      return FourVector<time_type, space_vec_type>(a.timeLike_ - b.timeLike_,
                                                   a.spaceLike_ - b.spaceLike_);
    }

    friend FourVector<time_type, space_vec_type> operator*(FourVector const& a,
                                                           double const b) {
      return FourVector<time_type, space_vec_type>(a.timeLike_ * b, a.spaceLike_ * b);
    }

    friend FourVector<time_type, space_vec_type> operator/(FourVector const& a,
                                                           double const b) {
      return FourVector<time_type, space_vec_type>(a.timeLike_ / b, a.spaceLike_ / b);
    }

    /**
     * streaming operator
     */
    friend std::ostream& operator<<(
        std::ostream& os, corsika::FourVector<TTimeType, TSpaceVecType> const& qv) {
      os << '(' << qv.timeLike_ << ", " << qv.spaceLike_ << ") ";
      return os;
    }
    /** @} */

  private:
    /**
     * This function is there to automatically remove the eventual
     * extra factor of "c" for the time-like quantity.
     */
    norm_square_type getTimeSquared() const;
  };

  /**
   * @typedef FourMomentum A FourVector with HEPEnergyType and MomentumVector.
   */
  typedef FourVector<HEPEnergyType, MomentumVector> FourMomentum;

} // namespace corsika

#include <corsika/detail/framework/geometry/FourVector.inl>
