/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

#include <functional>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

namespace corsika {
  namespace media {

    /**
     * Describes the composition of matter
     * Allowes and handles the creation of custom matter compositions.
     */

    class NuclearComposition {
    public:
      /**
       * Constructor
       *  The constructore takes a list of elements and a list which describe the relative
       *  amount. Booth lists need to have the same length and the sum all of fractions
       *  should be 1. Otherwise an exception is thrown.
       *  Example for air: Composition (N2,O2,Ar) = (78.084, 20.946, 0.934)
       *  Pfraction(Ar) = Ar/(2*N2 + 2*O2 + Ar) = 0.00469
       *  {Code::Nitrogen, Code::Oxygen, Code::Argon}, {0.78479, 0.21052, 0.00469}}
       *
       *  @param pComponents List of particle types.
       *  @param pFractions List of fractions how much each particle contributes. The sum
       *         needs to add up to 1.
       */
      NuclearComposition(std::vector<Code> const& pComponents,
                         std::vector<double> const& pFractions);

      /**
       * Returns a vector of the same length as elements in the material with the weighted
       * return of "func". The typical default application is for cross section weighted
       * with fraction in the material.
       *
       *  @tparam TFunction Type of functions for the weights. The type should be
       *          Code -> CrossSectionType.
       *  @param func Functions for reweighting specific elements.
       *  @retval returns the vector with weighted return types of func.
       */
      template <typename TFunction>
      auto getWeighted(TFunction func) const;

      /**
       * Sum all all relative composition weighted by func(element)
       *  This function sums all relative compositions given during this classes
       * construction. Each entry is weighted by the user defined function func given to
       * this function.
       *
       *  @tparam TFunction Type of functions for the weights. The type should be
       *          Code -> double.
       *  @param func Functions for reweighting specific elements.
       *  @retval returns the weighted sum with the type defined by the return type of
       * func.
       */
      template <typename TFunction>
      auto getWeightedSum(TFunction func) const -> decltype(func(std::declval<Code>()));

      /**
       * Number of elements in the composition array
       *  @retval returns the number of elements in the composition array.
       */
      size_t getSize() const;

      //! Returns a const reference to the fraction
      std::vector<double> const& getFractions() const;
      //! Returns a const reference to the fraction
      std::vector<Code> const& getComponents() const;
      double const getAverageMassNumber() const;

      template <class TRNG>
      Code sampleTarget(std::vector<CrossSectionType> const& sigma,
                        TRNG&& randomStream) const;

      // Note: when this class ever modifies its internal data, the hash
      // must be updated, too!
      size_t getHash() const;

      //! based on hash value
      bool operator==(media::NuclearComposition const& v) const;

    private:
      void updateHash();

      std::vector<double> const
          numberFractions_;                //!< relative fractions of number density
      std::vector<Code> const components_; //!< particle codes of consitutents

      double const avgMassNumber_;

      std::size_t hash_;
    };

  } // namespace media
} // namespace corsika

#include <corsika/detail/media/composition/NuclearComposition.inl>
