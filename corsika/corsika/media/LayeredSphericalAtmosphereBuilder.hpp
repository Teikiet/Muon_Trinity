/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/Point.hpp>
#include <corsika/media/Environment.hpp>
#include <corsika/media/interfaces/IMediumModel.hpp>
#include <corsika/media/composition/NuclearComposition.hpp>
#include <corsika/media/VolumeTreeNode.hpp>

// for detail namespace, NoExtraModelInner, NoExtraModel and traits
#include <corsika/detail/media/LayeredSphericalAtmosphereBuilder.hpp>

#include <memory>
#include <stack>
#include <tuple>
#include <type_traits>
#include <functional>

/**
 * @file LayeredSphericalAtmosphereBuilder.hpp
 */

namespace corsika {
  namespace media {
    /**
     * make_layered_spherical_atmosphere_builder.
     *
     * Helper class to create LayeredSphericalAtmosphereBuilder, the
     * extra environment models have to be passed as template-template
     * argument to make_layered_spherical_atmosphere_builder, the member
     * function `create` does then take an unspecified number of extra
     * parameters to internalize those models for all layers later
     * produced.
     */
    template <typename TMediumInterface = IMediumModel,
              template <typename> typename MExtraEnvirnoment = detail::NoExtraModel>
    struct make_layered_spherical_atmosphere_builder;

    /**
     * Helper class to setup concentric spheres of layered atmosphere
     * with spcified density profiles (exponential, linear, ...).
     *
     * This can be used most importantly to replicate CORSIKA7
     * atmospheres.
     *
     * Each layer by definition has a density profile and a (constant)
     * nuclear composition model.
     */

    template <typename TMediumInterface = IMediumModel,
              template <typename> typename TMediumModelExtra = detail::NoExtraModel,
              typename... TModelArgs>
    class LayeredSphericalAtmosphereBuilder {

      LayeredSphericalAtmosphereBuilder() = delete;
      LayeredSphericalAtmosphereBuilder(const LayeredSphericalAtmosphereBuilder&) =
          delete;
      LayeredSphericalAtmosphereBuilder(const LayeredSphericalAtmosphereBuilder&&) =
          delete;
      LayeredSphericalAtmosphereBuilder& operator=(
          const LayeredSphericalAtmosphereBuilder&) = delete;

      // friend, to allow construction
      template <typename, template <typename> typename>
      friend struct make_layered_spherical_atmosphere_builder;

    protected:
      LayeredSphericalAtmosphereBuilder(TModelArgs... args, Point const& center,
                                        LengthType const planetRadius)
          : center_(center)
          , planetRadius_(planetRadius)
          , additionalModelArgs_{args...} {}

    public:
      typedef typename VolumeTreeNode<TMediumInterface>::VTN_type volume_tree_node;
      typedef typename VolumeTreeNode<TMediumInterface>::VTNUPtr volume_tree_node_uptr;

      void setNuclearComposition(media::NuclearComposition const& composition);
      volume_tree_node* addExponentialLayer(GrammageType const b,
                                            LengthType const scaleHeight,
                                            LengthType const upperBoundary);
      void addLinearLayer(GrammageType const b, LengthType const scaleHeight,
                          LengthType const upperBoundary);

      void addTabularLayer(std::function<MassDensityType(LengthType)> const& funcRho,
                           unsigned int const nBins, LengthType const deltaHeight,
                           LengthType const upperBoundary);

      /**
       * Add a layer where parameters are interpolated
       * @param alt Vector of altitudes (measured from the Earth surface)
       * @param density Vector of densities at given altitudes
       */
      void addInterpolatedLayer(std::vector<LengthType> const& alt,
                                std::vector<MassDensityType> const& density);

      int getSize() const { return layers_.size(); }

      void assemble(media::Environment<TMediumInterface>& env);
      media::Environment<TMediumInterface> assemble();

      /**
       * Get the current planet radius.
       */
      LengthType getPlanetRadius() const { return planetRadius_; }

    private:
      void checkRadius(LengthType const r) const;

      std::unique_ptr<NuclearComposition> composition_;
      Point center_;
      LengthType previousRadius_{LengthType::zero()};
      LengthType planetRadius_;
      std::tuple<TModelArgs...> const additionalModelArgs_;

      std::stack<volume_tree_node_uptr> layers_; // innermost layer first

    }; // end class LayeredSphericalAtmosphereBuilder

  } // namespace media
} // namespace corsika

#include <corsika/detail/media/LayeredSphericalAtmosphereBuilder.inl>
