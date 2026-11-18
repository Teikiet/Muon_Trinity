/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/Logging.hpp>

#include <corsika/media/density_and_composition/FlatExponential.hpp>
#include <corsika/media/density_and_composition/HomogeneousMedium.hpp>
#include <corsika/media/density_and_composition/FlatInterpolatedMedium.hpp>
#include <corsika/media/density_and_composition/SlidingPlanarExponential.hpp>
#include <corsika/media/density_and_composition/SlidingPlanarTabular.hpp>

namespace corsika {
  namespace media {

    template <typename TMediumInterface, template <typename> typename TMediumModelExtra,
              typename... TModelArgs>
    inline void LayeredSphericalAtmosphereBuilder<
        TMediumInterface, TMediumModelExtra, TModelArgs...>::checkRadius(LengthType const
                                                                             r) const {
      if (r <= previousRadius_) {
        throw std::runtime_error("radius must be greater than previous");
      }
    }

    template <typename TMediumInterface, template <typename> typename TMediumModelExtra,
              typename... TModelArgs>
    inline void LayeredSphericalAtmosphereBuilder<TMediumInterface, TMediumModelExtra,
                                                  TModelArgs...>::
        setNuclearComposition(media::NuclearComposition const& composition) {
      composition_ = std::make_unique<NuclearComposition>(composition);
    }

    template <typename TMediumInterface, template <typename> typename TMediumModelExtra,
              typename... TModelArgs>
    inline typename LayeredSphericalAtmosphereBuilder<TMediumInterface, TMediumModelExtra,
                                                      TModelArgs...>::volume_tree_node*
    LayeredSphericalAtmosphereBuilder<
        TMediumInterface, TMediumModelExtra,
        TModelArgs...>::addExponentialLayer(GrammageType const b,
                                            LengthType const scaleHeight,
                                            LengthType const upperBoundary) {

      // outer radius
      auto const radius = planetRadius_ + upperBoundary;
      checkRadius(radius);
      previousRadius_ = radius;

      auto node = std::make_unique<media::VolumeTreeNode<TMediumInterface>>(
          std::make_unique<Sphere>(center_, radius));

      auto const rho0 = b / scaleHeight;

      if constexpr (detail::has_extra_models<TMediumModelExtra>::value) {
        // helper lambda in which the last 5 arguments to make_shared<...> are bound
        auto lastBound = [&](auto... argPack) {
          return std::make_shared<
              TMediumModelExtra<SlidingPlanarExponential<TMediumInterface>>>(
              argPack..., center_, rho0, -scaleHeight, *composition_, planetRadius_);
        };

        // now unpack the additional arguments
        auto model = std::apply(lastBound, additionalModelArgs_);
        node->setModelProperties(std::move(model));
      } else {
        node->template setModelProperties<SlidingPlanarExponential<TMediumInterface>>(
            center_, rho0, -scaleHeight, *composition_, planetRadius_);
      }

      layers_.push(std::move(node));
      return layers_.top().get();
    }

    template <typename TMediumInterface, template <typename> typename TMediumModelExtra,
              typename... TModelArgs>
    inline void LayeredSphericalAtmosphereBuilder<
        TMediumInterface, TMediumModelExtra,
        TModelArgs...>::addLinearLayer(GrammageType const b, LengthType const scaleHeight,
                                       LengthType const upperBoundary) {
      // outer radius
      auto const radius = planetRadius_ + upperBoundary;
      checkRadius(radius);
      previousRadius_ = radius;

      auto node = std::make_unique<media::VolumeTreeNode<TMediumInterface>>(
          std::make_unique<Sphere>(center_, radius));

      auto const rho0 = b / scaleHeight;

      if constexpr (detail::has_extra_models<TMediumModelExtra>::value) {
        // helper lambda in which the last 2 arguments to make_shared<...> are bound
        auto lastBound = [&](auto... argPack) {
          return std::make_shared<
              TMediumModelExtra<media::HomogeneousMedium<TMediumInterface>>>(
              argPack..., rho0, *composition_);
        };

        // now unpack the additional arguments
        auto model = std::apply(lastBound, additionalModelArgs_);
        node->setModelProperties(std::move(model));
      } else {
        node->template setModelProperties<media::HomogeneousMedium<TMediumInterface>>(
            rho0, *composition_);
      }

      layers_.push(std::move(node));
    }

    template <typename TMediumInterface, template <typename> typename TMediumModelExtra,
              typename... TModelArgs>
    inline void LayeredSphericalAtmosphereBuilder<TMediumInterface, TMediumModelExtra,
                                                  TModelArgs...>::
        addTabularLayer(std::function<MassDensityType(LengthType)> const& funcRho,
                        unsigned int const nBins, LengthType const deltaHeight,
                        LengthType const upperBoundary) {

      auto const radius = planetRadius_ + upperBoundary;
      checkRadius(radius);
      previousRadius_ = radius;

      auto node = std::make_unique<media::VolumeTreeNode<TMediumInterface>>(
          std::make_unique<Sphere>(center_, radius));

      if constexpr (detail::has_extra_models<TMediumModelExtra>::value) {
        // helper lambda in which the last 5 arguments to make_shared<...> are bound
        auto lastBound = [&](auto... argPack) {
          return std::make_shared<
              TMediumModelExtra<SlidingPlanarTabular<TMediumInterface>>>(
              argPack..., center_, funcRho, nBins, deltaHeight, *composition_,
              planetRadius_);
        };

        // now unpack the additional arguments
        auto model = std::apply(lastBound, additionalModelArgs_);
        node->setModelProperties(std::move(model));
      } else {
        node->template setModelProperties<SlidingPlanarTabular<TMediumInterface>>(
            center_, funcRho, nBins, deltaHeight, *composition_, planetRadius_);
      }
      layers_.push(std::move(node));
    }

    template <typename TMediumInterface, template <typename> typename TMediumModelExtra,
              typename... TModelArgs>
    inline void LayeredSphericalAtmosphereBuilder<TMediumInterface, TMediumModelExtra,
                                                  TModelArgs...>::
        addInterpolatedLayer(std::vector<LengthType> const& altitudes,
                             std::vector<MassDensityType> const& densities) {

      // Determine outer radius from the last entry of altitude vector
      auto const outer_radius = planetRadius_ + altitudes.at(altitudes.size() - 1);
      checkRadius(outer_radius);
      previousRadius_ = outer_radius;

      auto node = std::make_unique<VolumeTreeNode<TMediumInterface>>(
          std::make_unique<Sphere>(center_, outer_radius));

      if constexpr (detail::has_extra_models<TMediumModelExtra>::value) {
        // helper lambda in which the last 3 arguments to make_shared<...> are bound
        auto lastBound = [&](auto... argPack) {
          return std::make_shared<
              TMediumModelExtra<FlatInterpolatedMedium<TMediumInterface>>>(
              argPack..., center_, planetRadius_, altitudes, densities, *composition_);
        };

        // now unpack the additional arguments
        auto model = std::apply(lastBound, additionalModelArgs_);
        node->setModelProperties(std::move(model));
      } else {
        node->template setModelProperties<FlatInterpolatedMedium<TMediumInterface>>(
            center_, planetRadius_, altitudes, densities, *composition_);
      }

      layers_.push(std::move(node));
    }

    template <typename TMediumInterface, template <typename> typename TMediumModelExtra,
              typename... TModelArgs>
    inline media::Environment<TMediumInterface> LayeredSphericalAtmosphereBuilder<
        TMediumInterface, TMediumModelExtra, TModelArgs...>::assemble() {
      media::Environment<TMediumInterface> env;
      assemble(env);
      return env;
    }

    template <typename TMediumInterface, template <typename> typename TMediumModelExtra,
              typename... TModelArgs>
    inline void LayeredSphericalAtmosphereBuilder<
        TMediumInterface, TMediumModelExtra,
        TModelArgs...>::assemble(media::Environment<TMediumInterface>& env) {
      auto& universe = env.getUniverse();
      auto* outmost = universe.get();

      while (!layers_.empty()) {
        auto l = std::move(layers_.top());
        auto* tmp = l.get();
        outmost->addChild(std::move(l));
        layers_.pop();
        outmost = tmp;
      }
    }

    template <typename TMediumInterface, template <typename> typename MExtraEnvirnoment>
    struct make_layered_spherical_atmosphere_builder {
      template <typename... TArgs>
      static auto create(Point const& center, LengthType const planetRadius,
                         TArgs... args) {
        return LayeredSphericalAtmosphereBuilder<TMediumInterface, MExtraEnvirnoment,
                                                 TArgs...>{std::forward<TArgs>(args)...,
                                                           center, planetRadius};
      }
    };

  } // namespace media
} // namespace corsika
