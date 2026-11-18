#pragma once

#include <corsika/framework/process/InteractionProcess.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/core/Logging.hpp>
#include <corsika/framework/core/ParticleProperties.hpp>

namespace corsika::cherenkov {

  class CherenkovProcess : public corsika::InteractionProcess<CherenkovProcess> {
  private:
    double refractiveIndex_;
    double fLambdaMin;
    double fLambdaMax;
    double fBunchSize;
    mutable corsika::logging::Logger logger_;

  public:
    CherenkovProcess(double refractiveIndex,
                    double lambdaMin = 300.0,
                    double lambdaMax = 450.0,
                    double bunchSize = 1.0);

    template <typename TParticle>
    corsika::units::si::GrammageType 
    getInteractionLength(TParticle const& particle) const;

    template <typename TSecondaryView>
    void doInteraction(TSecondaryView& view) const;
  };

} // namespace corsika::cherenkov

#include <corsika/modules/cherenkov/Cherenkov.inl>
