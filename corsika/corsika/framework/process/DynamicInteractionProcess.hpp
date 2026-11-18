/*
 * (c) Copyright 2023 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/geometry/FourVector.hpp>
#include <corsika/framework/core/Logging.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/process/InteractionProcess.hpp>
#include <corsika/framework/process/ProcessTraits.hpp>

#include <boost/type_index.hpp>

#include <memory>

namespace corsika {

  /**
   * This class allows selecting/using different InteractionProcesses at runtime without
   * recompiling the process sequence. The implementation is based on the "type-erasure"
   * technique.
   *
   * @tparam TStack the stack type; has to match the one used by Cascade together with
   * the ProcessSequence containing the DynamicInteractionProcess.
   */
  template <typename TStack>
  class DynamicInteractionProcess
      : InteractionProcess<DynamicInteractionProcess<TStack>> {
  public:
    using stack_view_type = typename TStack::stack_view_type;

    DynamicInteractionProcess() = default;

    /**
     * Create new DynamicInteractionProcess. Calls to this instance will be forwared
     * to the process referred to by obj. The newly created DynamicInteractionProcess
     * shares ownership of the underlying process, so that you do not need to worry
     * about it going out of scope.
     */
    template <typename TInteractionProcess>
    DynamicInteractionProcess(std::shared_ptr<TInteractionProcess> obj)
        : concept_{std::make_unique<ConcreteModel<TInteractionProcess>>(std::move(obj))} {
      CORSIKA_LOG_DEBUG("creating DynamicInteractionProcess from {}",
                        boost::typeindex::type_id<TInteractionProcess>().pretty_name());

      // poor man's check while we don't have C++20 concepts
      static_assert(is_interaction_process_v<TInteractionProcess>);

      // since we only forward interaction process API, all other types of corsika
      // processes are prohibited
      static_assert(!is_continuous_process_v<TInteractionProcess>);
      static_assert(!is_decay_process_v<TInteractionProcess>);
      static_assert(!is_stack_process_v<TInteractionProcess>);
      static_assert(!is_cascade_equations_process_v<TInteractionProcess>);
      static_assert(!is_secondaries_process_v<TInteractionProcess>);
      static_assert(!is_boundary_process_v<TInteractionProcess>);
      static_assert(!is_cascade_equations_process_v<TInteractionProcess>);
    }

    /// forwards arguments to doInteraction() of wrapped instance
    void doInteraction(stack_view_type& view, Code projectileId, Code targetId,
                       FourMomentum const& projectileP4, FourMomentum const& targetP4) {
      return concept_->doInteraction(view, projectileId, targetId, projectileP4,
                                     targetP4);
    }

    /// forwards arguments to getCrossSection() of wrapped instance
    CrossSectionType getCrossSection(Code projectileId, Code targetId,
                                     FourMomentum const& projectileP4,
                                     FourMomentum const& targetP4) const {
      return concept_->getCrossSection(projectileId, targetId, projectileP4, targetP4);
    }

  private:
    struct IInteractionModel {
      virtual ~IInteractionModel() = default;
      virtual void doInteraction(stack_view_type&, Code, Code, FourMomentum const&,
                                 FourMomentum const&) = 0;
      virtual CrossSectionType getCrossSection(Code, Code, FourMomentum const&,
                                               FourMomentum const&) const = 0;
    };

    template <typename TModel>
    struct ConcreteModel final : IInteractionModel {
      ConcreteModel(std::shared_ptr<TModel> obj) noexcept
          : model_{std::move(obj)} {}

      void doInteraction(stack_view_type& view, Code projectileId, Code targetId,
                         FourMomentum const& projectileP4,
                         FourMomentum const& targetP4) override {
        return model_->doInteraction(view, projectileId, targetId, projectileP4,
                                     targetP4);
      }

      CrossSectionType getCrossSection(Code projectileId, Code targetId,
                                       FourMomentum const& projectileP4,
                                       FourMomentum const& targetP4) const override {
        return model_->getCrossSection(projectileId, targetId, projectileP4, targetP4);
      }

    private:
      std::shared_ptr<TModel> model_;
    };

    std::unique_ptr<IInteractionModel> concept_;
  };
} // namespace corsika
