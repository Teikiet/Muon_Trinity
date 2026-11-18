/*
 * (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/process/BaseProcess.hpp>
#include <corsika/framework/process/ProcessTraits.hpp>
#include <corsika/framework/process/BoundaryCrossingProcess.hpp>
#include <corsika/framework/process/ContinuousProcess.hpp>
#include <corsika/framework/process/ContinuousProcessStepLength.hpp>
#include <corsika/framework/process/ContinuousProcessIndex.hpp>
#include <corsika/framework/process/DecayProcess.hpp>
#include <corsika/framework/process/InteractionProcess.hpp>
#include <corsika/framework/process/ProcessReturn.hpp>
#include <corsika/framework/process/SecondariesProcess.hpp>
#include <corsika/framework/process/StackProcess.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

#include <cmath>
#include <limits>
#include <type_traits>

namespace corsika {

  template <typename TCondition, typename TSequence, typename USequence, int IndexStart,
            int IndexProcess1, int IndexProcess2>
  template <typename TParticle>
  inline ProcessReturn SwitchProcessSequence<
      TCondition, TSequence, USequence, IndexStart, IndexProcess1,
      IndexProcess2>::doBoundaryCrossing(TParticle& particle,
                                         typename TParticle::node_type const& from,
                                         typename TParticle::node_type const& to) {
    if (select_(particle)) {
      if constexpr (is_boundary_process_v<process1_type> ||
                    process1_type::is_process_sequence) {

        // interface checking on TSequence
        if constexpr (is_boundary_process_v<process1_type>) {

          static_assert(
              has_method_doBoundaryCrossing_v<TSequence, ProcessReturn, TParticle&>,
              "TDerived has no method with correct signature \"ProcessReturn "
              "doBoundaryCrossing(TParticle&, VolumeNode const&, VolumeNode const&)\" "
              "required for "
              "BoundaryCrossingProcess<TDerived>. ");
        }

        return A_.doBoundaryCrossing(particle, from, to);
      }
    } else {

      if constexpr (is_boundary_process_v<process2_type> ||
                    process2_type::is_process_sequence) {

        // interface checking on USequence
        if constexpr (is_boundary_process_v<process2_type>) {

          static_assert(
              has_method_doBoundaryCrossing_v<USequence, ProcessReturn, TParticle>,
              "TDerived has no method with correct signature \"ProcessReturn "
              "doBoundaryCrossing(TParticle&, VolumeNode const&, VolumeNode const&)\" "
              "required for "
              "BoundaryCrossingProcess<TDerived>. ");
        }

        return B_.doBoundaryCrossing(particle, from, to);
      }
    }
    return ProcessReturn::Ok;
  }

  template <typename TCondition, typename TSequence, typename USequence, int IndexStart,
            int IndexProcess1, int IndexProcess2>
  template <typename TParticle>
  inline ProcessReturn SwitchProcessSequence<TCondition, TSequence, USequence, IndexStart,
                                             IndexProcess1, IndexProcess2>::
      doContinuous(Step<TParticle>& step,
                   [[maybe_unused]] ContinuousProcessIndex const idLimit) {
    if (select_(step.getParticlePre())) {
      if constexpr (process1_type::is_process_sequence) {
        return A_.doContinuous(step, idLimit);
      }
      if constexpr (is_continuous_process_v<process1_type>) {

        //        static_assert(
        //            has_method_doContinuous_v<TSequence, ProcessReturn, TParticle&,
        //            TTrack&> ||
        //                has_method_doContinuous_v<TSequence, ProcessReturn, TParticle&,
        //                                          TTrack const&> ||
        //                has_method_doContinuous_v<TSequence, ProcessReturn, TParticle
        //                const&,
        //                                          TTrack const&>,
        //            "TDerived has no method with correct signature \"ProcessReturn "
        //            "doContinuous(TParticle[const]&,TTrack[const]&,bool)\" required for
        //            " "ContinuousProcess<TDerived>. ");

        return A_.doContinuous(
            step, idLimit == ContinuousProcessIndex(
                                 static_cast<void const*>(std::addressof(A_))));
      }
    } else {
      if constexpr (process2_type::is_process_sequence) {
        return B_.doContinuous(step, idLimit);
      }
      if constexpr (is_continuous_process_v<process2_type>) {

        // interface checking on USequence
        //        static_assert(
        //            has_method_doContinuous_v<USequence, ProcessReturn, TParticle&,
        //            TTrack&> ||
        //                has_method_doContinuous_v<USequence, ProcessReturn, TParticle&,
        //                                          TTrack const&> ||
        //                has_method_doContinuous_v<USequence, ProcessReturn, TParticle
        //                const&,
        //                                          TTrack const&>,
        //            "TDerived has no method with correct signature \"ProcessReturn "
        //            "doContinuous(TParticle [const]&,TTrack[const]&,bool)\" required for
        //            " "ContinuousProcess<TDerived>. ");

        return B_.doContinuous(
            step, idLimit == ContinuousProcessIndex(
                                 static_cast<void const*>(std::addressof(B_))));
      }
    }
    return ProcessReturn::Ok;
  }

  template <typename TCondition, typename TSequence, typename USequence, int IndexStart,
            int IndexProcess1, int IndexProcess2>
  template <typename TSecondaries>
  inline void
  SwitchProcessSequence<TCondition, TSequence, USequence, IndexStart, IndexProcess1,
                        IndexProcess2>::doSecondaries(TSecondaries& vS) {
    const auto& particle = vS.parent();
    if (select_(particle)) {
      if constexpr (is_secondaries_process_v<process1_type> ||
                    process1_type::is_process_sequence) {

        // interface checking on TSequence
        static_assert(
            has_method_doSecondaries_v<TSequence, void, TSecondaries&> ||
                has_method_doSecondaries_v<TSequence, void, TSecondaries const&>,
            "TDerived has no method with correct signature \"void "
            "doSecondaries(TStackView [const]&)\" required for "
            "SecondariesProcessProcess<TDerived>. ");

        A_.doSecondaries(vS);
      }
    } else {
      if constexpr (is_secondaries_process_v<process2_type> ||
                    process2_type::is_process_sequence) {

        // interface checking on USequence
        static_assert(
            has_method_doSecondaries_v<USequence, void, TSecondaries&> ||
                has_method_doSecondaries_v<USequence, void, TSecondaries const&>,
            "TDerived has no method with correct signature \"void "
            "doSecondaries(TStackView [const]&)\" required for "
            "SecondariesProcessProcess<TDerived>. ");

        B_.doSecondaries(vS);
      }
    }
  }

  template <typename TCondition, typename TSequence, typename USequence, int IndexStart,
            int IndexProcess1, int IndexProcess2>
  template <typename TParticle, typename TTrack>
  inline ContinuousProcessStepLength
  SwitchProcessSequence<TCondition, TSequence, USequence, IndexStart, IndexProcess1,
                        IndexProcess2>::getMaxStepLength(TParticle& particle,
                                                         TTrack& vTrack) {
    if (select_(particle)) {
      if constexpr (process1_type::is_process_sequence) {
        return A_.getMaxStepLength(particle, vTrack);
      }
      if constexpr (is_continuous_process_v<process1_type>) {

        // interface checking on TSequence
        static_assert(has_method_getMaxStepLength_v<TSequence, LengthType,
                                                    TParticle const&, TTrack const&>,
                      "TDerived has no method with correct signature \"LengthType "
                      "getMaxStepLength(TParticle const&, TTrack const&)\" required for "
                      "ContinuousProcess<TDerived>. ");

        return ContinuousProcessStepLength(
            A_.getMaxStepLength(particle, vTrack),
            ContinuousProcessIndex(static_cast<void const*>(std::addressof(A_))));
      }
    } else {
      if constexpr (process2_type::is_process_sequence) {
        return B_.getMaxStepLength(particle, vTrack);
      }
      if constexpr (is_continuous_process_v<process2_type>) {

        // interface checking on USequence
        static_assert(has_method_getMaxStepLength_v<USequence, LengthType,
                                                    TParticle const&, TTrack const&>,
                      "TDerived has no method with correct signature \"LengthType "
                      "getMaxStepLength(TParticle const&, TTrack const&)\" required for "
                      "ContinuousProcess<TDerived>. ");

        return ContinuousProcessStepLength(
            B_.getMaxStepLength(particle, vTrack),
            ContinuousProcessIndex(static_cast<void const*>(std::addressof(B_))));
      }
    }

    // if no other process in the sequence implements it
    return ContinuousProcessStepLength(std::numeric_limits<double>::infinity() * meter);
  }

  template <typename TCondition, typename TSequence, typename USequence, int IndexStart,
            int IndexProcess1, int IndexProcess2>
  template <typename TParticle>
  CrossSectionType SwitchProcessSequence<TCondition, TSequence, USequence, IndexStart,
                                         IndexProcess1, IndexProcess2>::
      getCrossSection(TParticle const& projectile, [[maybe_unused]] Code const targetId,
                      [[maybe_unused]] FourMomentum const& targetP4) const {

    if (select_(projectile)) {
      if constexpr (is_interaction_process_v<process1_type>) {
        bool constexpr has_signature_cx1 =
            has_method_getCrossSection_v<TSequence,        // process object
                                         CrossSectionType, // return type
                                         Code, Code,       // parameters
                                         FourMomentum const&, FourMomentum const&>;
        if constexpr (has_signature_cx1) {

          return A_.getCrossSection(projectile.getPID(), targetId,
                                    {projectile.getEnergy(), projectile.getMomentum()},
                                    targetP4);
        } else {
          return A_.getCrossSection(projectile, projectile.getPID(),
                                    {projectile.getEnergy(), projectile.getMomentum()});
        }
      } else if constexpr (process1_type::is_process_sequence) {
        return A_.getCrossSection(projectile, targetId, targetP4);
      }

    } else {
      if constexpr (is_interaction_process_v<process2_type>) {
        bool constexpr has_signature_cx1 =
            has_method_getCrossSection_v<USequence,        // process object
                                         CrossSectionType, // return type
                                         Code, Code,       // parameters
                                         FourMomentum const&, FourMomentum const&>;
        if constexpr (has_signature_cx1) {

          return B_.getCrossSection(projectile.getPID(), targetId,
                                    {projectile.getEnergy(), projectile.getMomentum()},
                                    targetP4);
        } else {
          return B_.getCrossSection(projectile, targetId, targetP4);
        }
      } else if constexpr (process2_type::is_process_sequence) {
        return B_.getCrossSection(projectile, targetId, targetP4);
      }
    }
    return CrossSectionType::zero(); // default value
  }

  template <typename TCondition, typename TSequence, typename USequence, int IndexStart,
            int IndexProcess1, int IndexProcess2>
  template <typename TSecondaryView, typename TRNG>
  inline ProcessReturn SwitchProcessSequence<
      TCondition, TSequence, USequence, IndexStart, IndexProcess1,
      IndexProcess2>::selectInteraction(TSecondaryView& view,
                                        FourMomentum const& projectileP4,
                                        media::NuclearComposition const& composition,
                                        TRNG& rng,
                                        [[maybe_unused]] CrossSectionType const cx_select,
                                        [[maybe_unused]] CrossSectionType cx_sum) {

    if (select_(view.parent())) {
      if constexpr (process1_type::is_process_sequence) {
        // if A_ is a process sequence --> check inside
        return A_.selectInteraction(view, projectileP4, composition, rng, cx_select,
                                    cx_sum);
      } else if constexpr (is_interaction_process_v<process1_type>) {

        auto const& projectile = view.parent();
        Code const projectileId = projectile.getPID();

        // get cross section vector for all material components
        // for selected process A
        bool constexpr has_signature_cx1 =
            has_method_getCrossSection_v<TSequence,        // process object
                                         CrossSectionType, // return type
                                         Code, Code,       // parameters
                                         FourMomentum const&, FourMomentum const&>;
        bool constexpr has_signature_cx2 = // needed for PROPOSAL interface
            has_method_getCrossSectionTemplate_v<
                TSequence,                   // process object
                CrossSectionType,            // return type
                decltype(projectile) const&, // template argument
                decltype(projectile) const&, // parameters
                Code, FourMomentum const&>;

        static_assert((has_signature_cx1 || has_signature_cx2),
                      "TSequence has no method with correct signature \"CrossSectionType "
                      "getCrossSection(Code, Code, FourMomentum const&, FourMomentum "
                      "const&)\" required by "
                      "InteractionProcess<TSequence>. ");

        std::vector<CrossSectionType> weightedCrossSections;
        if constexpr (has_signature_cx1) {
          /*std::vector<CrossSectionType> const*/ weightedCrossSections =
              composition.getWeighted([=](Code const targetId) -> CrossSectionType {
                FourMomentum const targetP4(get_mass(targetId),
                                            MomentumVector(0_GeV, 0_GeV, 0_GeV));
                return A_.getCrossSection(projectileId, targetId, projectileP4, targetP4);
              });

          cx_sum +=
              std::accumulate(weightedCrossSections.cbegin(),
                              weightedCrossSections.cend(), CrossSectionType::zero());
        } else { // this is for PROPOSAL
          cx_sum += A_.template getCrossSection(projectile, projectileId, projectileP4);
        }

        if (cx_select < cx_sum) {

          if constexpr (has_signature_cx1) {

            // now also sample targetId from weighted cross sections
            Code const targetId = composition.sampleTarget(weightedCrossSections, rng);
            FourMomentum const targetP4(get_mass(targetId),
                                        MomentumVector(0_GeV, 0_GeV, 0_GeV));

            // interface checking on TProcess1
            static_assert(
                has_method_doInteract_v<TSequence,       // process object
                                        void,            // return type
                                        TSecondaryView,  // template argument
                                        TSecondaryView&, // method parameters
                                        Code, Code, FourMomentum const&,
                                        FourMomentum const&>,
                "TSequence has no method with correct signature \"void "
                "doInteraction<TSecondaryView>(TSecondaryView&, "
                "Code, Code, FourMomentum const&, FourMomentum const&)\" required for "
                "InteractionProcess<TSequence>. ");

            A_.template doInteraction(view, projectileId, targetId, projectileP4,
                                      targetP4);
          } else { // this is for PROPOSAL
            A_.template doInteraction(view, projectileId, projectileP4);
          }

          return ProcessReturn::Interacted;
        } // end collision branch A
      }

    } else { // selection: end branch A, start branch B

      if constexpr (process2_type::is_process_sequence) {
        // if B_ is a process sequence --> check inside
        return B_.selectInteraction(view, projectileP4, composition, rng, cx_select,
                                    cx_sum);
      } else if constexpr (is_interaction_process_v<process2_type>) {

        auto const& projectile = view.parent();
        Code const projectileId = projectile.getPID();

        // get cross section vector for all material components, for selected process B
        bool constexpr has_signature_cx1 =
            has_method_getCrossSection_v<USequence,        // process object
                                         CrossSectionType, // return type
                                         Code, Code,       // parameters
                                         FourMomentum const&, FourMomentum const&>;
        bool constexpr has_signature_cx2 = // needed for PROPOSAL interface
            has_method_getCrossSectionTemplate_v<
                USequence,                   // process object
                CrossSectionType,            // return type
                decltype(projectile) const&, // template argument
                decltype(projectile) const&, // parameters
                Code, FourMomentum const&>;

        static_assert((has_signature_cx1 || has_signature_cx2),
                      "USequence has no method with correct signature \"CrossSectionType "
                      "getCrossSection(Code, Code, FourMomentum const&, FourMomentum "
                      "const&)\" required by "
                      "InteractionProcess<USequence>. ");

        std::vector<CrossSectionType> weightedCrossSections;
        if constexpr (has_signature_cx1) {
          /* std::vector<CrossSectionType> const*/ weightedCrossSections =
              composition.getWeighted([=](Code const targetId) -> CrossSectionType {
                FourMomentum const targetP4(get_mass(targetId),
                                            MomentumVector(0_GeV, 0_GeV, 0_GeV));
                return B_.getCrossSection(projectileId, targetId, projectileP4, targetP4);
              });

          cx_sum +=
              std::accumulate(weightedCrossSections.begin(), weightedCrossSections.end(),
                              CrossSectionType::zero());
        } else { // this is for PROPOSAL
          cx_sum += B_.template getCrossSection(projectile, projectileId, projectileP4);
        }

        // check if we should execute THIS process and then EXIT
        if (cx_select < cx_sum) {

          if constexpr (has_signature_cx1) {

            // now also sample targetId from weighted cross sections
            Code const targetId = composition.sampleTarget(weightedCrossSections, rng);
            FourMomentum const targetP4(get_mass(targetId),
                                        MomentumVector(0_GeV, 0_GeV, 0_GeV));

            // interface checking on TProcess2
            static_assert(
                has_method_doInteract_v<USequence,       // process object
                                        void,            // return type
                                        TSecondaryView,  // template argument
                                        TSecondaryView&, // method parameters
                                        Code, Code, FourMomentum const&,
                                        FourMomentum const&>,
                "USequence has no method with correct signature \"void "
                "doInteraction<TSecondaryView>(TSecondaryView&, "
                "Code, Code, FourMomentum const&, FourMomentum const&)\" required for "
                "InteractionProcess<USequence>. ");

            B_.doInteraction(view, projectileId, targetId, projectileP4, targetP4);
          } else { // this is for PROPOSAL
            B_.doInteraction(view, projectileId, projectileP4);
          }

          return ProcessReturn::Interacted;
        } // end collision in branch B
      }
    } // end branch B_

    return ProcessReturn::Ok;
  } // namespace corsika

  template <typename TCondition, typename TSequence, typename USequence, int IndexStart,
            int IndexProcess1, int IndexProcess2>
  template <typename TParticle>
  inline InverseTimeType
  SwitchProcessSequence<TCondition, TSequence, USequence, IndexStart, IndexProcess1,
                        IndexProcess2>::getInverseLifetime(TParticle&& particle) {

    if (select_(particle)) {
      if constexpr (is_decay_process_v<process1_type> ||
                    process1_type::is_process_sequence) {
        return A_.getInverseLifetime(particle);
      }

    } else {

      if constexpr (is_decay_process_v<process2_type> ||
                    process2_type::is_process_sequence) {
        return B_.getInverseLifetime(particle);
      }
    }
    return 0 / second; // default value
  }

  template <typename TCondition, typename TSequence, typename USequence, int IndexStart,
            int IndexProcess1, int IndexProcess2>
  // select decay process
  template <typename TSecondaryView>
  inline ProcessReturn SwitchProcessSequence<
      TCondition, TSequence, USequence, IndexStart, IndexProcess1,
      IndexProcess2>::selectDecay(TSecondaryView& view,
                                  [[maybe_unused]] InverseTimeType decay_inv_select,
                                  [[maybe_unused]] InverseTimeType decay_inv_sum) {
    if (select_(view.parent())) {
      if constexpr (process1_type::is_process_sequence) {
        // if A_ is a process sequence --> check inside
        ProcessReturn const ret = A_.selectDecay(view, decay_inv_select, decay_inv_sum);
        // if A_ did succeed, stop routine here (not checking other static branch B_)
        if (ret != ProcessReturn::Ok) { return ret; }
      } else if constexpr (is_decay_process_v<process1_type>) {
        // if this is not a ContinuousProcess --> evaluate probability
        decay_inv_sum += A_.getInverseLifetime(view.parent());
        // check if we should execute THIS process and then EXIT
        if (decay_inv_select < decay_inv_sum) {
          // more pedagogical: rndm_select < decay_inv_sum / decay_inv_tot

          // interface checking on TSequence
          static_assert(has_method_doDecay_v<TSequence, void, TSecondaryView&>,
                        "TDerived has no method with correct signature \"void "
                        "doDecay(TSecondaryView&)\" required for "
                        "DecayProcess<TDerived>. ");

          A_.doDecay(view);
          return ProcessReturn::Decayed;
        }
      } // end branch A_

    } else {

      if constexpr (process2_type::is_process_sequence) {
        // if B_ is a process sequence --> check inside
        return B_.selectDecay(view, decay_inv_select, decay_inv_sum);
      } else if constexpr (is_decay_process_v<process2_type>) {
        // if this is not a ContinuousProcess --> evaluate probability
        decay_inv_sum += B_.getInverseLifetime(view.parent());
        // check if we should execute THIS process and then EXIT
        if (decay_inv_select < decay_inv_sum) {

          // interface checking on TSequence
          static_assert(has_method_doDecay_v<USequence, void, TSecondaryView&>,
                        "TDerived has no method with correct signature \"void "
                        "doDecay(TSecondaryView&)\" required for "
                        "DecayProcess<TDerived>. ");

          B_.doDecay(view);
          return ProcessReturn::Decayed;
        }
      } // end branch B_
    }
    return ProcessReturn::Ok;
  }

} // namespace corsika
