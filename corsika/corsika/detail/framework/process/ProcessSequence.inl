/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>

#include <corsika/framework/process/BaseProcess.hpp>
#include <corsika/framework/process/BoundaryCrossingProcess.hpp>
#include <corsika/framework/process/ContinuousProcess.hpp>
#include <corsika/framework/process/ContinuousProcessStepLength.hpp>
#include <corsika/framework/process/ContinuousProcessIndex.hpp>
#include <corsika/framework/process/DecayProcess.hpp>
#include <corsika/framework/process/InteractionProcess.hpp>
#include <corsika/framework/process/ProcessReturn.hpp>
#include <corsika/framework/process/SecondariesProcess.hpp>
#include <corsika/framework/process/StackProcess.hpp>
#include <corsika/framework/process/CascadeEquationsProcess.hpp>
#include <corsika/framework/core/Step.hpp>

#include <cmath>
#include <limits>
#include <type_traits>

namespace corsika {

  template <typename TProcess1, typename TProcess2, int IndexStart, int IndexProcess1,
            int IndexProcess2>
  inline ProcessSequence<TProcess1, TProcess2, IndexStart, IndexProcess1,
                         IndexProcess2>::ProcessSequence(TProcess1 in_A, TProcess2 in_B)
      : A_(in_A)
      , B_(in_B) {

    // make sure only BaseProcess types TProcess1/2 are passed

    static_assert(is_process_v<TProcess1>,
                  "can only use process derived from BaseProcess in "
                  "ProcessSequence, for Process 1");
    static_assert(is_process_v<TProcess2>,
                  "can only use process derived from BaseProcess in "
                  "ProcessSequence, for Process 2");
  }

  template <typename TProcess1, typename TProcess2, int IndexStart, int IndexProcess1,
            int IndexProcess2>
  template <typename TParticle>
  inline ProcessReturn ProcessSequence<
      TProcess1, TProcess2, IndexStart, IndexProcess1,
      IndexProcess2>::doBoundaryCrossing(TParticle& particle,
                                         typename TParticle::node_type const& from,
                                         typename TParticle::node_type const& to) {

    ProcessReturn ret = ProcessReturn::Ok;

    if constexpr (is_process_v<process1_type>) { // to protect from further compiler
                                                 // errors if process1_type is invalid
      if constexpr (is_boundary_process_v<process1_type> ||
                    process1_type::is_process_sequence) {

        // interface checking on TProcess1
        static_assert(
            has_method_doBoundaryCrossing_v<TProcess1, ProcessReturn, TParticle>,
            "TDerived has no method with correct signature \"ProcessReturn "
            "doBoundaryCrossing(TParticle&, VolumeNode const&, VolumeNode const&)\" "
            "required for "
            "BoundaryCrossingProcess<TDerived>. ");

        ret |= A_.doBoundaryCrossing(particle, from, to);
      }
    }

    if constexpr (is_process_v<process2_type>) { // to protect from further compiler
                                                 // errors if process2_type is invalid
      if constexpr (is_boundary_process_v<process2_type> ||
                    process2_type::is_process_sequence) {

        // interface checking on TProcess2
        static_assert(
            has_method_doBoundaryCrossing_v<TProcess2, ProcessReturn, TParticle>,
            "TDerived has no method with correct signature \"ProcessReturn "
            "doBoundaryCrossing(TParticle&, VolumeNode const&, VolumeNode const&)\" "
            "required for "
            "BoundaryCrossingProcess<TDerived>. ");

        ret |= B_.doBoundaryCrossing(particle, from, to);
      }
    }

    return ret;
  }

  template <typename TProcess1, typename TProcess2, int IndexStart, int IndexProcess1,
            int IndexProcess2>
  template <typename TParticle>
  inline ProcessReturn
  ProcessSequence<TProcess1, TProcess2, IndexStart, IndexProcess1, IndexProcess2>::
      doContinuous(Step<TParticle>& step,
                   [[maybe_unused]] ContinuousProcessIndex const limitId) {
    ProcessReturn ret = ProcessReturn::Ok;

    if constexpr (is_process_v<process1_type>) { // to protect from further compiler
                                                 // errors if process1_type is invalid
      if constexpr (process1_type::is_process_sequence) {

        ret |= A_.doContinuous(step, limitId);
      } else if constexpr (is_continuous_process_v<process1_type>) {

        // interface checking on TProcess1
        //~ static_assert(
        //~ has_method_doContinuous_v<TProcess1, ProcessReturn, TParticle&, TTrack&> ||
        //~ has_method_doContinuous_v<TProcess1, ProcessReturn, TParticle&,
        //~ TTrack const&> ||
        //~ has_method_doContinuous_v<TProcess1, ProcessReturn, TParticle const&,
        //~ TTrack const&>,
        //~ "TDerived has no method with correct signature \"ProcessReturn "
        //~ "doContinuous(TParticle [const]&,TTrack [const]&,bool)\" required for "
        //~ "ContinuousProcess<TDerived>. ");

        ret |= A_.doContinuous(
            step, limitId == ContinuousProcessIndex(
                                 static_cast<void const*>(std::addressof(A_))));
      }
    }

    if constexpr (is_process_v<process2_type>) { // to protect from further compiler
                                                 // errors if process2_type is invalid
      if constexpr (process2_type::is_process_sequence) {
        ret |= B_.doContinuous(step, limitId);
      } else if constexpr (is_continuous_process_v<process2_type>) {

        // interface checking on TProcess2
        //~ static_assert(
        //~ has_method_doContinuous_v<TProcess2, ProcessReturn, TParticle&, TTrack&> ||
        //~ has_method_doContinuous_v<TProcess2, ProcessReturn, TParticle&,
        //~ TTrack const&> ||
        //~ has_method_doContinuous_v<TProcess2, ProcessReturn, TParticle const&,
        //~ TTrack const&>,
        //~ "TDerived has no method with correct signature \"ProcessReturn "
        //~ "doContinuous(TParticle [const]&,TTrack [const]&,bool)\" required for "
        //~ "ContinuousProcess<TDerived>. ");

        ret |= B_.doContinuous(
            step, limitId == ContinuousProcessIndex(
                                 static_cast<void const*>(std::addressof(B_))));
      }
    }

    return ret;
  }

  template <typename TProcess1, typename TProcess2, int IndexStart, int IndexProcess1,
            int IndexProcess2>
  template <typename TSecondaries>
  inline void
  ProcessSequence<TProcess1, TProcess2, IndexStart, IndexProcess1,
                  IndexProcess2>::doSecondaries([[maybe_unused]] TSecondaries& vS) {

    if constexpr (is_process_v<process1_type>) { // to protect from further compiler
                                                 // errors if process1_type is invalid
      if constexpr (is_secondaries_process_v<process1_type> ||
                    process1_type::is_process_sequence) {

        // interface checking on TProcess1
        static_assert(
            has_method_doSecondaries_v<TProcess1, void, TSecondaries&> ||
                has_method_doSecondaries_v<TProcess1, void, TSecondaries const&>,
            "TDerived has no method with correct signature \"void "
            "doSecondaries(TStackView [const]&)\" required for "
            "SecondariesProcessProcess<TDerived>. ");

        A_.doSecondaries(vS);
      }
    }

    if constexpr (is_process_v<process2_type>) { // to protect from further compiler
                                                 // errors if process2_type is invalid
      if constexpr (is_secondaries_process_v<process2_type> ||
                    process2_type::is_process_sequence) {

        // interface checking on TProcess2
        static_assert(
            has_method_doSecondaries_v<TProcess2, void, TSecondaries&> ||
                has_method_doSecondaries_v<TProcess2, void, TSecondaries const&>,
            "TDerived has no method with correct signature \"void "
            "doSecondaries(TStackView [const]&)\" required for "
            "SecondariesProcessProcess<TDerived>. ");

        B_.doSecondaries(vS);
      }
    }
  }

  template <typename TProcess1, typename TProcess2, int IndexStart, int IndexProcess1,
            int IndexProcess2>
  inline bool ProcessSequence<TProcess1, TProcess2, IndexStart, IndexProcess1,
                              IndexProcess2>::checkStep() {
    bool ret = false;
    if constexpr (is_process_v<process1_type>) { // to protect from further compiler
                                                 // errors if process1_type is invalid
      if constexpr (is_stack_process_v<process1_type> ||
                    (process1_type::is_process_sequence &&
                     !process1_type::is_switch_process_sequence)) {
        ret |= A_.checkStep();
      }
    }
    if constexpr (is_process_v<process2_type>) { // to protect from further compiler
                                                 // errors if process2_type is invalid
      if constexpr (is_stack_process_v<process2_type> ||
                    (process2_type::is_process_sequence &&
                     !process2_type::is_switch_process_sequence)) {
        ret |= B_.checkStep();
      }
    }
    return ret;
  }

  template <typename TProcess1, typename TProcess2, int IndexStart, int IndexProcess1,
            int IndexProcess2>
  template <typename TStack>
  inline void ProcessSequence<TProcess1, TProcess2, IndexStart, IndexProcess1,
                              IndexProcess2>::doStack(TStack& stack) {
    if constexpr (is_process_v<process1_type>) { // to protect from further compiler
                                                 // errors if process1_type is invalid
      if constexpr (is_stack_process_v<process1_type> ||
                    (process1_type::is_process_sequence &&
                     !process1_type::is_switch_process_sequence)) {

        // interface checking on TProcess1
        static_assert(has_method_doStack_v<TProcess1, void, TStack&> ||
                          has_method_doStack_v<TProcess1, void, TStack const&>,
                      "TDerived has no method with correct signature \"void "
                      "doStack(TStack [const]&)\" required for "
                      "StackProcess<TDerived>. ");

        if (A_.checkStep()) { A_.doStack(stack); }
      }
    }
    if constexpr (is_process_v<process2_type>) { // to protect from further compiler
                                                 // errors if process2_type is invalid
      if constexpr (is_stack_process_v<process2_type> ||
                    (process2_type::is_process_sequence &&
                     !process2_type::is_switch_process_sequence)) {

        // interface checking on TProcess1
        static_assert(has_method_doStack_v<TProcess2, void, TStack&> ||
                          has_method_doStack_v<TProcess2, void, TStack const&>,
                      "TDerived has no method with correct signature \"void "
                      "doStack(TStack [const]&)\" required for "
                      "StackProcess<TDerived>. ");

        if (B_.checkStep()) { B_.doStack(stack); }
      }
    }
  }

  template <typename TProcess1, typename TProcess2, int IndexStart, int IndexProcess1,
            int IndexProcess2>
  template <typename TParticle, typename TTrack>
  inline ContinuousProcessStepLength
  ProcessSequence<TProcess1, TProcess2, IndexStart, IndexProcess1,
                  IndexProcess2>::getMaxStepLength(TParticle&& particle,
                                                   TTrack&& vTrack) {
    // if no other process in the sequence implements it
    ContinuousProcessStepLength max_length(std::numeric_limits<double>::infinity() *
                                           meter);

    if constexpr (is_process_v<process1_type>) { // to protect from further compiler
                                                 // errors if process1_type is invalid
      if constexpr (process1_type::is_process_sequence) {
        ContinuousProcessStepLength const step = A_.getMaxStepLength(particle, vTrack);
        max_length = std::min(max_length, step);
      } else if constexpr (is_continuous_process_v<process1_type>) {

        // interface checking on TProcess1
        static_assert(has_method_getMaxStepLength_v<TProcess1, LengthType,
                                                    TParticle const&, TTrack const&>,
                      "TDerived has no method with correct signature \"LengthType "
                      "getMaxStepLength(TParticle const&, TTrack const&)\" required for "
                      "ContinuousProcess<TDerived>. ");

        ContinuousProcessStepLength const step(
            A_.getMaxStepLength(particle, vTrack),
            ContinuousProcessIndex(static_cast<void const*>(std::addressof(A_))));
        max_length = std::min(max_length, step);
      }
    }

    if constexpr (is_process_v<process2_type>) { // to protect from further compiler
                                                 // errors if process2_type is invalid
      if constexpr (process2_type::is_process_sequence) {
        ContinuousProcessStepLength const step = B_.getMaxStepLength(particle, vTrack);
        max_length = std::min(max_length, step);
      } else if constexpr (is_continuous_process_v<process2_type>) {

        // interface checking on TProcess2
        static_assert(has_method_getMaxStepLength_v<TProcess2, LengthType,
                                                    TParticle const&, TTrack const&>,
                      "TDerived has no method with correct signature \"LengthType "
                      "getMaxStepLength(TParticle const&, TTrack const&)\" required for "
                      "ContinuousProcess<TDerived>. ");

        ContinuousProcessStepLength const step(
            B_.getMaxStepLength(particle, vTrack),
            ContinuousProcessIndex(static_cast<void const*>(std::addressof(B_))));
        max_length = std::min(max_length, step);
      }
    }
    return max_length;
  }

  template <typename TProcess1, typename TProcess2, int IndexStart, int IndexProcess1,
            int IndexProcess2>
  template <typename TParticle>
  inline CrossSectionType
  ProcessSequence<TProcess1, TProcess2, IndexStart, IndexProcess1, IndexProcess2>::
      getCrossSection([[maybe_unused]] TParticle const& projectile,
                      [[maybe_unused]] Code const targetId,
                      [[maybe_unused]] FourMomentum const& targetP4) const {

    CrossSectionType tot = CrossSectionType::zero();

    if constexpr (is_process_v<process1_type>) { // to protect from further compiler
                                                 // errors if process1_type is invalid
      if constexpr (is_interaction_process_v<process1_type>) {

        bool constexpr has_signature_cx1 =
            has_method_getCrossSection_v<TProcess1,        // process object
                                         CrossSectionType, // return type
                                         Code, Code,       // parameters
                                         FourMomentum const&, FourMomentum const&>;

        if constexpr (has_signature_cx1) {
          tot += A_.getCrossSection(projectile.getPID(), targetId,
                                    {projectile.getEnergy(), projectile.getMomentum()},
                                    targetP4);
        } else { // for PROPOSAL
          tot += A_.getCrossSection(projectile, projectile.getPID(),
                                    {projectile.getEnergy(), projectile.getMomentum()});
        }

      } else if constexpr (process1_type::is_process_sequence) {
        tot += A_.getCrossSection(projectile, targetId, targetP4);
      }
    }
    if constexpr (is_process_v<process2_type>) { // to protect from further compiler
                                                 // errors if process2_type is invalid
      if constexpr (is_interaction_process_v<process2_type>) {

        bool constexpr has_signature_cx1 =
            has_method_getCrossSection_v<TProcess2,        // process object
                                         CrossSectionType, // return type
                                         Code, Code,       // parameters
                                         FourMomentum const&, FourMomentum const&>;

        if constexpr (has_signature_cx1) {
          tot += B_.getCrossSection(projectile.getPID(), targetId,
                                    {projectile.getEnergy(), projectile.getMomentum()},
                                    targetP4);
        } else { // for PROPOSAL
          tot += B_.getCrossSection(projectile, projectile.getPID(),
                                    {projectile.getEnergy(), projectile.getMomentum()});
        }

      } else if constexpr (process2_type::is_process_sequence) {
        tot += B_.getCrossSection(projectile, targetId, targetP4);
      }
    }
    return tot;
  }

  template <typename TProcess1, typename TProcess2, int IndexStart, int IndexProcess1,
            int IndexProcess2>
  template <typename TStack>
  void ProcessSequence<TProcess1, TProcess2, IndexStart, IndexProcess1,
                       IndexProcess2>::doCascadeEquations(TStack& stack) {

    if constexpr (is_process_v<process1_type>) { // to protect from further compiler
                                                 // errors if process1_type is invalid
      if constexpr (process1_type::is_process_sequence &&
                    !process1_type::is_switch_process_sequence) {

        A_.doCascadeEquations(stack);

      } else if constexpr (is_cascade_equations_process_v<process1_type>) {

        // interface checking on TProcess1
        static_assert(has_method_doCascadeEquations_v<TProcess1,
                                                      void,     // return type
                                                      TStack&>, // parameter
                      "TDerived has no method with correct signature \"void "
                      "doCascadeEquations(TStack&)\" required for "
                      "CascadeEquationsProcess<TDerived>. ");

        A_.doCascadeEquations(stack);
      }
    }

    if constexpr (is_process_v<process2_type>) { // to protect from further compiler
                                                 // errors if process2_type is invalid
      if constexpr (process2_type::is_process_sequence &&
                    !process2_type::is_switch_process_sequence) {

        B_.doCascadeEquations(stack);

      } else if constexpr (is_cascade_equations_process_v<process2_type>) {

        // interface checking on TProcess2
        static_assert(has_method_doCascadeEquations_v<TProcess2,
                                                      void,     // return type
                                                      TStack&>, // parameter
                      "TDerived has no method with correct signature \"void "
                      "doCascadeEquations(TStack&)\" required for "
                      "CascadeEquationsProcess<TDerived>. ");

        B_.doCascadeEquations(stack);
      }
    }
  }

  template <typename TProcess1, typename TProcess2, int IndexStart, int IndexProcess1,
            int IndexProcess2>
  void ProcessSequence<TProcess1, TProcess2, IndexStart, IndexProcess1,
                       IndexProcess2>::initCascadeEquations() {

    if constexpr (is_process_v<process1_type>) { // to protect from further compiler
                                                 // errors if process1_type is invalid
      if constexpr ((process1_type::is_process_sequence &&
                     !process1_type::is_switch_process_sequence) ||
                    is_cascade_equations_process_v<process1_type>) {
        A_.initCascadeEquations();
      }
    }

    if constexpr (is_process_v<process2_type>) { // to protect from further compiler
                                                 // errors if process2_type is invalid
      if constexpr ((process2_type::is_process_sequence &&
                     !process2_type::is_switch_process_sequence) ||
                    is_cascade_equations_process_v<process2_type>) {
        B_.initCascadeEquations();
      }
    }
  } // namespace corsika

  template <typename TProcess1, typename TProcess2, int IndexStart, int IndexProcess1,
            int IndexProcess2>
  template <typename TSecondaryView, typename TRNG>
  inline ProcessReturn
  ProcessSequence<TProcess1, TProcess2, IndexStart, IndexProcess1, IndexProcess2>::
      selectInteraction(TSecondaryView&& view, FourMomentum const& projectileP4,
                        [[maybe_unused]] media::NuclearComposition const& composition,
                        [[maybe_unused]] TRNG&& rng,
                        [[maybe_unused]] CrossSectionType const cx_select,
                        [[maybe_unused]] CrossSectionType cx_sum) {

    // TODO: add check for cx_select > cx_tot

    if constexpr (is_process_v<process1_type>) { // to protect from further compiler
                                                 // errors if process1_type is invalid
      if constexpr (process1_type::is_process_sequence) {
        // if A is a process sequence --> check inside
        ProcessReturn const ret =
            A_.selectInteraction(view, projectileP4, composition, rng, cx_select, cx_sum);
        // if A_ did succeed, stop routine. Not checking other static branch B_.
        if (ret != ProcessReturn::Ok) { return ret; }
      } else if constexpr (is_interaction_process_v<process1_type>) {

        auto const& projectile = view.parent();
        Code const projectileId = projectile.getPID();

        // get cross section vector for all material components
        // for selected process A

        bool constexpr has_signature_cx1 =
            has_method_getCrossSection_v<TProcess1,        // process object
                                         CrossSectionType, // return type
                                         Code, Code,       // parameters
                                         FourMomentum const&, FourMomentum const&>;
        bool constexpr has_signature_cx2 = // needed for PROPOSAL interface
            has_method_getCrossSectionTemplate_v<
                TProcess1,                   // process object
                CrossSectionType,            // return type
                decltype(projectile) const&, // template argument
                decltype(projectile) const&, // parameters
                Code, FourMomentum const&>;

        static_assert((has_signature_cx1 || has_signature_cx2),
                      "TProcess1 has no method with correct signature \"CrossSectionType "
                      "getCrossSection(Code, Code, FourMomentum const&, FourMomentum "
                      "const&)\" required by "
                      "InteractionProcess<TProcess1>. ");

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

        // check if we should execute THIS process and then EXIT
        if (cx_select < cx_sum) {

          if constexpr (has_signature_cx1) {
            // now also sample targetId from weighted cross sections
            Code const targetId = composition.sampleTarget(weightedCrossSections, rng);
            FourMomentum const targetP4(get_mass(targetId),
                                        MomentumVector(0_GeV, 0_GeV, 0_GeV));

            // interface checking on TProcess1
            static_assert(
                has_method_doInteract_v<TProcess1,       // process object
                                        void,            // return type
                                        TSecondaryView,  // template argument
                                        TSecondaryView&, // method parameters
                                        Code, Code, FourMomentum const&,
                                        FourMomentum const&>,
                "TProcess1 has no method with correct signature \"void "
                "doInteraction<TSecondaryView>(TSecondaryView&, "
                "Code, Code, FourMomentum const&, FourMomentum const&)\" required for "
                "InteractionProcess<TProcess1>. ");

            A_.template doInteraction(view, projectileId, targetId, projectileP4,
                                      targetP4);

          } else { // this is for PROPOSAL
            A_.template doInteraction(view, projectileId, projectileP4);
          }

          return ProcessReturn::Interacted;
        }
      }
    } // end branch A

    if constexpr (is_process_v<process2_type>) { // to protect from further compiler
                                                 // errors if process2_type is invalid

      if constexpr (process2_type::is_process_sequence) {
        // if B_ is a process sequence --> check inside
        return B_.selectInteraction(view, projectileP4, composition, rng, cx_select,
                                    cx_sum);
      } else if constexpr (is_interaction_process_v<process2_type>) {

        auto const& projectile = view.parent();
        Code const projectileId = projectile.getPID();

        // get cross section vector for all material components, for selected process B
        bool constexpr has_signature_cx1 =
            has_method_getCrossSection_v<TProcess2,        // process object
                                         CrossSectionType, // return type
                                         Code, Code,       // parameters
                                         FourMomentum const&, FourMomentum const&>;
        bool constexpr has_signature_cx2 = // needed for PROPOSAL interface
            has_method_getCrossSectionTemplate_v<
                TProcess2,                    // process object
                CrossSectionType,             // return type
                decltype(*projectile) const&, // template argument
                decltype(*projectile) const&, // parameters
                Code,                         // parameters
                FourMomentum const&>;
        static_assert((has_signature_cx1 || has_signature_cx2),
                      "TProcess2 has no method with correct signature \"CrossSectionType "
                      "getCrossSection(Code, Code, FourMomentum const&, FourMomentum "
                      "const&)\" required by "
                      "InteractionProcess<TProcess1>. ");

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
                has_method_doInteract_v<TProcess2,       // process object
                                        void,            // return type
                                        TSecondaryView,  // template argument
                                        TSecondaryView&, // method parameters
                                        Code, Code, FourMomentum const&,
                                        FourMomentum const&>,
                "TProcess1 has no method with correct signature \"void "
                "doInteraction<TSecondaryView>(TSecondaryView&, "
                "Code, Code, FourMomentum const&, FourMomentum const&)\" required for "
                "InteractionProcess<TProcess2>. ");

            B_.doInteraction(view, projectileId, targetId, projectileP4, targetP4);
          } else { // this is for PROPOSAL
            B_.doInteraction(view, projectileId, projectileP4);
          }
          return ProcessReturn::Interacted;
        }
      }
    } // end branch B_
    return ProcessReturn::Ok;
  }

  template <typename TProcess1, typename TProcess2, int IndexStart, int IndexProcess1,
            int IndexProcess2>
  template <typename TParticle>
  inline InverseTimeType
  ProcessSequence<TProcess1, TProcess2, IndexStart, IndexProcess1,
                  IndexProcess2>::getInverseLifetime(TParticle&& particle) {

    InverseTimeType tot = 0 / second; // default value

    if constexpr (is_process_v<process1_type>) { // to protect from further compiler
                                                 // errors if process1_type is invalid
      if constexpr (is_decay_process_v<process1_type> ||
                    process1_type::is_process_sequence) {
        tot += A_.getInverseLifetime(particle);
      }
    }
    if constexpr (is_process_v<process2_type>) { // to protect from further compiler
                                                 // errors if process2_type is invalid
      if constexpr (is_decay_process_v<process2_type> ||
                    process2_type::is_process_sequence) {
        tot += B_.getInverseLifetime(particle);
      }
    }
    return tot;
  }

  template <typename TProcess1, typename TProcess2, int IndexStart, int IndexProcess1,
            int IndexProcess2>
  // select decay process
  template <typename TSecondaryView>
  inline ProcessReturn ProcessSequence<
      TProcess1, TProcess2, IndexStart, IndexProcess1,
      IndexProcess2>::selectDecay(TSecondaryView&& view,
                                  [[maybe_unused]] InverseTimeType decay_inv_select,
                                  [[maybe_unused]] InverseTimeType decay_inv_sum) {

    // TODO: add check for decay_inv_select>decay_inv_tot

    if constexpr (is_process_v<process1_type>) { // to protect from further compiler
                                                 // errors if process1_type is invalid
      if constexpr (process1_type::is_process_sequence) {
        // if A_ is a process sequence --> check inside
        ProcessReturn const ret = A_.selectDecay(view, decay_inv_select, decay_inv_sum);
        // if A_ did succeed, stop routine here (not checking other static branch B_)
        if (ret != ProcessReturn::Ok) { return ret; }
      } else if constexpr (is_decay_process_v<process1_type>) {
        // if this is not a ContinuousProcess --> evaluate probability
        decay_inv_sum += A_.getInverseLifetime(view.parent());
        // check if we should execute THIS process and then EXIT
        if (decay_inv_select < decay_inv_sum) { // more pedagogical: rndm_select <
                                                // decay_inv_sum / decay_inv_tot
          // interface checking on TProcess1
          static_assert(has_method_doDecay_v<TProcess1, void, TSecondaryView&>,
                        "TDerived has no method with correct signature \"void "
                        "doDecay(TSecondaryView&)\" required for "
                        "DecayProcess<TDerived>. ");

          A_.doDecay(view);
          return ProcessReturn::Decayed;
        }
      } // end branch A_
    }

    if constexpr (is_process_v<process2_type>) { // to protect from further compiler
                                                 // errors if process2_type is invalid
      if constexpr (process2_type::is_process_sequence) {
        // if B_ is a process sequence --> check inside
        return B_.selectDecay(view, decay_inv_select, decay_inv_sum);
      } else if constexpr (is_decay_process_v<process2_type>) {
        // if this is not a ContinuousProcess --> evaluate probability
        decay_inv_sum += B_.getInverseLifetime(view.parent());
        // check if we should execute THIS process and then EXIT
        if (decay_inv_select < decay_inv_sum) {

          // interface checking on TProcess1
          static_assert(has_method_doDecay_v<TProcess2, void, TSecondaryView&>,
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

  /**
   * traits marker to identify objects containing any StackProcesses.
   */
  namespace detail {
    // need helper alias to achieve this:
    template <
        typename TProcess1, typename TProcess2,
        typename = typename std::enable_if_t<
            contains_stack_process_v<TProcess1> || is_stack_process_v<TProcess1> ||
                contains_stack_process_v<TProcess2> || is_stack_process_v<TProcess2>,
            int>>
    using enable_if_stack = ProcessSequence<TProcess1, TProcess2>;
  } // namespace detail

  template <typename TProcess1, typename TProcess2>
  struct contains_stack_process<detail::enable_if_stack<TProcess1, TProcess2>>
      : std::true_type {};

} // namespace corsika
