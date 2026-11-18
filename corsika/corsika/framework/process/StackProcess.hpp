/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/process/BaseProcess.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

#include <corsika/detail/framework/process/StackProcess.hpp> // for extra traits, method/interface checking

namespace corsika {

  /**
   * @ingroup Processes
   * @{
   *
   * Process to act on the entire particle stack.
   *
   * Create a new StackProcess, e.g. for XYModel, via:
   * @code{.cpp}
   * class XYModel : public StackProcess<XYModel> {};
   * @endcode
   *
   * and provide the necessary interface method:
   * @code{.cpp}
   * template <typename TStack>
   * void XYModel::doStack(TStack&);
   * @endcode
   *
   * Where, of course, Stack is the valid
   * class to access particles on the Stack. This methods does
   * not need to be templated, they could use the types
   * e.g. `corsika::setup::Stack` directly -- but by the cost of
   * loosing all flexibility otherwise provided.
   *
   * A StackProcess has only one constructor `StackProcess::StackProcess(unsigned int
   * const nStep)` where nStep (@f$n_{step}@f$) is the number of steps of the cascade
   * stepping after which the stack process should be run. Good values are on the order of
   * 1000, which will not compromise run time in the end, but provide all the benefits of
   * the StackProcess.
   *
   * The number of *steps* during the cascade processing after
   * which a StackProcess is going to be executed is determined from `iStep_` and `nStep_`
   * using the modulo: @f$ !(iStep\_ \; \% \; nStep\_) @f$. And `iStep_` is increased
   * for each evaluation (step).
   */

  template <typename TDerived>
  class StackProcess : public BaseProcess<TDerived> {
  private:
  protected:
  public:
    StackProcess() = delete;
    StackProcess(const unsigned int nStep)
        : nStep_(nStep) {}

    static bool const is_stack_process = true;

    //! return the current Cascade step counter
    auto getStep() const { return iStep_; }

    //! check if current step is where StackProcess should be executed, this also
    //! increases the internal step counter implicitly
    bool checkStep() { return !((++iStep_) % nStep_); }

  private:
    /**
     * @name Execution control and counters.
     * @{
     */
    unsigned int nStep_ = 0;
    unsigned long int iStep_ = 0;
    //! @}
  };

  /**
   * ProcessTraits specialization to flag StackProcess objects.
   */
  template <typename TProcess>
  struct is_stack_process<
      TProcess,
      std::enable_if_t<std::is_base_of_v<StackProcess<typename std::decay_t<TProcess>>,
                                         typename std::decay_t<TProcess>>>>
      : std::true_type {};

  //! @}

} // namespace corsika
