/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/Logging.hpp>
#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/stack/Stack.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

namespace corsika {

  /**
   * CombinedParticleInterface can be used to combine the data of several StackData
   * objects.
   *
   * You may combine two StackData object, see class CombinedStackImpl
   * below, into one Stack, using a combined StackIterator (aka
   * CombinedParticleInterface) interface class.
   *
   * This allows to add specific information to a given Stack, could
   * be special information on a subset of entries
   * (e.g. NuclearStackExtension) or also (multi) thinning weights for
   * all particles.
   *
   * Many Stacks can be combined into more complex object.
   *
   * The two sub-stacks must both provide their independent
   * ParticleInterface classes.
   *
   */
  template <template <typename> typename TParticleInterfaceA,
            template <typename> class TParticleInterfaceB, typename TStackIterator>
  struct CombinedParticleInterface
      : public TParticleInterfaceB<TParticleInterfaceA<TStackIterator>> {

    typedef CombinedParticleInterface<TParticleInterfaceA, TParticleInterfaceB,
                                      TStackIterator>
        pi_c_type;
    typedef TParticleInterfaceA<TStackIterator> pi_a_type;
    typedef TParticleInterfaceB<TParticleInterfaceA<TStackIterator>> pi_b_type;

  protected:
    using pi_b_type::getIndex;     // choose B, A would also work
    using pi_b_type::getStackData; // choose B, A would also work

  public:
    /**
     * @name wrapper for user functions
     * @{
     *
     * In this set of functions we call the user-provide
     * TParticleInterface setParticleData(...) methods, either with
     * parent particle reference, or w/o.
     *
     * There is one implicit assumption here: if only one data tuple
     * is provided for setParticleData, the data is passed on to
     * TParticleInterfaceA and the TParticleInterfaceB is
     * default-initialized. There are many occasions where this is the
     * desired behaviour, e.g. for thinning etc.
     *
     */

    template <typename... TArgs1>
    void setParticleData(std::tuple<TArgs1...> const vA);

    template <typename... TArgs1, typename... TArgs2>
    void setParticleData(std::tuple<TArgs1...> const vA, std::tuple<TArgs2...> const vB);

    template <typename... TArgs1>
    void setParticleData(pi_a_type& p, std::tuple<TArgs1...> const vA);

    template <typename... TArgs1, typename... TArgs2>
    void setParticleData(pi_c_type& p, std::tuple<TArgs1...> const vA,
                         std::tuple<TArgs2...> const vB);
    ///@}

    std::string asString() const;

  protected:
  };

  /**
   * @class CombinedStackImpl
   *
   * Memory implementation of a combined data stack.
   *
   * The two stack data user objects Stack1Impl and Stack2Impl are
   * merged into one consistent Stack container object providing
   * access to the combined number of data entries.
   */
  template <typename Stack1Impl, typename Stack2Impl>
  struct CombinedStackImpl : public Stack1Impl, public Stack2Impl {

  public:
    void clear();

    unsigned int getSize() const { return Stack1Impl::getSize(); }
    unsigned int getCapacity() const { return Stack1Impl::getCapacity(); }

    /**
     *   Function to copy particle at location i1 in stack to i2
     */
    void copy(const unsigned int i1, const unsigned int i2);

    /**
     *   Function to copy particle at location i2 in stack to i1
     */
    void swap(const unsigned int i1, const unsigned int i2);

    void incrementSize();

    void decrementSize();

  }; // end class CombinedStackImpl

  /**
   * Helper template alias `CombinedStack` to construct new combined
   * stack from two stack data objects and a particle readout interface.
   *
   * Note that the Stack2Impl provides only /additional/ data to
   * Stack1Impl. This is important (see above) since tuple data for
   * initialization are forwarded to Stack1Impl (first).
   */

  template <typename Stack1Impl, typename Stack2Impl, template <typename> typename _PI,
            template <typename T1, template <class> class T2> class _MSecondaryProducer =
                DefaultSecondaryProducer>
  using CombinedStack =
      Stack<CombinedStackImpl<Stack1Impl, Stack2Impl>, _PI, _MSecondaryProducer>;

} // namespace corsika

#include <corsika/detail/framework/stack/CombinedStack.inl>
