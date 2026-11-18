/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/stack/ParticleBase.hpp>

#include <boost/type_index.hpp>

namespace corsika::history {
  template <typename T, template <typename> typename TParticleInterface>
  class HistorySecondaryProducer; // forward decl.
}

namespace corsika {

  template <typename TStackData, template <typename> typename TParticleInterface,
            template <class T1, template <class> class T2> class MSecondaryProducer>
  class Stack; // forward decl

  template <typename TStackData, template <typename> typename TParticleInterface,
            template <class T1, template <class> class T2> class MSecondaryProducer>
  class SecondaryView; // forward decl

  template <typename TStackData, template <typename> typename TParticleInterface,
            template <class T1, template <class> class T2> class MSecondaryProducer,
            typename StackType>
  class ConstStackIteratorInterface; // forward decl

  /**
   * The StackIteratorInterface is the main interface to iterator over
   * particles on a stack.
   *
   * At the same time StackIteratorInterface is a
   * Particle object by itself, thus there is no difference between
   * type and ref_type for convenience of the physicist.
   *
   * This allows to write code like:
   * @code
   * for (auto& p : theStack) {
   *    p.setEnergy(newEnergy);
   * }
   * @endcode
   *
   * The template argument Stack determines the type of Stack object
   * the data is stored in. A pointer to the Stack object is part of
   * the StackIteratorInterface. In addition to Stack the iterator only knows
   * the index index_ in the Stack data.
   *
   * The template argument `TParticleInterface` acts as a policy to provide
   * readout function of Particle data from the stack. The TParticleInterface
   * class must know how to retrieve information from the Stack data
   * for a particle entry at any index index_.
   *
   * The TParticleInterface class must be written and provided by the
   * user, it contains methods like <code> auto getData() const {
   * return getStackData().getData(getIndex()); }</code>, where
   * StackIteratorInterface::getStackData() return a reference to the
   * object storing the particle data of type TStackData. And
   * StackIteratorInterface::getIndex() provides the iterator index to
   * be readout. The TStackData is another user-provided class to
   * store data and must implement functions compatible with
   * TParticleInterface, in this example TStackData::getData(const unsigned int
   * vIndex).
   *
   * For two examples see stack_example.cc, or the
   * corsika::sibyll::SibStack class.
   */

  template <typename TStackData, template <typename> typename TParticleInterface,
            template <typename T1, template <class> class T2> class MSecondaryProducer,
            typename TStackType =
                Stack<TStackData, TParticleInterface, MSecondaryProducer>>
  class StackIteratorInterface
      : public TParticleInterface<StackIteratorInterface<
            TStackData, TParticleInterface, MSecondaryProducer, TStackType>> {

  public:
    typedef TParticleInterface<corsika::StackIteratorInterface<
        TStackData, TParticleInterface, MSecondaryProducer, TStackType>>
        particle_interface_type;

    typedef TStackType stack_type;
    typedef TStackData stack_data_type;

    // it is not allowed to create a "dangling" stack iterator
    StackIteratorInterface() = delete; //! \todo check rule of five

    StackIteratorInterface(StackIteratorInterface&& rhs)
        : index_(std::move(rhs.index_))
        , data_(std::move(rhs.data_)) {}

    StackIteratorInterface(StackIteratorInterface const& vR)
        : index_(vR.index_)
        , data_(vR.data_) {}

    StackIteratorInterface& operator=(StackIteratorInterface const& vR) {
      if (&vR != this) {
        index_ = vR.index_;
        data_ = vR.data_;
      }
      return *this;
    }

    /**
     *  Iterator must always point to data, with an index.
     *
     *    @param data reference to the stack [rw]
     *    @param index index on stack
     */
    StackIteratorInterface(stack_type& data, unsigned int const index)
        : index_(index)
        , data_(&data) {}

    /**
     *  Constructor that also sets new values on particle data object.
     *
     *  @param data reference to the stack [rw].
     *  @param index index on stack.
     *  @param args variadic list of data to initialize stack entry, this must be
     *        consistent with the definition of the user-provided
     *        particle_interface_type::setParticleData(...) function.
     */
    template <typename... TArgs>
    StackIteratorInterface(stack_type& data, unsigned int const index,
                           const TArgs... args)
        : index_(index)
        , data_(&data) {

      (**this).setParticleData(args...);
    }

    /**
     * Constructor that also sets new values on particle data object, including reference
     * to parent particle.
     *
     *  @param data reference to the stack [rw]
     *  @param index index on stack
     *  @param parent to parent particle [rw]. This can be used for thinning, particle
     *         counting, history, etc.
     *   @param args variadic list of data to initialize stack entry, this must be
     *        consistent with the definition of the user-provided
     *         particle_interface_type::setParticleData(...) function.
     */
    template <typename... TArgs>
    StackIteratorInterface(stack_type& data, unsigned int const index,
                           StackIteratorInterface& parent, const TArgs... args)
        : index_(index)
        , data_(&data) {

      (**this).setParticleData(*parent, args...);
    }

    bool isErased() const { return getStack().isErased(*this); }

  public:
    /**
     * @name Iterator interface
     * @{
     */
    StackIteratorInterface& operator++() {
      do {
        ++index_;
      } while (
          getStack().isErased(*this)); // this also check the allowed bounds of index_
      return *this;
    }
    StackIteratorInterface operator++(int) {
      StackIteratorInterface tmp(*this);
      do {
        ++index_;
      } while (
          getStack().isErased(*this)); // this also check the allowed bounds of index_
      return tmp;
    }
    StackIteratorInterface operator+(int delta) const {
      return StackIteratorInterface(*data_, index_ + delta);
    }
    bool operator==(StackIteratorInterface const& rhs) const {
      return index_ == rhs.index_;
    }
    bool operator!=(StackIteratorInterface const& rhs) const {
      return index_ != rhs.index_;
    }
    bool operator==(const ConstStackIteratorInterface<
                    TStackData, TParticleInterface, MSecondaryProducer,
                    stack_type>& rhs) const; // implemented below
    bool operator!=(const ConstStackIteratorInterface<
                    TStackData, TParticleInterface, MSecondaryProducer,
                    stack_type>& rhs) const; // implemented below

    /**
     * Convert iterator to value type, where value type is the user-provided particle
     * readout class.
     */
    particle_interface_type& operator*() {
      return static_cast<particle_interface_type&>(*this);
    }

    /**
     * Convert iterator to const value type, where value type is the user-provided
     * particle readout class.
     */
    particle_interface_type const& operator*() const {
      return static_cast<particle_interface_type const&>(*this);
    }
    //! @}

  protected:
    /**
     * @name Stack data access
     * @{
     */
    /// Get current particle index
    unsigned int getIndex() const { return index_; }
    /// Get current particle Stack object
    stack_type& getStack() { return *data_; }
    /// Get current particle const Stack object
    stack_type const& getStack() const { return *data_; }
    /// Get current user particle TStackData object
    TStackData& getStackData() { return data_->getStackData(); }
    /// Get current const user particle TStackData object
    TStackData const& getStackData() const { return data_->getStackData(); }
    /// Get data index as mapped in Stack class
    unsigned int getIndexFromIterator() const {
      return data_->getIndexFromIterator(index_);
    }
    //! @}

    // friends are needed for access to protected methods
    friend class Stack<TStackData, TParticleInterface,
                       MSecondaryProducer>; // for access to getIndex for Stack
    friend class Stack<TStackData&, TParticleInterface,
                       MSecondaryProducer>;            // for access to getIndex
                                                       // SecondaryView : public Stack
    friend class ParticleBase<StackIteratorInterface>; // for access to getStackDataType

    template <typename T1,                     // best fix this to: TStackData,
              template <typename> typename M1, // best fix this to: TParticleInterface,
              template <typename T, template <typename> typename T3> typename M2>
    friend class SecondaryView; // access grant for SecondaryView

    template <typename T, template <typename> typename TParticleInterface_>
    friend class corsika::history::HistorySecondaryProducer;

    friend class ConstStackIteratorInterface<TStackData, TParticleInterface,
                                             MSecondaryProducer, stack_type>;

  protected:
    unsigned int index_ = 0;

  private:
    stack_type* data_ = 0; // info: Particles and StackIterators become invalid when
                           // parent Stack is copied or deleted!

  }; // end class StackIterator

  /**
   * This is the iterator class for const-access to stack data.
   *
   * The const counterpart of StackIteratorInterface, which is used
   * for read-only iterator access on particle stack:
   *
   * @code
   * for (auto const& p : theStack) { E += p.getEnergy(); }
   * @endcode
   *
   * See documentation of StackIteratorInterface for more details:
   * \sa StackIteratorInterface.
   */

  template <typename TStackData, template <typename> typename TParticleInterface,
            template <typename T1, template <class> class T2> class MSecondaryProducer,
            typename TStackType =
                Stack<TStackData, TParticleInterface, MSecondaryProducer>>
  class ConstStackIteratorInterface
      : public TParticleInterface<ConstStackIteratorInterface<
            TStackData, TParticleInterface, MSecondaryProducer, TStackType>> {

  public:
    typedef TParticleInterface<ConstStackIteratorInterface<
        TStackData, TParticleInterface, MSecondaryProducer, TStackType>>
        particle_interface_type;

    typedef TStackType stack_type;
    typedef TStackData stack_data_type;

    // we don't want to allow dangling iterators to exist
    ConstStackIteratorInterface() = delete; //! \todo check rule of five

  public:
    ConstStackIteratorInterface(ConstStackIteratorInterface&& rhs)
        : index_(std::move(rhs.index_))
        , data_(std::move(rhs.data_)) {}

    ConstStackIteratorInterface(stack_type const& data, unsigned int const index)
        : index_(index)
        , data_(&data) {}

    bool isErased() const { return getStack().isErased(*this); }

    /**
     * @name Iterator interface
     * @{
     */
    ConstStackIteratorInterface& operator++() {
      do {
        ++index_;
      } while (
          getStack().isErased(*this)); // this also check the allowed bounds of index_
      return *this;
    }
    ConstStackIteratorInterface operator++(int) {
      ConstStackIteratorInterface tmp(*this);
      do {
        ++index_;
      } while (
          getStack().isErased(*this)); // this also check the allowed bounds of index_
      return tmp;
    }
    ConstStackIteratorInterface operator+(int const delta) const {
      return ConstStackIteratorInterface(*data_, index_ + delta);
    }
    bool operator==(ConstStackIteratorInterface const& rhs) const {
      return index_ == rhs.index_;
    }
    bool operator!=(ConstStackIteratorInterface const& rhs) const {
      return index_ != rhs.index_;
    }
    bool operator==(
        StackIteratorInterface<stack_data_type, TParticleInterface, MSecondaryProducer,
                               stack_type> const& rhs) const {
      return index_ == rhs.index_;
    }
    bool operator!=(
        StackIteratorInterface<stack_data_type, TParticleInterface, MSecondaryProducer,
                               stack_type> const& rhs) const {
      return index_ != rhs.index_;
    }

    particle_interface_type const& operator*() const {
      return static_cast<particle_interface_type const&>(*this);
    }
    //! @}

  protected:
    /**
     * @name Stack data access
     * Only the const versions for read-only access
     * @{
     */
    unsigned int getIndex() const { return index_; }
    stack_type const& getStack() const { return *data_; }
    stack_data_type const& getStackData() const { return data_->getStackData(); }
    /// Get data index as mapped in Stack class
    unsigned int getIndexFromIterator() const {
      return data_->getIndexFromIterator(index_);
    }
    //! @}

    // friends are needed for access to protected methods
    friend class Stack<stack_data_type, TParticleInterface,
                       MSecondaryProducer>; // for access to GetIndex for Stack
    friend class Stack<stack_data_type&, TParticleInterface,
                       MSecondaryProducer>; // for access to GetIndex

    friend class ParticleBase<ConstStackIteratorInterface>; // for access to GetStackData

    template <typename T1,                     // best fix to: stack_data_type,
              template <typename> typename M1, // best fix to: TParticleInterface,
              template <class T2, template <class> class T3> class M2>
    friend class SecondaryView; // access for SecondaryView

    friend class StackIteratorInterface<stack_data_type, TParticleInterface,
                                        MSecondaryProducer, stack_type>;

    template <typename T, template <typename> typename TParticleInterface_>
    friend class corsika::history::HistorySecondaryProducer;

  protected:
    unsigned int index_ = 0;

  private:
    stack_type const* data_ = 0; // info: Particles and StackIterators become invalid when
                                 // parent Stack is copied or deleted!

  }; // end class ConstStackIterator

} // namespace corsika
