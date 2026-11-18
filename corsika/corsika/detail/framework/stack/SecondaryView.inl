/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/stack/Stack.hpp>
#include <corsika/framework/core/Logging.hpp>

#include <stdexcept>
#include <vector>

namespace corsika {
  template <typename TStackDataType, template <typename> typename TParticleInterface,
            template <typename T1, template <class> class T2> class MSecondaryProducer>
  template <typename... Args>
  inline typename SecondaryView<TStackDataType, TParticleInterface,
                                MSecondaryProducer>::stack_view_iterator
  SecondaryView<TStackDataType, TParticleInterface, MSecondaryProducer>::addSecondary(
      const Args... v) {
    CORSIKA_LOG_TRACE("SecondaryView::addSecondary(Args&&)");
    stack_view_iterator proj = getProjectile(); // make this const
    return addSecondary(proj, v...);
  }

  template <typename TStackDataType, template <typename> typename TParticleInterface,
            template <typename T1, template <class> class T2> class MSecondaryProducer>
  inline typename SecondaryView<TStackDataType, TParticleInterface,
                                MSecondaryProducer>::stack_view_iterator
  SecondaryView<TStackDataType, TParticleInterface, MSecondaryProducer>::begin() {
    unsigned int i = 0;
    for (; i < getSize(); ++i) {

      if (!isErased(i)) break;
    }
    return stack_view_iterator(*this, i + 1);
  }

  template <typename TStackDataType, template <typename> typename TParticleInterface,
            template <typename T1, template <class> class T2> class MSecondaryProducer>
  inline typename SecondaryView<TStackDataType, TParticleInterface,
                                MSecondaryProducer>::const_stack_view_iterator
  SecondaryView<TStackDataType, TParticleInterface, MSecondaryProducer>::begin() const {
    unsigned int i = 0;
    for (; i < getSize(); ++i) {
      if (!isErased(i)) break;
    }

    return const_stack_view_iterator(*this, i + 1);
  }

  template <typename TStackDataType, template <typename> typename TParticleInterface,
            template <typename T1, template <class> class T2> class MSecondaryProducer>
  inline typename SecondaryView<TStackDataType, TParticleInterface,
                                MSecondaryProducer>::const_stack_view_iterator
  SecondaryView<TStackDataType, TParticleInterface, MSecondaryProducer>::cbegin() const {
    unsigned int i = 0;
    for (; i < getSize(); ++i) {
      if (!isErased(i)) break;
    }

    return const_stack_view_iterator(*this, i + 1);
  }

  template <typename TStackDataType, template <typename> typename TParticleInterface,
            template <typename T1, template <class> class T2> class MSecondaryProducer>
  inline typename SecondaryView<TStackDataType, TParticleInterface,
                                MSecondaryProducer>::stack_view_iterator
  SecondaryView<TStackDataType, TParticleInterface, MSecondaryProducer>::last() {
    unsigned int i = 0;
    for (; i < getSize(); ++i) {
      if (!isErased(getSize() - 1 - i)) break;
    }
    return stack_view_iterator(*this, getSize() - 1 - i + 1);
  }

  template <typename TStackDataType, template <typename> typename TParticleInterface,
            template <typename T1, template <class> class T2> class MSecondaryProducer>
  inline typename SecondaryView<TStackDataType, TParticleInterface,
                                MSecondaryProducer>::const_stack_view_iterator
  SecondaryView<TStackDataType, TParticleInterface, MSecondaryProducer>::last() const {
    unsigned int i = 0;
    for (; i < getSize(); ++i) {
      if (!isErased(getSize() - 1 - i)) break;
    }
    return stack_view_iterator(*this, getSize() - 1 - i + 1);
  }

  template <typename TStackDataType, template <typename> typename TParticleInterface,
            template <typename T1, template <class> class T2> class MSecondaryProducer>
  inline typename SecondaryView<TStackDataType, TParticleInterface,
                                MSecondaryProducer>::const_stack_view_iterator
  SecondaryView<TStackDataType, TParticleInterface, MSecondaryProducer>::clast() const {
    unsigned int i = 0;
    for (; i < getSize(); ++i) {

      if (!isErased(getSize() - 1 - i)) break;
    }
    return const_stack_view_iterator(*this, getSize() - 1 - i + 1);
  }

  template <typename TStackDataType, template <typename> typename TParticleInterface,
            template <typename T1, template <class> class T2> class MSecondaryProducer>
  inline void SecondaryView<TStackDataType, TParticleInterface, MSecondaryProducer>::swap(
      stack_view_iterator a, stack_view_iterator b) {

    CORSIKA_LOG_TRACE("View::swap");
    inner_stack_.swap(getIndexFromIterator(a.getIndex()),
                      getIndexFromIterator(b.getIndex()));
  }

  template <typename TStackDataType, template <typename> typename TParticleInterface,
            template <typename T1, template <class> class T2> class MSecondaryProducer>
  inline void SecondaryView<TStackDataType, TParticleInterface, MSecondaryProducer>::copy(
      stack_view_iterator a, stack_view_iterator b) {

    CORSIKA_LOG_TRACE("View::copy");
    inner_stack_.copy(getIndexFromIterator(a.getIndex()),
                      getIndexFromIterator(b.getIndex()));
  }

  template <typename TStackDataType, template <typename> typename TParticleInterface,
            template <typename T1, template <class> class T2> class MSecondaryProducer>
  inline void SecondaryView<TStackDataType, TParticleInterface, MSecondaryProducer>::copy(
      const_stack_view_iterator a, stack_view_iterator b) {

    CORSIKA_LOG_TRACE("View::copy");
    inner_stack_.copy(getIndexFromIterator(a.getIndex()),
                      getIndexFromIterator(b.getIndex()));
  }

  template <typename TStackDataType, template <typename> typename TParticleInterface,
            template <typename T1, template <class> class T2> class MSecondaryProducer>
  inline void SecondaryView<TStackDataType, TParticleInterface,
                            MSecondaryProducer>::erase(stack_view_iterator p) {

    CORSIKA_LOG_TRACE("SecondaryView::Delete");
    if (isEmpty()) { /*error*/
      throw std::runtime_error("Stack, cannot delete entry since size is zero");
    }
    if (isErased(p.getIndex() - 1)) { /*error*/
      throw std::runtime_error("Stack, cannot delete entry since already deleted");
    }
    inner_stack_.erase(getIndexFromIterator(p.getIndex()));
    inner_stack_reference_type::nDeleted_++; // also count in SecondaryView
  }

  template <typename TStackDataType, template <typename> typename TParticleInterface,
            template <typename T1, template <class> class T2> class MSecondaryProducer>
  inline bool SecondaryView<TStackDataType, TParticleInterface,
                            MSecondaryProducer>::purgeLastIfDeleted() {
    CORSIKA_LOG_TRACE("SecondaryView::purgeLastIfDeleted");
    if (!isErased(getSize() - 1))
      return false; // the last particle is not marked for deletion. Do nothing.
    inner_stack_.purge(getIndexFromIterator(getSize()));
    inner_stack_reference_type::nDeleted_--;
    indices_.pop_back();
    return true;
  }

  template <typename TStackDataType, template <typename> typename TParticleInterface,
            template <typename T1, template <class> class T2> class MSecondaryProducer>
  inline void
  SecondaryView<TStackDataType, TParticleInterface, MSecondaryProducer>::purge() {
    unsigned int iStack = 0;
    unsigned int size = getSize();
    while (iStack < size) {

      if (isErased(iStack)) {
        inner_stack_.purge(iStack);
        indices_.erase(indices_.begin() + iStack);
      }
      size = getSize();
      iStack++;
    }
    inner_stack_reference_type::nDeleted_ = 0;
  }

  template <typename TStackDataType, template <typename> typename TParticleInterface,
            template <typename T1, template <class> class T2> class MSecondaryProducer>
  inline std::string SecondaryView<TStackDataType, TParticleInterface,
                                   MSecondaryProducer>::asString() const {
    std::string str(fmt::format("size {}\n", getSize()));
    // we make our own begin/end since we want ALL entries
    std::string new_line = "     ";
    for (unsigned int iPart = 0; iPart != getSize(); ++iPart) {
      const_stack_view_iterator itPart(*this, iPart);
      str += fmt::format(
          "{}{}{}", new_line, itPart.asString(),
          (inner_stack_.deleted_[getIndexFromIterator(itPart.getIndex())] ? " [deleted]"
                                                                          : ""));
      new_line = "\n     ";
    }
    return str;
  }

  template <typename TStackDataType, template <typename> typename TParticleInterface,
            template <typename T1, template <class> class T2> class MSecondaryProducer>
  template <typename... Args>
  inline typename SecondaryView<TStackDataType, TParticleInterface,
                                MSecondaryProducer>::stack_view_iterator
  SecondaryView<TStackDataType, TParticleInterface, MSecondaryProducer>::addSecondary(
      stack_view_iterator& proj, const Args... v) {
    CORSIKA_LOG_TRACE("SecondaryView::addSecondary(stack_view_iterator&, Args&&)");
    // make space on stack
    inner_stack_reference_type::getStackData().incrementSize();
    inner_stack_.deleted_.push_back(false);
    // get current number of secondaries on stack
    const unsigned int idSec = getSize();
    // determine index on (inner) stack where new particle will be located
    const unsigned int index = inner_stack_reference_type::getStackData().getSize() - 1;
    indices_.push_back(index);
    // NOTE: "+1" is since "0" is special marker here for PROJECTILE, see
    // getIndexFromIterator
    auto sec = stack_view_iterator(*this, idSec + 1, proj, v...);
    MSecondaryProducer<TStackDataType, TParticleInterface>::new_secondary(sec);
    return sec;
  }

} // namespace corsika
