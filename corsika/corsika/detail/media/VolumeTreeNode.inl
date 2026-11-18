/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/geometry/IVolume.hpp>
#include <corsika/media/interfaces/IMediumModel.hpp>

namespace corsika {
  namespace media {

    template <typename IModelProperties>
    inline bool VolumeTreeNode<IModelProperties>::contains(Point const& p) const {
      return geoVolume_->contains(p);
    }

    template <typename IModelProperties>
    inline VolumeTreeNode<IModelProperties> const*
    VolumeTreeNode<IModelProperties>::excludes(Point const& p) const {
      auto exclContainsIter =
          std::find_if(excludedNodes_.cbegin(), excludedNodes_.cend(),
                       [&](auto const& s) { return bool(s->contains(p)); });

      return exclContainsIter != excludedNodes_.cend() ? *exclContainsIter : nullptr;
    }

    /** returns a pointer to the sub-VolumeTreeNode which is "responsible" for the given
     * \class Point \p p, or nullptr iff \p p is not contained in this volume.
     */
    template <typename IModelProperties>
    inline VolumeTreeNode<IModelProperties> const*
    VolumeTreeNode<IModelProperties>::getContainingNode(Point const& p) const {
      if (!contains(p)) { return nullptr; }

      if (auto const childContainsIter =
              std::find_if(childNodes_.cbegin(), childNodes_.cend(),
                           [&](auto const& s) { return bool(s->contains(p)); });
          childContainsIter == childNodes_.cend()) // not contained in any of the children
      {
        if (auto const exclContainsIter = excludes(p)) // contained in any excluded nodes
        {
          return exclContainsIter->getContainingNode(p);
        } else {
          return this;
        }
      } else {
        return (*childContainsIter)->getContainingNode(p);
      }
    }

    template <typename IModelProperties>
    inline VolumeTreeNode<IModelProperties>*
    VolumeTreeNode<IModelProperties>::getContainingNode(Point const& p) {
      // see Scott Meyers, Effective C++ 3rd ed., Item 3
      return const_cast<VTN_type*>(std::as_const(*this).getContainingNode(p));
    }

    template <typename IModelProperties>
    inline void VolumeTreeNode<IModelProperties>::addChildToContainingNode(
        Point const& p, VTNUPtr pChild) {
      VolumeTreeNode<IModelProperties>* node = getContainingNode(p);
      if (!node) {
        CORSIKA_LOG_ERROR("Adding child at {} failed!. No containing node", p);
        throw std::runtime_error(
            "Failed adding child node. No parent at chosen location");
      }
      node->addChild(std::move(pChild));
    }

    template <typename IModelProperties>
    template <typename TCallable, bool preorder>
    inline void VolumeTreeNode<IModelProperties>::walk(TCallable func) const {
      if constexpr (preorder) { func(*this); }

      std::for_each(childNodes_.begin(), childNodes_.end(),
                    [&](auto const& v) { v->walk(func); });

      if constexpr (!preorder) { func(*this); };
    }

    template <typename IModelProperties>
    inline void VolumeTreeNode<IModelProperties>::addChild(
        typename VolumeTreeNode<IModelProperties>::VTNUPtr pChild) {
      pChild->parentNode_ = this;
      childNodes_.push_back(std::move(pChild));
      // It is a bad idea to return an iterator to the inserted element
      // because it might get invalidated when the vector needs to grow
      // later and the caller won't notice.
    }

    template <typename IModelProperties>
    inline void VolumeTreeNode<IModelProperties>::excludeOverlapWith(
        typename VolumeTreeNode<IModelProperties>::VTNUPtr const& pNode) {
      excludedNodes_.push_back(pNode.get());
    }

  } // namespace media
} // namespace corsika
