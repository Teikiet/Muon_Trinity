/*
 * (c) Copyright 2025 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/Vector.hpp>
#include <corsika/framework/stack/Stack.hpp>
#include <corsika/modules/epos-lhcr/ParticleConversion.hpp>

#include <epos-lhcr-public.hpp>

namespace corsika::EPOS_LHCR {

  typedef corsika::Vector<hepmomentum_d> MomentumVector;

  class EposStackData {

  public:
    void dump() const {}

    void clear() { ::EPOS_LHCR::cptl_->nptl = 0; }
    unsigned int getSize() const { return ::EPOS_LHCR::cptl_->nptl; }
    unsigned int getCapacity() const { return ::EPOS_LHCR::mxptl; }

    void setId(const unsigned int i, const int v) { ::EPOS_LHCR::cptl_->idptl[i] = v; }
    void setEnergy(const unsigned int i, const HEPEnergyType v) {
      ::EPOS_LHCR::cptl_->pptl[i][3] = v / 1_GeV;
    }
    void setMass(const unsigned int i, const HEPMassType v) {
      ::EPOS_LHCR::cptl_->pptl[i][4] = v / 1_GeV;
    }
    void setMomentum(const unsigned int i, const MomentumVector& v) {
      auto tmp = v.getComponents();
      for (int idx = 0; idx < 3; ++idx)
        ::EPOS_LHCR::cptl_->pptl[i][idx] = tmp[idx] / 1_GeV;
    }
    void setState(const unsigned int i, const int v) {
      ::EPOS_LHCR::cptl_->istptl[i] = v;
    }

    int getId(const unsigned int i) const { return ::EPOS_LHCR::cptl_->idptl[i]; }
    int getState(const unsigned int i) const { return ::EPOS_LHCR::cptl_->istptl[i]; }
    HEPEnergyType getEnergy(const int i) const {
      return ::EPOS_LHCR::cptl_->pptl[i][3] * 1_GeV;
    }
    HEPEnergyType getMass(const unsigned int i) const {
      return ::EPOS_LHCR::cptl_->pptl[i][4] * 1_GeV;
    }
    MomentumVector getMomentum(const unsigned int i) const {
      QuantityVector<hepmomentum_d> components = {::EPOS_LHCR::cptl_->pptl[i][0] * 1_GeV,
                                                  ::EPOS_LHCR::cptl_->pptl[i][1] * 1_GeV,
                                                  ::EPOS_LHCR::cptl_->pptl[i][2] * 1_GeV};
      return MomentumVector(components);
    }

    MomentumVector getMomentum(const unsigned int i, const CoordinateSystem& CS) const {
      QuantityVector<hepmomentum_d> components = {::EPOS_LHCR::cptl_->pptl[i][0] * 1_GeV,
                                                  ::EPOS_LHCR::cptl_->pptl[i][1] * 1_GeV,
                                                  ::EPOS_LHCR::cptl_->pptl[i][2] * 1_GeV};
      return MomentumVector(CS, components);
    }

    QuantityVector<hepmomentum_d> getRawMomentum(const unsigned int i) const {
      return QuantityVector<hepmomentum_d>{::EPOS_LHCR::cptl_->pptl[i][0] * 1_GeV,
                                           ::EPOS_LHCR::cptl_->pptl[i][1] * 1_GeV,
                                           ::EPOS_LHCR::cptl_->pptl[i][2] * 1_GeV};
    }

    void copy(const unsigned int i1, const unsigned int i2) {
      ::EPOS_LHCR::cptl_->idptl[i2] = ::EPOS_LHCR::cptl_->idptl[i1];
      ::EPOS_LHCR::cptl_->iorptl[i2] = ::EPOS_LHCR::cptl_->iorptl[i1];
      ::EPOS_LHCR::cptl_->jorptl[i2] = ::EPOS_LHCR::cptl_->jorptl[i1];
      ::EPOS_LHCR::cptl_->istptl[i2] = ::EPOS_LHCR::cptl_->istptl[i1];
      ::EPOS_LHCR::cptl_->ityptl[i2] = ::EPOS_LHCR::cptl_->ityptl[i1];
      for (unsigned int i = 0; i < 5; ++i)
        ::EPOS_LHCR::cptl_->pptl[i2][i] = ::EPOS_LHCR::cptl_->pptl[i1][i];
      for (unsigned int i = 0; i < 2; ++i) {
        ::EPOS_LHCR::cptl_->tivptl[i2][i] = ::EPOS_LHCR::cptl_->tivptl[i1][i];
        ::EPOS_LHCR::cptl_->ifrptl[i2][i] = ::EPOS_LHCR::cptl_->ifrptl[i1][i];
      }
      for (unsigned int i = 0; i < 4; ++i) {
        ::EPOS_LHCR::cptl_->xorptl[i2][i] = ::EPOS_LHCR::cptl_->xorptl[i1][i];
        ::EPOS_LHCR::cptl_->ibptl[i2][i] = ::EPOS_LHCR::cptl_->ibptl[i1][i];
      }
    }

    void swap(const unsigned int i1, const unsigned int i2) {
      std::swap(::EPOS_LHCR::cptl_->idptl[i2], ::EPOS_LHCR::cptl_->idptl[i1]);
      std::swap(::EPOS_LHCR::cptl_->iorptl[i2], ::EPOS_LHCR::cptl_->iorptl[i1]);
      std::swap(::EPOS_LHCR::cptl_->jorptl[i2], ::EPOS_LHCR::cptl_->jorptl[i1]);
      std::swap(::EPOS_LHCR::cptl_->istptl[i2], ::EPOS_LHCR::cptl_->istptl[i1]);
      std::swap(::EPOS_LHCR::cptl_->ityptl[i2], ::EPOS_LHCR::cptl_->ityptl[i1]);
      for (unsigned int i = 0; i < 5; ++i)
        std::swap(::EPOS_LHCR::cptl_->pptl[i2][i], ::EPOS_LHCR::cptl_->pptl[i1][i]);
      for (unsigned int i = 0; i < 2; ++i) {
        std::swap(::EPOS_LHCR::cptl_->tivptl[i2][i], ::EPOS_LHCR::cptl_->tivptl[i1][i]);
        std::swap(::EPOS_LHCR::cptl_->ifrptl[i2][i], ::EPOS_LHCR::cptl_->ifrptl[i1][i]);
      }
      for (unsigned int i = 0; i < 4; ++i) {
        std::swap(::EPOS_LHCR::cptl_->xorptl[i2][i], ::EPOS_LHCR::cptl_->xorptl[i1][i]);
        std::swap(::EPOS_LHCR::cptl_->ibptl[i2][i], ::EPOS_LHCR::cptl_->ibptl[i1][i]);
      }
    }

    void incrementSize() { ::EPOS_LHCR::cptl_->nptl++; }
    void decrementSize() {
      if (::EPOS_LHCR::cptl_->nptl > 0) { ::EPOS_LHCR::cptl_->nptl--; }
    }
  };

  template <typename StackIteratorInterface>
  class ParticleInterface : public corsika::ParticleBase<StackIteratorInterface> {

    using corsika::ParticleBase<StackIteratorInterface>::getStackData;
    using corsika::ParticleBase<StackIteratorInterface>::getIndex;

  public:
    void setParticleData(const int vID, // corsika::EPOS_LHCR::EposCode vID,
                         const HEPEnergyType vE, const MomentumVector& vP,
                         const HEPMassType vM) {
      setPID(vID);
      setEnergy(vE);
      setMomentum(vP);
      setMass(vM);
      setState(0);
    }

    void setParticleData(ParticleInterface<StackIteratorInterface>& /*parent*/,
                         const int vID, //  corsika::EPOS_LHCR::EposCode vID,
                         const HEPEnergyType vE, const MomentumVector& vP,
                         const HEPMassType vM) {
      setPID(vID);
      setEnergy(vE);
      setMomentum(vP);
      setMass(vM);
      setState(0);
    }

    void setEnergy(const HEPEnergyType v) { getStackData().setEnergy(getIndex(), v); }

    HEPEnergyType getEnergy() const { return getStackData().getEnergy(getIndex()); }

    bool isFinal() const { return getStackData().getState(getIndex()) == 0; }

    void setMass(const HEPMassType v) { getStackData().setMass(getIndex(), v); }

    HEPEnergyType getMass() const { return getStackData().getMass(getIndex()); }

    void setPID(const int v) { getStackData().setId(getIndex(), v); }

    corsika::EPOS_LHCR::EposCode getPID() const {
      return static_cast<corsika::EPOS_LHCR::EposCode>(getStackData().getId(getIndex()));
    }

    void setState(const int v) { getStackData().setState(getIndex(), v); }

    corsika::EPOS_LHCR::EposCode getState() const {
      return static_cast<corsika::EPOS_LHCR::EposCode>(
          getStackData().getState(getIndex()));
    }

    MomentumVector getMomentum() const { return getStackData().getMomentum(getIndex()); }

    MomentumVector getMomentum(const CoordinateSystem& CS) const {
      return getStackData().getMomentum(getIndex(), CS);
    }

    QuantityVector<hepmomentum_d> getRawMomentum() const {
      return getStackData().getRawMomentum(getIndex());
    }

    void setMomentum(const MomentumVector& v) {
      getStackData().setMomentum(getIndex(), v);
    }
  };

  typedef corsika::Stack<EposStackData, ParticleInterface> EposStack;

} // end namespace corsika::EPOS_LHCR
