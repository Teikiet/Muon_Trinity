/*
 * (c) Copyright 2023 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <Tauola/TauolaParticle.h>
#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

namespace corsika::tauola {

  /**
   * @class TauolaInterfaceParticle
   *
   *
   * Since TAUOLA requires a concrete instance of the
   * TauolaParticle class in order to decay, we implement a
   * simple class that is used to perform the Tau decay and store
   * the tree of decay products.
   *
   * This is a relatively lightweight class and has methods
   * to easily convert to/from the CORSIKA8 particle representation.
   *
   */
  class TauolaInterfaceParticle final : public Tauolapp::TauolaParticle {

    int pdg_{0};     ///< The PDG ID of this particle.
    int status_{0};  ///< The status code of this particle.
    double mass_{0}; ///< The mass of this particle.

    // we store the four-vector as independent components in unit of GeV.
    double Px_{0.}; ///< The x-component of the four momentum.
    double Py_{0.}; ///< The y-component of the four momentum.
    double Pz_{0.}; ///< The z-component of the four momentum.
    double E_{0.};  ///< The energy of the four momentum.

    std::vector<TauolaParticle*> mothers_;   ///< The mothers of this particle.
    std::vector<TauolaParticle*> daughters_; ///< The daughters of this particle.
    std::vector<TauolaParticle*> created_;   /// The particles created by this particle.

  public:
    /**
     * Returns the x component of the momentum four vector.
     */
    double getPx() final override { return Px_; };

    /**
     * Returns the y component of the momentum four vector.
     */
    double getPy() final override { return Py_; };

    /**
     * Returns the z component of the momentum four vector.
     */
    double getPz() final override { return Pz_; };

    /**
     * Returns the energy component of the momentum four vector.
     */
    double getE() final override { return E_; };

    /**
     * Set the x component of the momentum four vector.
     */
    void setPx(double px) final override { Px_ = px; };

    /**
     * Set the y component of the momentum four vector.
     */
    void setPy(double py) final override { Py_ = py; };

    /**
     * Set the z component of the momentum four vector.
     */
    void setPz(double pz) final override { Pz_ = pz; };

    /**
     * Set the x component of the momentum four vector.
     */
    void setE(double e) final override { E_ = e; };

    /**
     * Set the mothers of this particle.
     *
     * @param mothers    A vector of TauolaInterfaceParticles to assign as mothers.
     */
    void setMothers(std::vector<TauolaParticle*> mothers) final override {

      // set our mother records to point the new list of mothers
      mothers_ = mothers;

      // add this particle as a daughter to all the mothers
      for (auto it = mothers.begin(); it < mothers.end(); ++it) {

        // get a TauolaInterfaceParticle reference
        auto particle{dynamic_cast<TauolaInterfaceParticle*>(*it)};

        // and add us to the list of daughters
        particle->addDaughter(this);
      }
    };

    /**
     * Set the daughters of this particle.
     *
     * @param daughters    A vector of TauolaInterfaceParticles to assign as daughters.
     *
     */
    void setDaughters(std::vector<TauolaParticle*> daughters) final override {

      // set our daughter list to the list of daughters
      daughters_ = daughters;
    };

    /**
     * Returns the mothers of this particle.
     *
     * @returns The mothers of this TauolaInterfaceParticle.
     *
     */
    std::vector<TauolaParticle*> getMothers() final override { return mothers_; };

    /**
     * Returns the daughters of this particle.
     *
     * @returns The daughters of this TauolaInterfaceParticle.
     *
     */
    std::vector<TauolaParticle*> getDaughters() final override { return daughters_; };

    /**
     * Returns the particles created by this particle.
     *
     * @returns The TauolaInterfaceParticle's created by this particle.
     *
     */
    std::vector<TauolaParticle*> getCreated() { return created_; };

    /**
     * Set the particles created by this particle.
     *
     * @param created    A vector of TauolaInterfaceParticles.
     *
     */
    void setCreated(std::vector<TauolaParticle*> created) { created_ = created; };

    /**
     * Set the PDG ID of this particle.
     *
     * @param pdg_id   The new particle's PDG ID.
     */
    void setPdgID(int pdg_id) final override { pdg_ = pdg_id; };

    /**
     * Set the mass of this particle.
     *
     * @param mass   The new particle mass (Gev/c^2).
     *
     */
    void setMass(double mass) final override { mass_ = mass; };

    /**
     * Set the status of this particle.
     *
     * @param status   The new particle status to set.
     *
     */
    void setStatus(int status) final override { status_ = status; };

    /**
     * Get the PDG ID of this particle.
     *
     * @returns The PDG ID of the current particle.
     *
     */
    int getPdgID() final override { return pdg_; };

    /**
     * Get the status of this particle.
     *
     * @returns The status of the current particle.
     *
     */
    int getStatus() final override { return status_; };

    /**
     * Get the barcode of this particle.
     *
     * This is unused in this implementation.
     *
     */
    int getBarcode() final override { return -1; };

    /**
     * Find the particle mass from Tauolapp::amxx
     */
    static HEPEnergyType findMassTauola(Code pid) {
      switch (pid) {
        case Code::TauPlus:
        case Code::TauMinus:
          return Tauolapp::parmas_.amtau * 1_GeV;

        case Code::RhoPlus:
        case Code::RhoMinus:
          return Tauolapp::parmas_.amro * 1_GeV;

        case Code::KStarPlus:
        case Code::KStarMinus:
          return Tauolapp::parmas_.amkst * 1_GeV;

        default:
          return get_mass(pid);
      }
    }

    /**
     * Add a daughter to this particle.
     *
     * @param particle    The daughter particle to store in our tree.
     */
    void addDaughter(TauolaParticle* particle) { daughters_.push_back(particle); };

    /**
     * Create a new particle with the given parameters.
     *
     * @param pdg_id   The new particle's PDG ID.
     * @param status   The new particle' status.
     * @param mass     The new particle's mass (GeV/c^2).
     * @param px       The new particle's x-momentum (GeV/c).
     * @param py       The new particle's y-momentum (GeV/c).
     * @param pz       The new particle's z-momentum (GeV/c).
     * @param e        The new particle's energy (GeV).
     *
     * @returns A new TauolaInterfaceParticle instance.
     *
     */
    TauolaParticle* createNewParticle(int pdg_id, int status, double mass, double px,
                                      double py, double pz, double e) final override {

      // create a new particle
      auto particle{new TauolaInterfaceParticle(pdg_id, status, mass, px, py, pz, e)};

      // save it onto our vector of created particles
      created_.push_back(particle);

      // and return it to the caller
      return particle;
    };

    /**
     * Construct a new TauolaInterfaceParticle instance.
     *
     * @param pdg_id   The particle's PDG ID.
     * @param status   The particle' status.
     * @param mass     The particle's mass (GeV/c^2).
     * @param px       The particle's x-momentum (GeV/c).
     * @param py       The particle's y-momentum (GeV/c).
     * @param pz       The particle's z-momentum (GeV/c).
     * @param e        The particle's energy (GeV).
     *
     */
    TauolaInterfaceParticle(const int pdg, const int status, const double mass,
                            const double px, const double py, const double pz,
                            const double e)
        : pdg_(pdg)
        , status_(status)
        , mass_(mass)
        , Px_(px)
        , Py_(py)
        , Pz_(pz)
        , E_(e){};

    /**
     * Copy an existing TauolaInterfaceParticle instance.
     *
     * @param other    An existing TauolaInterfaceParticle instance.
     *
     */
    TauolaInterfaceParticle(TauolaInterfaceParticle& other)
        : pdg_(other.getPdgID())
        , status_(other.getStatus())
        , mass_(other.getMass())
        , Px_(other.getPx())
        , Py_(other.getPy())
        , Pz_(other.getPz())
        , E_(other.getE()){};

    /**
     * Copy-assign from an existing TauolaInterfaceParticle instance.
     */
    TauolaInterfaceParticle& operator=(TauolaInterfaceParticle other) {

      // set the copyable members
      pdg_ = other.getPdgID();
      status_ = other.getStatus();
      mass_ = other.getMass();
      Px_ = other.getPx();
      Py_ = other.getPy();
      Pz_ = other.getPz();
      E_ = other.getE();

      // get references to the particle lists of other
      auto other_mothers{other.getMothers()};
      auto other_daughters{other.getDaughters()};
      auto other_created{other.getCreated()};

      // assign our values to other so they are properly deleted
      other.setMothers(this->mothers_);
      other.setDaughters(this->daughters_);
      other.setCreated(this->created_);

      // and assign other's values to us
      this->setMothers(other_mothers);
      this->setDaughters(other_daughters);
      this->setCreated(other_created);

      // and we have copy assigned.
      return *this;
    };

    /**
     * Construct a TauolaInterfaceParticle from a CORSIKA projectile.
     */
    template <typename TProjectile>
    static std::unique_ptr<TauolaInterfaceParticle> from_projectile(
        const TProjectile& projectile) {
      using namespace units::si;

      // get the particles PDG ID
      const auto pdg{static_cast<int>(get_PDG(projectile.getPID()))};

      // get the momentum in the lab frame
      auto const& labMomentum{projectile.getMomentum()};

      // get the typed mass
      const auto mass{findMassTauola(projectile.getPID())};

      // get the particle parameters in GeV (the units that we set for TAUOLA)
      const double px{labMomentum.getComponents()[0] / 1_GeV};
      const double py{labMomentum.getComponents()[1] / 1_GeV};
      const double pz{labMomentum.getComponents()[2] / 1_GeV};

      // compute the total energy of the particle
      const double E{sqrt(labMomentum.getSquaredNorm() + mass * mass) / 1_GeV};
      // const double E{projectile.getEnergy()/1_GeV};

      // we start each particle with a stable status
      const auto status{Tauolapp::TauolaParticle::STABLE};

      // construct the TaulaInterfaceParticle as a unique pointer
      return std::make_unique<TauolaInterfaceParticle>(pdg, status, mass / 1_GeV, px, py,
                                                       pz, E);
    }

    /*
     * Free all the memory of mothers, daughters, and created particles.
     */
    ~TauolaInterfaceParticle() {

      // delete all the particles that we created
      while (created_.size() != 0) {
        auto temp{created_.back()};
        created_.pop_back();
        delete temp;
      }
    };

    /**
     * Print this particle to std::out.
     */
    void print() final override {
      char buf[256];
      sprintf(buf, "P: %6i %2i | %11.4e %11.4e %11.4e %11.4e | %11.4e |\n", pdg_, status_,
              Px_, Py_, Pz_, E_, mass_);

      cout << buf;
    };
  };

} // namespace corsika::tauola
