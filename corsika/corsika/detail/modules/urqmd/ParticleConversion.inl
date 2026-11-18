/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/modules/urqmd/ParticleConversion.hpp>

#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>

#include <urqmd.hpp>

namespace corsika::urqmd {

  inline bool canInteract(Code const vCode) {
    // According to the manual, UrQMD can use all mesons, baryons and nucleons
    // which are modeled also as input particles. I think it is safer to accept
    // only the usual long-lived species as input.
    // TODO: Charmed mesons should be added to the list, too

    static std::array constexpr validProjectileCodes{
        Code::Proton,  Code::AntiProton, Code::Neutron, Code::AntiNeutron, Code::PiPlus,
        Code::PiMinus, Code::KPlus,      Code::KMinus,  Code::K0Short,     Code::K0Long};

    return std::find(std::cbegin(validProjectileCodes), std::cend(validProjectileCodes),
                     vCode) != std::cend(validProjectileCodes);
  }

  inline std::pair<int, int> convertToUrQMD(Code const code) {
    static const std::map<int, std::pair<int, int>> mapPDGToUrQMD{
        // data mostly from github.com/afedynitch/ParticleDataTool
        {22, {100, 0}},      // photon
        {111, {101, 0}},     // pi0
        {211, {101, 2}},     // pi+
        {-211, {101, -2}},   // pi-
        {321, {106, 1}},     // K+
        {-321, {-106, -1}},  // K-
        {311, {106, -1}},    // K0
        {-311, {-106, 1}},   // K0bar
        {2212, {1, 1}},      // p
        {2112, {1, -1}},     // n
        {-2212, {-1, -1}},   // pbar
        {-2112, {-1, 1}},    // nbar
        {221, {102, 0}},     // eta
        {213, {104, 2}},     // rho+
        {-213, {104, -2}},   // rho-
        {113, {104, 0}},     // rho0
        {323, {108, 2}},     // K*+
        {-323, {108, -2}},   // K*-
        {313, {108, 0}},     // K*0
        {-313, {-108, 0}},   // K*0-bar
        {223, {103, 0}},     // omega
        {333, {109, 0}},     // phi
        {3222, {40, 2}},     // Sigma+
        {3212, {40, 0}},     // Sigma0
        {3112, {40, -2}},    // Sigma-
        {3322, {49, 0}},     // Xi0
        {3312, {49, -1}},    // Xi-
        {3122, {27, 0}},     // Lambda0
        {2224, {17, 4}},     // Delta++
        {2214, {17, 2}},     // Delta+
        {2114, {17, 0}},     // Delta0
        {1114, {17, -2}},    // Delta-
        {3224, {41, 2}},     // Sigma*+
        {3214, {41, 0}},     // Sigma*0
        {3114, {41, -2}},    // Sigma*-
        {3324, {50, 0}},     // Xi*0
        {3314, {50, -1}},    // Xi*-
        {3334, {55, 0}},     // Omega-
        {411, {133, 2}},     // D+
        {-411, {133, -2}},   // D-
        {421, {133, 0}},     // D0
        {-421, {-133, 0}},   // D0-bar
        {441, {107, 0}},     // etaC
        {431, {138, 1}},     // Ds+
        {-431, {138, -1}},   // Ds-
        {433, {139, 1}},     // Ds*+
        {-433, {139, -1}},   // Ds*-
        {413, {134, 1}},     // D*+
        {-413, {134, -1}},   // D*-
        {10421, {134, 0}},   // D*0
        {-10421, {-134, 0}}, // D*0-bar
        {443, {135, 0}},     // jpsi
    };

    return mapPDGToUrQMD.at(static_cast<int>(get_PDG(code)));
  }

  inline Code convertFromUrQMD(int vItyp, int vIso3) {
    int const pdgInt =
        ::urqmd::pdgid_(vItyp, vIso3); // use the conversion function provided by UrQMD
    if (pdgInt == 0) {                 // ::urqmd::pdgid_ returns 0 on error
      throw std::runtime_error("UrQMD pdgid() returned 0");
    }
    auto const pdg = static_cast<PDGCode>(pdgInt);
    return convert_from_PDG(pdg);
  }

} // namespace corsika::urqmd
