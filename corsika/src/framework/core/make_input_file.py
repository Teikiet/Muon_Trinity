#!/usr/bin/env python3

from collections import OrderedDict
import particle
import hepunits
from particle import Particle
import numpy as np
import os
import pickle
import re

ABS_PATH_HERE = os.path.abspath(os.path.dirname(__file__))

# These particles will be renamed to their more common ones
# i.e. gamma --> photon
# The comments give what `particle` would assign the name as
special_names = {
    0: "Unknown",
    11: "Electron",
    -11: "Positron",
    22: "Photon",
    113: "Rho0", # Rho_770_0
    213: "RhoPlus", # Rho_770_Plus
    -213: "RhoMinus", # Rho_770_Minus
    223: "Omega", # Omega_782
    130: "K0Long",
    310: "K0Short",
    313: "KStar0", # KStar_892_0
    -313: "KStar0Bar", # KStar_892_0Bar
    323: "KStarPlus", # KStar_892_Plus
    -323: "KStarMinus", # KStar_892_Minus
    331: "EtaPrime", # EtaPrime_958
    333: "Phi", # Phi_1020
    413: "DStarPlus", # DStar_2010_Plus
    -413: "DStarMinus", # DStar_2010_Minus
    423: "DStar0", # DStar_2007_0
    -423: "DStar0Bar", # DStar_2007_0Bar
    441: "EtaC", # EtaC1s
    443: "Jpsi", # Jpsi1s
    2112: "Neutron",
    -2112: "AntiNeutron",
    2212: "Proton",
    -2212: "AntiProton",
    2114: "Delta0", # Delta_1232_0
    -2114: "Delta0Bar", # Delta_1232_0Bar
}

# For making the conversion between pdg and the interal C8 enums,
# values > this number will be done using a map while values less
# than this will use an array. This number should encapsulate
# the pdg ids of common air shower particles
max_conversion_array_id = 500


def c_identifier_camel(name, remove_numbers=True):
    # Converts from snake case to camel case

    temp = str(name)
    if "(" in temp and ")" in temp:
        temp = temp.split()

    temp = temp.replace("st_", "_star_")
    temp = temp.replace("*", "_star_")
    temp = temp.replace("p_", "_prime_")
    temp = temp.replace("_pp", "_plus_plus")
    temp = temp.replace("_mm", "_minus_minus")
    parts = temp.split("_")
    for i in range(len(parts)):
        if parts[i].isdigit() and int(parts[i]) > 9:
            if i == len(parts) - 1:
                parts[i] = "_" + parts[i]
            else:
                parts[i] = f"_{parts[i]}_"
        parts[i] = parts[i][0].upper() + parts[i][1:].lower()
    temp = "".join(parts)

    return temp


##########################
### Read in the data files
##########################

# file containing a list of pdg ids
desired_ids = np.loadtxt(ABS_PATH_HERE + "/PDG_Ids.txt")
# file of known nuclei containing columns pdg, name, A, Z
nuclear_data = np.loadtxt(ABS_PATH_HERE + "/NuclearData.txt", dtype=str)


######################################################
### Create the particle information from supplied pdgs
######################################################

particle_db = {}

particle_db[special_names[0]] = {
    "name": "void",
    "antiName": special_names[0],
    "pdg": 0,
    "mass": 0,
    "charge": 0,
    "lifetime": np.inf,
    "ngc_code": 0,
    "isHadron": False,
    "isNucleus": False,
}
counter = len(particle_db)

for pdgid in desired_ids:
    part = Particle.from_pdgid(pdgid)
    name = part.programmatic_name
    name = name.replace("_plus", "+").replace("_minus", "-")

    if pdgid in special_names:
        c_id = special_names[pdgid]
        # print("Replacing", c_identifier_camel(part.programmatic_name), "with", c_id)
    else:
        c_id = c_identifier_camel(part.programmatic_name, abs(int(pdgid)) < 6000)

    if c_id in particle_db.keys():
        print("CID", c_id, "already in")
        print("Current ", pdgid, part.programmatic_name, c_id)
        print("Existing", particle_db[c_id])
        exit()

    if pdgid in[10333]:
        print(pdgid, part.name, part.programmatic_name, c_id)

    # Get the name of the antiparticle for this pdg
    anti_id = int(-1 * pdgid)
    if anti_id in special_names:
        antiName = special_names[anti_id]
    elif particle.PDGID(anti_id).is_valid and part.pdgid.has_fundamental_anti:
        anti_parti = Particle.from_pdgid(anti_id)
        antiName = c_identifier_camel(anti_parti.programmatic_name, abs(int(pdgid)) < 6000)
    else:
        antiName = c_id

    mass = 0.0 if part.mass is None else part.mass

    particle_db[c_id] = {
        "name": c_id,
        "antiName": antiName,
        "pdg": int(pdgid),
        "mass": mass,
        "charge": int(part.charge),
        "lifetime": part.lifetime,
        "ngc_code": counter,
        "isHadron": part.pdgid.is_hadron,
        "isNucleus": False,
    }

    counter += 1

nuclear_db = {}

for i in range(len(nuclear_data)):
    pdg, name, A, Z = nuclear_data[i]
    c_id = c_identifier_camel(name)

    L = 0
    I = 0

    nuclear_db[c_id] = {
        "name": name,
        "antiName": "Anti" + c_id,
        "pdg": int(pdg),
        "mass": Particle.from_pdgid(pdg).mass,
        "charge": Particle.from_pdgid(pdg).charge,
        "lifetime": Particle.from_pdgid(pdg).lifetime,
        "ngc_code": int(f"10{L:01d}{int(Z):03d}{int(A):03d}{I:01d}"),
        "A": int(A),
        "Z": int(Z),
        "isNucleus": True,
        "isHadron": True,
    }


######################################################
### Create the output file
######################################################

with open("GeneratedParticleProperties.inc", "w") as file:

    print(
        f"// generated by {__file__}\n"
        "// MANUAL EDITS ON OWN RISK. THEY WILL BE OVERWRITTEN. \n"
        "\n"
        "#pragma once"
        "\n"
        "namespace corsika {\n"
        "/** @ingroup Particles \n"
        "    @{ \n"
        "  */ \n",
        file=file,
    )

    ###############################
    ### Code values (enums)
    ###############################

    print("//! @cond EXCLUDE_DOXY", file=file)
    print("enum class Code : CodeIntType {", file=file)
    print('  FirstParticle = 1, // if you want to loop over particles, you want to start with "1"', file=file)
    for key, val in particle_db.items():
        print(f"  {key} = {val['ngc_code']},", file=file)

    print("  Nucleus = 1000000000,", file=file)
    for key, vals in nuclear_db.items():
        print(f"  {key} = {int(vals['pdg'])},", file=file)

    print("}; //! @endcond\n", file=file)

    ###############################
    ### PDGCode values
    ###############################

    print("//! @cond EXCLUDE_DOXY", file=file)
    print("enum class PDGCode : PDGCodeIntType {", file=file)
    for key, val in particle_db.items():
        print(f"  {key} = {int(val['pdg'])},", file=file)
    print("}; //! @endcond\n", file=file)

    ###############################
    ### Code names
    ###############################

    print("/** @cond EXCLUDE_DOXY */\n", file=file)
    print("namespace particle::detail {", file=file)

    print(f"static constexpr std::size_t size = {len(particle_db)};\n", file=file)

    print("constexpr std::initializer_list<Code> all_particles = {", file=file)
    for key, val in particle_db.items():
        print(f"  Code::{key},", file=file)
    print("};\n", file=file)

    ###############################
    ### Masses
    ###############################

    print("static constexpr std::array<HEPMassType const, size> masses = {", file=file)
    for key, val in particle_db.items():
        print(f"  {val['mass'] / hepunits.GeV:.6e} * 1e9 * electronvolt, // {key}", file=file)
    print("};\n", file=file)

    ###############################
    ### Propagation thresholds
    ###############################

    print("static std::array<HEPEnergyType, size> propagation_thresholds = {", file=file)
    for key, val in particle_db.items():
        print(f"  1e9 * electronvolt, // {key}", file=file)
    print("};\n", file=file)

    ###############################
    ### Production thresholds
    ###############################

    print("static HEPEnergyType threshold_nuclei = 0_eV;", file=file)
    print("static std::array<HEPEnergyType, size> production_thresholds = {", file=file)
    for key, val in particle_db.items():
        print(f"  1e6 * electronvolt, // {key}", file=file)
    print("};\n", file=file)

    ###############################
    ### PDG Names
    ###############################

    print("static constexpr std::array<PDGCode, size> pdg_codes = {", file=file)
    for key, val in particle_db.items():
        print(f"  PDGCode::{key},", file=file)
    print("};\n", file=file)

    ##################################
    ### Short, computer friendly names
    ##################################

    print("static constexpr std::array<std::string_view, size> names = {", file=file)
    for key, val in particle_db.items():
        print(f"  \"{val['name']}\",", file=file)
    print("};\n", file=file)

    ###############################
    ### Charges
    ###############################

    print("static constexpr std::array<int32_t, size> electric_charges = {", file=file)
    for key, val in particle_db.items():
        print(f"  {val['charge']}, // {key}", file=file)
    print("};\n", file=file)

    ###############################
    ### Lifetimes
    ###############################

    print("static constexpr std::array<double const, size> lifetime = {", file=file)
    for key, val in particle_db.items():
        time = val["lifetime"]
        if time == 0.0 or time is None:
            time = 3.335641e-27
        elif np.isinf(time):
            time = "std::numeric_limits<double>::infinity()"
        else:
            time = f"{time / hepunits.second:.6e}"

        print(f"  {time}, // {key}", file=file)
    print("};\n", file=file)

    ###############################
    ### Is hadron bool
    ###############################

    print("static constexpr std::array<bool, size> isHadron = {", file=file)
    for key, val in particle_db.items():
        print(f"  {str(val['isHadron']).lower()}, // {key}", file=file)
    print("};\n", file=file)

    ####################################
    ###### Conversion array
    ####################################

    # For all known particles, add them to a conversion table list
    # such that they are at pdg + N where N is `max_conversion_array_id`
    # but will be updated so that the conversion has the minimum size
    conversion_table = [None] * (max_conversion_array_id * 2 + 1)
    for key, val in particle_db.items():
        pdg = val["pdg"]
        ngc = val["ngc_code"]

        if abs(pdg) <= max_conversion_array_id:
            if conversion_table[int(pdg + max_conversion_array_id)]:
                raise Exception(f"table entry {pdg} already occupied")
            else:
                conversion_table[int(pdg + max_conversion_array_id)] = key

    first_index = 0
    while conversion_table[first_index] is None:
        first_index += 1

    conversion_table = conversion_table[first_index : len(conversion_table) - first_index]
    max_conversion_array_id = (len(conversion_table) - 1) // 2

    print("////////////////////////  pdg --> ngc conversion  /////////////////////////", file=file)
    print("// This is a lookup table for which you can get the C8 particle enum coding", file=file)
    print(f"// by using this_array[pdg + {max_conversion_array_id}]", file=file)
    print("// this is inteded for faster lookup than using the map, which is", file=file)
    print(f"// used for |pdg| values larger than {max_conversion_array_id}, instead\n", file=file)

    print(f"static std::array<Code, {int(max_conversion_array_id * 2 + 1)}> constexpr conversionArray {{", file=file)

    for val in conversion_table:
        print(f"  Code::{val if val else 'Unknown'},", file=file)
    print("};\n", file=file)

    ###############################
    ### Make the conversion map
    ###############################

    print("//////////////////////  pdg --> ngc conversion  //////////////////////", file=file)
    print("// This is a lookup table to convert from pdg to the C8 interal coding", file=file)
    print(f"// for |pdg| values larger than {max_conversion_array_id}\n", file=file)

    print("static std::map<PDGCode, Code> const conversionMap {", file=file)
    for key, val in particle_db.items():
        if abs(val["pdg"]) > max_conversion_array_id:
            print(f"  {{PDGCode::{key}, Code::{key}}},", file=file)
    print("};\n", file=file)

    print("}//end namespace particle::detail", file=file)
    print("/** @endcond */", file=file)
    print("/** @} */", file=file)
    print("} // end namespace corsika", file=file)


with open("GeneratedParticleClasses.inc", "w") as file:
    print(
        f"// generated by {__file__}\n"
        "// MANUAL EDITS ON OWN RISK. THEY WILL BE OVERWRITTEN. \n"
        "\n"
        "#pragma once"
        "\n"
        "namespace corsika {\n"
        "/** @ingroup Particles \n"
        "    @{ \n"
        "  */ \n\n"
        "// list of C++ classes to access particle properties\n"
        "/** @defgroup ParticleClasses \n"
        "    @{ */\n"
        "\n",
        file=file,
    )

    for key, val in particle_db.items():
        time = val["lifetime"]
        if time == 0.0 or time is None:
            time = "3.335641e-27 s"
        elif np.isinf(time):
            time = "inf"
        else:
            time = f"{time / hepunits.second:.6e} s"
        print(
            f"/** @class {key}\n"
            "\n"
            " * Particle properties are taken from PDG database:<br>\n"
            f" *  - pdg={int(val['pdg'])}\n"
            f" *  - mass={val['mass'] / hepunits.GeV:.6f} GeV \n"
            f" *  - charge={int(val['charge'])} \n"
            f" *  - lifetime={time}\n"
            f" *  - name={key}\n"
            f" *  - anti={val['antiName']}\n"
            "*/\n",
            file=file,
        )

        print(
            f"class {key} {{\n"
            "  /** @cond EXCLUDE_DOXY */\n"
            "  public:\n"
            f"    {key}() = delete;\n"
            f"    static constexpr Code code{{Code::{key}}};\n"
            f"    static constexpr Code anti_code{{Code::{val['antiName']}}};\n"
            "    static constexpr HEPMassType mass{corsika::get_mass(code)};\n"
            "    static constexpr ElectricChargeType charge{corsika::get_charge(code)};\n"
            "    static constexpr int charge_number{corsika::get_charge_number(code)};\n"
            "    static constexpr std::string_view name{corsika::get_name(code)};\n"
            "    static constexpr bool is_nucleus{corsika::is_nucleus(code)};\n"
            "  private:\n"
            "    static constexpr CodeIntType TypeIndex = static_cast<CodeIntType>(code);\n"
            "  /** @endcond */\n"
            "};\n",
            file=file,
        )

    for key, val in nuclear_db.items():
        print(
            f"/** @class {key}\n"
            "\n"
            " * Particle properties are taken from PDG database:<br>\n"
            f" *  - pdg={int(vals['pdg'])}\n"
            f" *  - mass={vals['mass'] / hepunits.GeV:.6f} GeV \n"
            f" *  - charge={int(vals['charge'])} \n"
            f" *  - name={key}\n"
            f" *  - anti=Unknown\n"
            f" *  - nuclear A={int(vals['A'])}\n"
            f" *  - nuclear Z={int(vals['Z'])}\n"
            "*/\n",
            file=file,
        )

        print(
            f"class {key} {{\n"
            "  /** @cond EXCLUDE_DOXY */\n"
            "  public:\n"
            f"    {key}() = delete;\n"
            f"    static constexpr Code code{{Code::{key}}};\n"
            "    static constexpr Code anti_code{Code::Unknown};\n"
            "    static constexpr HEPMassType mass{corsika::get_mass(code)};\n"
            "    static constexpr ElectricChargeType charge{corsika::get_charge(code)};\n"
            "    static constexpr int charge_number{corsika::get_charge_number(code)};\n"
            "    static constexpr std::string_view name{corsika::get_name(code)};\n"
            "    static constexpr bool is_nucleus{corsika::is_nucleus(code)};\n"
            "    static constexpr int nucleus_A{corsika::get_nucleus_A(code)};\n"
            "    static constexpr int nucleus_Z{corsika::get_nucleus_Z(code)};\n"
            "  private:\n"
            "    static constexpr CodeIntType TypeIndex = static_cast<CodeIntType>(code);\n"
            "  /** @endcond */\n"
            "};\n",
            file=file,
        )

    print("  //! @}", file=file)
    print("", file=file)
    print("/** @} */", file=file)
    print("} // end namespace corsika", file=file)

with open("particle_db.pkl", "wb") as file:
    common_db = OrderedDict(list(particle_db.items()) + list(nuclear_db.items()))
    pickle.dump(common_db, file)
