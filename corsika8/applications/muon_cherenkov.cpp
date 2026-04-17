/*
 * (c) Copyright 2018 CORSIKA Project, corsika-project@lists.kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

/* clang-format off */
// InteractionCounter used boost/histogram, which
// fails if boost/type_traits have been included before. Thus, we have
// to include it first...
#include <corsika/framework/process/InteractionCounter.hpp>
/* clang-format on */
#include <corsika/framework/core/Cascade.hpp>
#include <corsika/framework/core/EnergyMomentumOperations.hpp>
#include <corsika/framework/core/Logging.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/PhysicalGeometry.hpp>
#include <corsika/framework/geometry/Plane.hpp>
#include <corsika/framework/geometry/Sphere.hpp>
#include <corsika/framework/process/DynamicInteractionProcess.hpp>
#include <corsika/framework/process/ProcessSequence.hpp>
#include <corsika/framework/process/SwitchProcessSequence.hpp>
#include <corsika/framework/random/RNGManager.hpp>
#include <corsika/framework/random/PowerLawDistribution.hpp>
#include <corsika/framework/utility/CorsikaFenv.hpp>
#include <corsika/framework/utility/SaveBoostHistogram.hpp>

#include <corsika/modules/writers/EnergyLossWriter.hpp>
#include <corsika/modules/writers/InteractionWriter.hpp>
#include <corsika/modules/writers/LongitudinalWriter.hpp>
#include <corsika/modules/writers/ProductionWriter.hpp>
#include <corsika/modules/writers/PrimaryWriter.hpp>
#include <corsika/modules/writers/SubWriter.hpp>
#include <corsika/output/OutputManager.hpp>

#include <corsika/media/CORSIKA7Atmospheres.hpp>
#include <corsika/media/Environment.hpp>
#include <corsika/media/GeomagneticModel.hpp>
#include <corsika/media/GladstoneDaleRefractiveIndex.hpp>
#include <corsika/media/HomogeneousMedium.hpp>
#include <corsika/media/IMagneticFieldModel.hpp>
#include <corsika/media/LayeredSphericalAtmosphereBuilder.hpp>
#include <corsika/media/MediumPropertyModel.hpp>
#include <corsika/media/NuclearComposition.hpp>
#include <corsika/media/ShowerAxis.hpp>
#include <corsika/media/UniformMagneticField.hpp>

#include <corsika/modules/BetheBlochPDG.hpp>
#include <corsika/modules/Epos.hpp>
#include <corsika/modules/EposLhcr.hpp>
#include <corsika/modules/ObservationPlane.hpp>
#include <corsika/modules/PROPOSAL.hpp>
#include <corsika/modules/ParticleCut.hpp>
#include <corsika/modules/Pythia8.hpp>
#include <corsika/modules/QGSJetII.hpp>
#include <corsika/modules/QGSJetIII.hpp>
#include <corsika/modules/Sibyll.hpp>
#include <corsika/modules/Sophia.hpp>
#include <corsika/modules/StackInspector.hpp>
#include <corsika/modules/thinning/EMThinning.hpp>
#include <corsika/modules/LongitudinalProfile.hpp>
#include <corsika/modules/ProductionProfile.hpp>

// for ICRC2023
#ifdef WITH_FLUKA
#include <corsika/modules/FLUKA.hpp>
#else
#include <corsika/modules/UrQMD.hpp>
#endif
#include <corsika/modules/TAUOLA.hpp>

#include <corsika/modules/radio/CoREAS.hpp>
#include <corsika/modules/radio/RadioProcess.hpp>
#include <corsika/modules/radio/ZHS.hpp>
#include <corsika/modules/radio/observers/Observer.hpp>
#include <corsika/modules/radio/observers/TimeDomainObserver.hpp>
#include <corsika/modules/radio/detectors/ObserverCollection.hpp>
#include <corsika/modules/radio/propagators/TabulatedFlatAtmospherePropagator.hpp>

#include <corsika/setup/SetupStack.hpp>
#include <corsika/setup/SetupTrajectory.hpp>
#include <corsika/setup/SetupC7trackedParticles.hpp>

#include <boost/filesystem.hpp>

#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>

#include <corsika/modules/cherenkov/Cherenkov.hpp>
#include <fstream>
#include <sstream>

#include <chrono>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <limits>
#include <string>

// Eventio format support for Cherenkov data
namespace eventio_writer {
  // Sync marker for eventio blocks
  static const std::array<uint8_t, 4> SYNC_MARKER = {0x37, 0x8a, 0x1f, 0xd4};

  // IACT type codes  
  enum class IACTType : uint32_t {
    RunHeader = 1200,
    InputCard = 1212,
    AtmosphericProfile = 1216,
    TelescopeDefinition = 1201,
    EventHeader = 1202,
    ArrayOffsets = 1203,
    Longitudinal = 1211,
    TelescopeData = 1204,
    Photons = 1205,
    CameraLayout = 1206,
    TriggerTime = 1207,
    PhotoElectrons = 1208,
    EventEnd = 1209,
    RunEnd = 1210,   
  };

  // Helper to write little-endian integers
  inline void write_int32(std::ofstream& file, uint32_t val) {
    file.write(reinterpret_cast<const char*>(&val), 4);
  }

  inline void write_int16(std::ofstream& file, int16_t val) {
    file.write(reinterpret_cast<const char*>(&val), 2);
  }

  inline void write_float32(std::ofstream& file, float val) {
    file.write(reinterpret_cast<const char*>(&val), 4);
  }

  inline void write_bytes(std::ofstream& file, const std::vector<uint8_t>& data) {
    file.write(reinterpret_cast<const char*>(data.data()), data.size());
  }

  // Convert a 4-char string (e.g. "RUNH") to a float via memcpy
  inline float str_to_float(const char* s) {
    float result;
    std::memcpy(&result, s, 4);
    return result;
  }

  // Write eventio toplevel object header:
  //   4 bytes: sync marker
  //   4 bytes: type_word = type | (user << 16) | (extended << 17) | (version << 20)
  //   4 bytes: id (int32)
  //   4 bytes: length_word = content_size_bytes | (only_subobjects << 30)
  //   [4 bytes: extension_word if extended]
  // Total: 16 bytes (or 20 if extended)
  inline void write_toplevel_header(std::ofstream& file, IACTType type, int32_t id,
                                    uint64_t content_size_bytes,
                                    bool only_subobjects = false,
                                    uint32_t version = 0) {
    file.write(reinterpret_cast<const char*>(SYNC_MARKER.data()), 4);
    
    // Check if extended format is needed (content > 30 bits)
    bool extended = (content_size_bytes > 0x3FFFFFFF);
    
    uint32_t type_word = static_cast<uint32_t>(type) | (version << 20);
    if (extended) {
      type_word |= (1U << 17);  // Set extended bit (bit 17 per eventio spec)
    }
    write_int32(file, type_word);
    write_int32(file, static_cast<uint32_t>(id));
    
    uint32_t length_word = static_cast<uint32_t>(content_size_bytes & 0x3FFFFFFF);
    if (only_subobjects) {
      length_word |= (1U << 30);
    }
    write_int32(file, length_word);
    
    if (extended) {
      // Write upper bits of content_size
      uint32_t extension_word = static_cast<uint32_t>(content_size_bytes >> 30);
      write_int32(file, extension_word);
    }
  }

  // Write eventio sub-object header (no sync marker):
  //   4 bytes: type_word = type | (user << 16) | (extended << 17) | (version << 20)
  //   4 bytes: id (int32)
  //   4 bytes: length_word
  //   [4 bytes: extension_word if extended]
  // Total: 12 bytes (or 16 if extended)
  inline void write_subobject_header(std::ofstream& file, IACTType type, int32_t id,
                                     uint64_t content_size_bytes,
                                     bool only_subobjects = false,
                                     uint32_t version = 0) {
    // Check if extended format is needed (content > 30 bits)
    bool extended = (content_size_bytes > 0x3FFFFFFF);
    
    uint32_t type_word = static_cast<uint32_t>(type) | (version << 20);
    if (extended) {
      type_word |= (1U << 17);  // Set extended bit (bit 17 per eventio spec)
    }
    write_int32(file, type_word);
    write_int32(file, static_cast<uint32_t>(id));
    
    uint32_t length_word = static_cast<uint32_t>(content_size_bytes & 0x3FFFFFFF);
    if (only_subobjects) {
      length_word |= (1U << 30);
    }
    write_int32(file, length_word);
    
    if (extended) {
      // Write upper bits of content_size
      uint32_t extension_word = static_cast<uint32_t>(content_size_bytes >> 30);
      write_int32(file, extension_word);
    }
  }

  // PhotonBunch structure matching CHASM format
  struct PhotonBunch {
    float x;              // impact x (cm)
    float y;              // impact y (cm)
    float cx;             // direction cosine x
    float cy;             // direction cosine y
    float time;           // arrival time
    float zem;            // emission height (cm)
    float n_photons;      // number of photons
    float wavelength;     // wavelength (nm)
  };

  // Write a single photon bunch
  inline void write_photon_bunch(std::ofstream& file, const PhotonBunch& bunch) {
    write_float32(file, bunch.x);
    write_float32(file, bunch.y);
    write_float32(file, bunch.cx);
    write_float32(file, bunch.cy);
    write_float32(file, bunch.time);
    write_float32(file, bunch.zem);
    write_float32(file, bunch.n_photons);
    write_float32(file, bunch.wavelength);
  }

  // RunHeader block - matches CORSIKA 7 / eventio reference format
  // Reference: type=1200, id=102, content_size=1096 bytes
  class RunHeaderBlock {
  public:
    static void write(std::ofstream& file, int run_no = 1) {
      // content: 4 (user_len) + 4 (count) + 273*4 (floats) = 1100 bytes
      // But reference uses content_size=1096, so: 4 (count) + 273*4 = 1096
      const uint32_t content_size = 1096;  // matches reference
      write_toplevel_header(file, IACTType::RunHeader, /*id=*/102, content_size);
      
      // Count of floats
      write_int32(file, 273);
      
      // Word 1: RUNH marker
      write_float32(file, str_to_float("RUNH"));
      
      // Word 2: Run number
      write_float32(file, static_cast<float>(run_no));
      
      // Word 3: Date (placeholder)
      write_float32(file, 1.0f);
      
      // Words 4-273: Fill with zeros (273-3 = 270 more floats)
      for (int i = 3; i < 273; ++i) {
        write_float32(file, 0.0f);
      }
    }
  };

  // InputCard block - matches reference: type=1212, id=0, content_size=2748
  // But eventio only needs n_strings=0 (4 bytes) to parse successfully
  // Reference has 2748 bytes of content but we can use just 4
  class InputCardBlock {
  public:
    static void write(std::ofstream& file) {
      const uint32_t content_size = 4;  // just n_strings field
      write_toplevel_header(file, IACTType::InputCard, /*id=*/0, content_size);
      write_int32(file, 0);  // n_strings = 0
    }
  };

  // EventHeader block - matches reference: type=1202, id=1, content_size=1096
  class EventHeaderBlock {
  public:
    static void write(std::ofstream& file, double zenith_rad, double azimuth_rad,
                     double primary_energy, int event_no = 1) {
      const uint32_t content_size = 1096;
      write_toplevel_header(file, IACTType::EventHeader, /*id=*/event_no, content_size);
      
      // Count of floats
      write_int32(file, 273);
      
      // Word 1: EVTH marker
      write_float32(file, str_to_float("EVTH"));
      
      // Word 2: Event number
      write_float32(file, static_cast<float>(event_no));
      
      // Word 3: Particle ID
      write_float32(file, 6.0f);  // muon
      
      // Word 4: Total energy (GeV)
      write_float32(file, static_cast<float>(primary_energy));
      
      // Word 5-6: zeros (placeholder)
      write_float32(file, 0.0f);
      write_float32(file, 0.0f);
      
      // Word 7: first interaction height
      write_float32(file, 0.0f);
      
      // Word 8-10: momentum
      write_float32(file, 0.0f);
      write_float32(file, 0.0f);
      write_float32(file, 0.0f);
      
      // Word 11: zenith (rad)
      write_float32(file, static_cast<float>(zenith_rad));
      
      // Word 12: azimuth (rad)
      write_float32(file, static_cast<float>(azimuth_rad));
      
      // Words 13-273: Fill with zeros
      for (int i = 12; i < 273; ++i) {
        write_float32(file, 0.0f);
      }
    }
  };

  // TelescopeDefinition block - reference: type=1201, id=0, content_size=20
  // 4 bytes ntel + 4*4 bytes (x,y,z,r) per telescope
  class TelescopeDefBlock {
  public:
    static void write(std::ofstream& file, double tel_x, double tel_y, double tel_z,
                     double tel_radius, int ntel = 1) {
      const uint32_t content_size = 4 + ntel * 16;  // 4 (ntel) + ntel * 4 floats
      write_toplevel_header(file, IACTType::TelescopeDefinition, /*id=*/0, content_size);
      write_int32(file, static_cast<uint32_t>(ntel));
      write_float32(file, static_cast<float>(tel_x));
      write_float32(file, static_cast<float>(tel_y));
      write_float32(file, static_cast<float>(tel_z));
      write_float32(file, static_cast<float>(tel_radius));
    }
  };

  // Container for all photon bunches from one event
  struct CherenkovEvent {
    std::vector<PhotonBunch> bunches;
    double zenith_rad;
    double azimuth_rad;
    double primary_energy;
    int primary_pdg;
    double tel_x, tel_y, tel_z, tel_radius;
    int event_no;
    int run_no;
    
    CherenkovEvent() : event_no(1), run_no(1) {}
  };

  // Main eventio writer for Cherenkov data
  class CherenkovEventioWriter {
  private:
    CherenkovEvent event;
    std::string filename;

  public:
    CherenkovEventioWriter(const std::string& fname) : filename(fname) {}

    void add_photon(float x, float y, float cx, float cy, float time, 
                   float zem, float n_photons, float wavelength) {
      PhotonBunch bunch{x, y, cx, cy, time, zem, n_photons, wavelength};
      event.bunches.push_back(bunch);
    }

    void set_event_params(double zen, double az, double energy, int pdg,
                         double tx, double ty, double tz, double tr) {
      event.zenith_rad = zen;
      event.azimuth_rad = az;
      event.primary_energy = energy;
      event.primary_pdg = pdg;
      event.tel_x = tx;
      event.tel_y = ty;
      event.tel_z = tz;
      event.tel_radius = tr;
    }

    void write() {
      std::ofstream file(filename, std::ios::binary);
      if (!file) {
        throw std::runtime_error("Cannot open eventio file: " + filename);
      }

      // Write RunHeader (CORSIKA 7 format)
      RunHeaderBlock::write(file, event.run_no);

      // Write InputCard (required by eventio format)
      InputCardBlock::write(file);

      // Write TelescopeDefinition (CHASM-compatible)
      TelescopeDefBlock::write(file, event.tel_x, event.tel_y, event.tel_z,
                              event.tel_radius, 1);

      // Write EventHeader with shower parameters
      EventHeaderBlock::write(file, event.zenith_rad, event.azimuth_rad,
                             event.primary_energy, event.event_no);

      // Write ArrayOffsets (required by eventio between EventHeader and TelescopeData)
      // Format: n_arrays(int32) + time_offset(float) + x[n_arrays](float) + y[n_arrays](float)
      // Version 0: columns are (x, y) per array
      {
        const int n_arrays = 1;  // one reuse
        const uint32_t content_size = 4 + 4 + n_arrays * 4 + n_arrays * 4;  // 16 bytes
        write_toplevel_header(file, IACTType::ArrayOffsets, /*id=*/0, content_size);
        write_int32(file, n_arrays);
        write_float32(file, 0.0f);  // time_offset
        write_float32(file, 0.0f);  // x offset for array 0
        write_float32(file, 0.0f);  // y offset for array 0
      }

      // Write TelescopeData containing the photon bunches
      // TelescopeData is ALWAYS required (one per reuse in ArrayOffsets)
      {
        if (!event.bunches.empty()) {
          // Photons sub-object content: arr(2) + tel(2) + n_photons(4) + n_bunches(4) + bunches
          const uint64_t n_bunches = event.bunches.size();
          const uint64_t photons_content_size = 4 + 4 + 4 + n_bunches * 8 * 4;
          
          // Subobject header is 12 bytes normally, 16 bytes if extended (content > 30 bits)
          const uint64_t photons_header_size = (photons_content_size > 0x3FFFFFFF) ? 16 : 12;
          
          // TelescopeData content = Photons sub-object header + Photons content
          const uint64_t telescope_content_size = photons_header_size + photons_content_size;

          // TelescopeData: toplevel, only_subobjects=true
          write_toplevel_header(file, IACTType::TelescopeData, /*id=*/0,
                                telescope_content_size, /*only_subobjects=*/true);

          // Photons: sub-object (no sync marker)
          write_subobject_header(file, IACTType::Photons, /*id=*/0, photons_content_size);

          // Photons content: arr(int16) + tel(int16) + n_photons(float) + n_bunches(int32)
          write_int16(file, 0);  // array number
          write_int16(file, 0);  // telescope number
          float total_photons = 0;
          for (const auto& b : event.bunches) total_photons += b.n_photons;
          write_float32(file, total_photons);
          write_int32(file, static_cast<uint32_t>(n_bunches));

          // Write all photon bunches (8 floats each)
          for (const auto& bunch : event.bunches) {
            write_photon_bunch(file, bunch);
          }
        } else {
          // Empty TelescopeData (no sub-objects, content_size=0)
          write_toplevel_header(file, IACTType::TelescopeData, /*id=*/0,
                                0, /*only_subobjects=*/true);
        }
      }

      // Write EventEnd: type=1209, id=event_no, content_size=1096
      // Same structure as EventHeader: count(273) + 273 floats
      {
        const uint32_t content_size = 1096;
        write_toplevel_header(file, IACTType::EventEnd, /*id=*/event.event_no, content_size);
        write_int32(file, 273);
        write_float32(file, str_to_float("EVTE"));
        write_float32(file, static_cast<float>(event.event_no));
        for (int i = 2; i < 273; ++i) {
          write_float32(file, 0.0f);
        }
      }

      // Write RunEnd: type=1210, id=102, content_size=16
      // content: count(3) + RUNE + run_no + n_events
      {
        const uint32_t content_size = 16;
        write_toplevel_header(file, IACTType::RunEnd, /*id=*/102, content_size);
        write_int32(file, 3);
        write_float32(file, str_to_float("RUNE"));
        write_float32(file, static_cast<float>(event.run_no));
        write_float32(file, 1.0f);  // n_events
      }
    }
  };
}

using namespace corsika;
using namespace std;

using EnvironmentInterface =
    IRefractiveIndexModel<IMediumPropertyModel<IMagneticFieldModel<IMediumModel>>>;
using EnvType = Environment<EnvironmentInterface>;
using StackType = setup::Stack<EnvType>;
using TrackingType = setup::Tracking;
using Particle = StackType::particle_type;

//
// This is the main example script which runs EAS with fairly standard settings
// w.r.t. what was implemented in CORSIKA 7. Users may want to change some of the
// specifics (observation altitude, magnetic field, energy cuts, etc.), but this
// example is the most physics-complete one and should be used for full simulations
// of particle cascades in air
//

long registerRandomStreams(long seed) {
  RNGManager<>::getInstance().registerRandomStream("cascade");
  RNGManager<>::getInstance().registerRandomStream("qgsjet");
  RNGManager<>::getInstance().registerRandomStream("qgsjetIII");
  RNGManager<>::getInstance().registerRandomStream("sibyll");
  RNGManager<>::getInstance().registerRandomStream("sophia");
  RNGManager<>::getInstance().registerRandomStream("epos");
  RNGManager<>::getInstance().registerRandomStream("epos-lhcr");
  RNGManager<>::getInstance().registerRandomStream("pythia");
  RNGManager<>::getInstance().registerRandomStream("urqmd");
  RNGManager<>::getInstance().registerRandomStream("fluka");
  RNGManager<>::getInstance().registerRandomStream("proposal");
  RNGManager<>::getInstance().registerRandomStream("thinning");
  RNGManager<>::getInstance().registerRandomStream("tauola");
  RNGManager<>::getInstance().registerRandomStream("primary_particle");
  if (seed == 0) {
    std::random_device rd;
    seed = rd();
    CORSIKA_LOG_INFO("random seed (auto) {}", seed);
  } else {
    CORSIKA_LOG_INFO("random seed {}", seed);
  }
  RNGManager<>::getInstance().setSeed(seed);
  return seed;
}

template <typename T>
using MyExtraEnv =
    GladstoneDaleRefractiveIndex<MediumPropertyModel<UniformMagneticField<T>>>;

int main(int argc, char** argv) {

  auto const run_start = std::chrono::steady_clock::now();

  // the main command line description
  CLI::App app{"Simulate standard (downgoing) showers with CORSIKA 8."};

  CORSIKA_LOG_INFO(
      "Please cite the following papers when using CORSIKA 8:\n"
      " - \"Towards a Next Generation of CORSIKA: A Framework for the Simulation of "
      "Particle Cascades in Astroparticle Physics\", Comput. Softw. Big Sci. 3 (2019) "
      "2, https://doi.org/10.1007/s41781-018-0013-0\n"
      " - \"Simulating radio emission from particle cascades with CORSIKA 8\", "
      "Astropart. Phys. 166 (2025) 103072, "
      "https://doi.org/10.1016/j.astropartphys.2024.103072");

  //////// Primary options ////////

  // some options that we want to fill in
  int A, Z, nevent = 0;
  std::vector<double> cli_energy_range;

  // the following section adds the options to the parser

  // we start by definining a sub-group for the primary ID
  auto opt_Z = app.add_option("-Z", Z, "Atomic number for primary")
                   ->check(CLI::Range(0, 26))
                   ->group("Primary");
  auto opt_A = app.add_option("-A", A, "Atomic mass number for primary")
                   ->needs(opt_Z)
                   ->check(CLI::Range(1, 58))
                   ->group("Primary");
  app.add_option("-p,--pdg",
                 "PDG code for primary (p=2212, gamma=22, e-=11, nu_e=12, mu-=13, "
                 "nu_mu=14, tau=15, nu_tau=16).")
      ->excludes(opt_A)
      ->excludes(opt_Z)
      ->group("Primary");
  app.add_option("-E,--energy", "Primary energy in GeV")->default_val(0);
  app.add_option("--energy_range", cli_energy_range,
                 "Low and high values that define the range of the primary energy in GeV")
      ->expected(2)
      ->check(CLI::PositiveNumber)
      ->group("Primary");
  app.add_option("--eslope", "Spectral index for sampling energies, dN/dE = E^eSlope")
      ->default_val(-1.0)
      ->group("Primary");
  app.add_option("-z,--zenith", "Primary zenith angle (deg)")
      ->default_val(0.)
      ->check(CLI::Range(0., 92.))
      ->group("Primary");
  app.add_option("-a,--azimuth", "Primary azimuth angle (deg)")
      ->default_val(0.)
      ->check(CLI::Range(0., 360.))
      ->group("Primary");

  //////// Config options ////////

  app.add_option("--emcut",
                 "Min. kin. energy of photons, electrons and "
                 "positrons in tracking (GeV)")
      ->default_val(0.5e-3)
      ->check(CLI::Range(0.000001, 1.e13))
      ->group("Config");
  app.add_option("--hadcut", "Min. kin. energy of hadrons in tracking (GeV)")
      ->default_val(0.3)
      ->check(CLI::Range(0.02, 1.e13))
      ->group("Config");
  app.add_option("--mucut", "Min. kin. energy of muons in tracking (GeV)")
      ->default_val(0.3)
      ->check(CLI::Range(0.000001, 1.e13))
      ->group("Config");
  app.add_option("--taucut", "Min. kin. energy of tau leptons in tracking (GeV)")
      ->default_val(0.3)
      ->check(CLI::Range(0.000001, 1.e13))
      ->group("Config");
  app.add_option("--max-deflection-angle",
                 "maximal deflection angle in tracking in radians")
      ->default_val(0.2)
      ->check(CLI::Range(1.e-8, 1.))
      ->group("Config");
  bool track_neutrinos = false;
  app.add_flag("--track-neutrinos", track_neutrinos, "switch on tracking of neutrinos")
      ->group("Config");

  //////// Misc options ////////

  app.add_option("--neutrino-interaction-type",
                 "charged (CC) or neutral current (NC) or both")
      ->default_val("both")
      ->check(CLI::IsMember({"neutral", "NC", "charged", "CC", "both"}))
      ->group("Misc.");
  app.add_option("--observation-level",
                 "Height above earth radius of the observation level (in m)")
      ->default_val(0.)
      ->check(CLI::Range(-1.e3, 1.e5))
      ->group("Config");
  app.add_option("--injection-height",
                 "Height above earth radius of the injection point (in m)")
      ->default_val(112.75e3)
      ->check(CLI::Range(-1.e3, 1.e6))
      ->group("Config");
  app.add_option("-N,--nevent", nevent, "The number of events/showers to run.")
      ->default_val(1)
      ->check(CLI::PositiveNumber)
      ->group("Library/Output");
  app.add_option("-f,--filename", "Filename for output library.")
      ->required()
      ->default_val("corsika_library")
      ->check(CLI::NonexistentPath)
      ->group("Library/Output");
  bool compressOutput = false;
  app.add_flag("--compress", compressOutput, "Compress the output directory to a tarball")
      ->group("Library/Output");
  app.add_option("-s,--seed", "The random number seed.")
      ->default_val(0)
      ->check(CLI::NonNegativeNumber)
      ->group("Misc.");
  bool force_interaction = false;
  app.add_flag("--force-interaction", force_interaction,
               "Force the location of the first interaction.")
      ->group("Misc.");
  bool force_decay = false;
  app.add_flag("--force-decay", force_decay, "Force the primary to immediately decay")
      ->group("Misc.");
  bool disable_interaction_hists = false;
  app.add_flag("--disable-interaction-histograms", disable_interaction_hists,
               "Store interaction histograms")
      ->group("Misc.");
  app.add_option("-v,--verbosity", "Verbosity level: warn, info, debug, trace.")
      ->default_val("info")
      ->check(CLI::IsMember({"warn", "info", "debug", "trace"}))
      ->group("Misc.");
  app.add_option("-M,--hadronModel", "High-energy hadronic interaction model")
      ->default_val("SIBYLL-2.3d")
      ->check(CLI::IsMember(
          {"SIBYLL-2.3d", "QGSJet-II.04", "QGSJet-III", "EPOS-LHC-R", "Pythia8"}))
      ->group("Misc.");
  app.add_option("-T,--hadronModelTransitionEnergy",
                 "Transition between high-/low-energy hadronic interaction "
                 "model in GeV")
      ->default_val(std::pow(10, 1.9)) // 79.4 GeV
      ->check(CLI::NonNegativeNumber)
      ->group("Misc.");

  //////// Thinning options ////////

  app.add_option("--emthin",
                 "fraction of primary energy at which thinning of EM particles starts")
      ->default_val(1.e-6)
      ->check(CLI::Range(0., 1.))
      ->group("Thinning");
  app.add_option("--max-weight",
                 "maximum weight for thinning of EM particles (0 to select Kobal's "
                 "optimum times 0.5)")
      ->default_val(0)
      ->check(CLI::NonNegativeNumber)
      ->group("Thinning");
  bool multithin = false;
  app.add_flag("--multithin", multithin, "keep thinned particles (with weight=0)")
      ->group("Thinning");
  app.add_option("--ring", "concentric ring of star shape pattern of observers")
      ->default_val(0)
      ->check(CLI::Range(0, 20))
      ->group("Radio");
  app.add_option("--cherenkov-projection",
                 "Cherenkov photon projection mode: ground (fixed z-plane), shower (perpendicular to shower axis), or telescope (circular detector)")
      ->default_val("ground")
      ->check(CLI::IsMember({"ground", "shower", "telescope"}))
      ->group("Misc.");
  app.add_option("--cherenkov-output-format",
                 "Output format for Cherenkov data: csv or eventio (CHASM-compatible for GrOptics)")
      ->default_val("csv")
      ->check(CLI::IsMember({"csv", "eventio"}))
      ->group("Misc.");
  app.add_option("--telescope-x",
                 "Telescope surface x-coordinate in meters (relative to observation height at origin)")
      ->default_val(0.0)
      ->group("Misc.");
  app.add_option("--telescope-y",
                 "Telescope surface y-coordinate in meters (relative to observation height at origin)")
      ->default_val(0.0)
      ->group("Misc.");
  app.add_option("--telescope-z",
                 "Telescope surface z-coordinate in meters (relative to observation height at origin)")
      ->default_val(0.0)
      ->group("Misc.");
  app.add_option("--telescope-radius",
                 "Telescope surface radius in meters (for telescope projection mode)")
      ->default_val(100.0)
      ->check(CLI::PositiveNumber)
      ->group("Misc.");
  app.add_option("--telescope-azimuth",
                 "Telescope pointing azimuth angle in degrees (0=+x, 90=+y, 180=-x, 270=-y). If not specified, uses shower azimuth.")
      ->default_val(-1.0)
      ->group("Misc.");
  app.add_option("--telescope-zenith",
                 "Telescope pointing zenith angle in degrees (0=down, 90=horizontal, 180=up). If not specified, uses shower zenith.")
      ->default_val(-1.0)
      ->group("Misc.");
  app.add_option("--atmosphere-file",
                 "Path to atmospheric model file (columns: Alt[km], rho, thick, n-1)")
      ->default_val("/uufs/chpc.utah.edu/common/home/u1520754/corsika/modules/data/CHERENKOV/atmosphere/atmprof9.dat")
      ->group("Misc.");
  bool disable_dispersion = false;
  app.add_flag("--disable-dispersion", disable_dispersion,
               "Disable wavelength-dependent group-velocity dispersion for benchmarking")
      ->group("Misc.");
  bool enable_curvature = false;
  app.add_flag("--enable-curvature", enable_curvature,
               "Enable curved photon path correction for large zenith angles (>60°)")
      ->group("Misc.");
  app.add_option("--atmosphere-absorption-file",
                 "Path to atmospheric absorption table file for Cherenkov photons. "
                 "When specified, applies Monte Carlo absorption based on the table. "
                 "File format: text file with header 'n_wavelengths n_heights', "
                 "followed by wavelength array [nm], height array [m], and ecoeff matrix. "
                 "Use convert_npz_to_dat.py to convert CHASM abstable.npz files.")
      ->default_val("")
      ->group("Cherenkov");
  // parse the command line options into the variables
  CLI11_PARSE(app, argc, argv);

  if (app.count("--verbosity")) {
    auto const loglevel = app["--verbosity"]->as<std::string>();
    if (loglevel == "warn") {
      logging::set_level(logging::level::warn);
    } else if (loglevel == "info") {
      logging::set_level(logging::level::info);
    } else if (loglevel == "debug") {
      logging::set_level(logging::level::debug);
    } else if (loglevel == "trace") {
#ifndef _C8_DEBUG_
      CORSIKA_LOG_ERROR("trace log level requires a Debug build.");
      return 1;
#endif
      logging::set_level(logging::level::trace);
    }
  }

  // check that we got either PDG or A/Z
  // this can be done with option_groups but the ordering
  // gets all messed up
  if (app.count("--pdg") == 0) {
    if ((app.count("-A") == 0) || (app.count("-Z") == 0)) {
      CORSIKA_LOG_ERROR("If --pdg is not provided, then both -A and -Z are required.");
      return 1;
    }
  }

  // initialize random number sequence(s)
  auto seed = registerRandomStreams(app["--seed"]->as<long>());

  /* === START: SETUP ENVIRONMENT AND ROOT COORDINATE SYSTEM === */
  EnvType env;
  CoordinateSystemPtr const& rootCS = env.getCoordinateSystem();
  Point const center{rootCS, 0_m, 0_m, 0_m};
  Point const surface_{rootCS, 0_m, 0_m, constants::EarthRadius::Mean};
  GeomagneticModel wmm(center, corsika_data("GeoMag/WMM.COF"));

  // build an atmosphere with Keilhauer's parametrization of the
  // US standard atmosphere into `env`
    create_5layer_atmosphere<EnvironmentInterface, MyExtraEnv>(
      env, AtmosphereId::USStdBK, center, 1.000327, surface_, Medium::AirDry1Atm,
      // Set geomagnetic field: Bx=28.6 µT (horizontal), Bz=46.4 µT downward => -46.4 µT in +z-up convention
      MagneticFieldVector{rootCS, 28.6_uT, 0_uT, -46.4_uT});

  /* === END: SETUP ENVIRONMENT AND ROOT COORDINATE SYSTEM === */

  /* === START: CONSTRUCT PRIMARY PARTICLE === */

  // parse the primary ID as a PDG or A/Z code
  Code beamCode;

  // check if we want to use a PDG code instead
  if (app.count("--pdg") > 0) {
    beamCode = convert_from_PDG(PDGCode(app["--pdg"]->as<int>()));
  } else {
    // check manually for proton and neutrons
    if ((A == 1) && (Z == 1))
      beamCode = Code::Proton;
    else if ((A == 1) && (Z == 0))
      beamCode = Code::Neutron;
    else
      beamCode = get_nucleus_code(A, Z);
  }

  HEPEnergyType eMin = 0_GeV;
  HEPEnergyType eMax = 0_GeV;
  // check the particle energy parameters
  if (app["--energy"]->as<double>() > 0.0) {
    eMin = app["--energy"]->as<double>() * 1_GeV;
    eMax = app["--energy"]->as<double>() * 1_GeV;
  } else if (cli_energy_range.size()) {
    if (cli_energy_range[0] > cli_energy_range[1]) {
      CORSIKA_LOG_WARN(
          "Energy range lower bound is greater than upper bound. swapping...");
      eMin = cli_energy_range[1] * 1_GeV;
      eMax = cli_energy_range[0] * 1_GeV;
    } else {
      eMin = cli_energy_range[0] * 1_GeV;
      eMax = cli_energy_range[1] * 1_GeV;
    }
  } else {
    CORSIKA_LOG_CRITICAL(
        "Must set either the (--energy) flag or the (--energy_range) flag to "
        "positive value(s)");
    return 0;
  }

  // direction of the shower in (theta, phi) space
  auto const thetaRad = app["--zenith"]->as<double>() / 180. * M_PI;
  auto const phiRad = app["--azimuth"]->as<double>() / 180. * M_PI;

  auto const [nx, ny, nz] = std::make_tuple(sin(thetaRad) * cos(phiRad),
                                            sin(thetaRad) * sin(phiRad), -cos(thetaRad));
  auto propDir = DirectionVector(rootCS, {nx, ny, nz});
  /* === END: CONSTRUCT PRIMARY PARTICLE === */

  /* === START: CONSTRUCT GEOMETRY === */
  auto const observationHeight =
      app["--observation-level"]->as<double>() * 1_m + constants::EarthRadius::Mean;
  auto const injectionHeight =
      app["--injection-height"]->as<double>() * 1_m + constants::EarthRadius::Mean;
  auto const t = -observationHeight * cos(thetaRad) +
                 sqrt(-static_pow<2>(sin(thetaRad) * observationHeight) +
                      static_pow<2>(injectionHeight));
  Point const showerCore{rootCS, 0_m, 0_m, observationHeight};
  Point const injectionPos =
      showerCore + DirectionVector{rootCS,
                                   {-sin(thetaRad) * cos(phiRad),
                                    -sin(thetaRad) * sin(phiRad), cos(thetaRad)}} *
                       t;
  
  // we make the axis much longer than the inj-core distance since the
  // profile will go beyond the core, depending on zenith angle
  ShowerAxis const showerAxis{injectionPos, (showerCore - injectionPos) * 1.2, env};
  auto const dX = 10_g / square(1_cm); // Binning of the writers along the shower axis
  /* === END: CONSTRUCT GEOMETRY === */

// Cherenkov process setup
// Use observation plane as "ground" in rootCS:
// Pass absolute ground Z coordinate for projection calculations
  double const groundZ_m_abs = observationHeight / 1_m;  // absolute coordinate

// Base output path from -f/--filename
  auto const outBase = app["--filename"]->as<std::string>();
  boost::filesystem::path basePath(outBase);

// Directory where we want the Cherenkov output:
//   - if -f is "muon_test"          -> "."
//   - if -f is "/home/user/test_muon" -> "/home/user"
  boost::filesystem::path cherDir = basePath.parent_path();
  if (cherDir.empty()) {
    cherDir = ".";  // current working directory
  }

  boost::filesystem::create_directories(cherDir);

  // Get output format preference
  auto const cherOutputFormat = app["--cherenkov-output-format"]->as<std::string>();
  bool useEventio = (cherOutputFormat == "eventio");

  // Get primary particle info for eventio header
  int primary_pdg = 0;
  if (app.count("--pdg") > 0) {
    primary_pdg = app["--pdg"]->as<int>();
  } else {
    // Approximate PDG codes for nuclei (not exact but functional)
    primary_pdg = 1000 * Z + A;  // CORSIKA nucleus code as approximation
  }
  
  double primary_energy = (app["--energy"]->as<double>() > 0.0) 
    ? app["--energy"]->as<double>() 
    : (cli_energy_range.size() ? cli_energy_range[0] : 1.0);

  // Setup output based on format
  std::shared_ptr<std::ofstream> cherenkovCsvPtr;
  std::shared_ptr<eventio_writer::CherenkovEventioWriter> eventioWriterPtr;
  std::size_t cherenkovHitCount = 0;

  if (useEventio) {
    // Prepare eventio file (CHASM-compatible)
    boost::filesystem::path eventioFile = cherDir / "cherenkov_hits.dat";
    eventioWriterPtr = std::make_shared<eventio_writer::CherenkovEventioWriter>(
        eventioFile.string());
    
    // Set event parameters
    double tel_x = app["--telescope-x"]->as<double>();
    double tel_y = app["--telescope-y"]->as<double>();
    double tel_z_relative = app["--telescope-z"]->as<double>();
    double tel_z = observationHeight / 1_m + tel_z_relative;
    double tel_radius = app["--telescope-radius"]->as<double>();
    
    eventioWriterPtr->set_event_params(
        thetaRad, phiRad, primary_energy, primary_pdg,
        tel_x, tel_y, tel_z, tel_radius);
    
    CORSIKA_LOG_INFO("Writing Cherenkov hits to eventio file: {}", 
                     eventioFile.string());
  } else {
    // Prepare CSV file
    boost::filesystem::path cherFile = cherDir / "cherenkov_hits.csv";
    cherenkovCsvPtr = std::make_shared<std::ofstream>(cherFile.string());
    
    if (!(*cherenkovCsvPtr)) {
      CORSIKA_LOG_INFO("ERROR: could not open {}", cherFile.string());
    } else {
      CORSIKA_LOG_INFO("Writing Cherenkov hits to CSV: {}", cherFile.string());
      (*cherenkovCsvPtr) << "x_m,y_m,z_m,"
                          "wavelength_nm,z_production_m,"
                          "dir_x,dir_y,dir_z,"
                          "dir_ground_x,dir_ground_y,dir_ground_z,"
                          "time_ns,weight,n_minus_1,path_length_m\n";
    }
  }

// Define recorder lambda
  corsika::cherenkov::CherenkovProcess::HitRecorder cherenkovRecorder;
  if (useEventio && eventioWriterPtr) {
    cherenkovRecorder =
        [eventioWriterPtr, &cherenkovHitCount](
            double x_m, double y_m, double z_m,
            double lambda_nm, double z_prod_m,
            double dir_x, double dir_y, double dir_z,
            double dir_ground_x, double dir_ground_y, double dir_ground_z,
            double time_ns, double weight, double n_prod, double path_length_m) {
          ++cherenkovHitCount;

          (void)z_m;
          (void)dir_z;
          (void)dir_ground_x;
          (void)dir_ground_y;
          (void)dir_ground_z;
          (void)n_prod;
          (void)path_length_m;
          
          // Convert to eventio format (CHASM-compatible)
          // x, y in cm (multiply by 100)
          // zem in cm
          // time in ns
          // wavelength in nm
          // cx, cy use shower/projection coordinates (dir_x, dir_y)
          
          float x_cm = static_cast<float>(x_m * 100.0);
          float y_cm = static_cast<float>(y_m * 100.0);
          float cx = static_cast<float>(dir_x);
          float cy = static_cast<float>(dir_y);
          float time_ns_f = static_cast<float>(time_ns);
          float zem_cm = static_cast<float>(z_prod_m * 100.0);
          float n_photons_f = static_cast<float>(weight);
          float wavelength_f = static_cast<float>(lambda_nm);
          
          eventioWriterPtr->add_photon(x_cm, y_cm, cx, cy, time_ns_f,
                                       zem_cm, n_photons_f, wavelength_f);
        };
  } else if (cherenkovCsvPtr && *cherenkovCsvPtr) {
    cherenkovRecorder =
        [cherenkovCsvPtr, &cherenkovHitCount](
            double x_m, double y_m, double z_m,
            double lambda_nm, double z_prod_m,
            double dir_x, double dir_y, double dir_z,
            double dir_ground_x, double dir_ground_y, double dir_ground_z,
            double time_ns, double weight, double n_prod, double path_length_m) {
          ++cherenkovHitCount;
          (*cherenkovCsvPtr) << x_m << "," << y_m << "," << z_m << ","
                             << lambda_nm << "," << z_prod_m << ","
                             << dir_x << "," << dir_y << "," << dir_z << ","
                             << dir_ground_x << "," << dir_ground_y << "," 
                             << dir_ground_z << ","
                             << time_ns << "," << weight << "," 
                             << (n_prod - 1.0) << ","
                             << path_length_m << "\n";
        };
  }


// Instantiate Cherenkov process
  auto projModeStr = app["--cherenkov-projection"]->as<std::string>();
  corsika::cherenkov::ProjectionMode projMode;
  if (projModeStr == "shower") {
    projMode = corsika::cherenkov::ProjectionMode::ShowerAxis;
  } else if (projModeStr == "telescope") {
    projMode = corsika::cherenkov::ProjectionMode::Telescope;
  } else {
    projMode = corsika::cherenkov::ProjectionMode::Ground;
  }
  
  // Parse telescope options (z relative to observation level)
  double const tel_x = app["--telescope-x"]->as<double>();
  double const tel_y = app["--telescope-y"]->as<double>();
  double const tel_z_relative = app["--telescope-z"]->as<double>(); // relative to observation level
  double const tel_radius = app["--telescope-radius"]->as<double>();

  // Parse telescope pointing direction (azimuth and zenith in degrees)
  double tel_azimuth_deg = app["--telescope-azimuth"]->as<double>();
  double tel_zenith_deg = app["--telescope-zenith"]->as<double>();
  
  // Convert to radians (use shower direction if not specified, i.e., value is -1)
  double const tel_azimuth_rad = (tel_azimuth_deg >= 0.0) ? 
    (tel_azimuth_deg * M_PI / 180.0) : 
    phiRad;  // phiRad is the shower azimuth in radians
  double const tel_zenith_rad = (tel_zenith_deg >= 0.0) ? 
    (tel_zenith_deg * M_PI / 180.0) : 
    thetaRad;  // thetaRad is the shower zenith in radians

  // Convert telescope Z to absolute coordinates (same reference as shower core)
  double const tel_z = observationHeight / 1_m + tel_z_relative;
  
  // Shower core coordinates in observation plane reference frame
  double const shower_core_x = showerCore.getCoordinates().getX() / 1_m;
  double const shower_core_y = showerCore.getCoordinates().getY() / 1_m;
  double const shower_core_z = showerCore.getCoordinates().getZ() / 1_m;

  // Load atmosphere file and build n(z) function (z absolute in meters)
  // Using cubic spline interpolation on ln(n-1) for maximum smoothness
  std::string atmPath = app["--atmosphere-file"]->as<std::string>();
  std::ifstream atmFile(atmPath);
  std::vector<double> atmAlt_km;           // altitude in km
  std::vector<double> atmNminus1;          // n-1 values
  std::vector<double> atmLnNminus1;        // log(n-1) values for interpolation
  std::vector<double> atmSplineCoeffs;     // cubic spline coefficients (second derivatives)
  
  if (atmFile) {
    std::string line;
    while (std::getline(atmFile, line)) {
      // skip comments and empty lines
      if (line.empty() || line[0] == '#') continue;
      std::istringstream iss(line);
      double alt_km=0, rho=0, thick=0, nminus1=0;
      if (iss >> alt_km >> rho >> thick >> nminus1) {
        atmAlt_km.push_back(alt_km);
        atmNminus1.push_back(nminus1);
        // Pre-compute ln(n-1) for cubic spline interpolation
        // Handle edge case: very small (n-1) at high altitude
        if (nminus1 > 1e-10) {
          atmLnNminus1.push_back(std::log(nminus1));
        } else {
          // For very high altitudes, treat as vacuum (n ≈ 1)
          atmLnNminus1.push_back(std::log(1e-10));
        }
      }
    }
    // Data should already be sorted, but ensure it
    if (!atmAlt_km.empty()) {
      std::vector<std::size_t> indices(atmAlt_km.size());
      std::iota(indices.begin(), indices.end(), 0);
      std::sort(indices.begin(), indices.end(),
                [&](std::size_t a, std::size_t b) { return atmAlt_km[a] < atmAlt_km[b]; });
      
      std::vector<double> tmp_alt(atmAlt_km.size());
      std::vector<double> tmp_nminus1(atmNminus1.size());
      std::vector<double> tmp_lnNminus1(atmLnNminus1.size());
      for (std::size_t i = 0; i < indices.size(); ++i) {
        tmp_alt[i] = atmAlt_km[indices[i]];
        tmp_nminus1[i] = atmNminus1[indices[i]];
        tmp_lnNminus1[i] = atmLnNminus1[indices[i]];
      }
      atmAlt_km = tmp_alt;
      atmNminus1 = tmp_nminus1;
      atmLnNminus1 = tmp_lnNminus1;
      
      // Compute natural cubic spline coefficients (second derivatives) for ln(n-1)
      // This provides C² continuity for smoother interpolation
      std::size_t n = atmAlt_km.size();
      if (n >= 3) {
        atmSplineCoeffs.resize(n, 0.0);
        std::vector<double> h(n-1), alpha(n-1), l(n), mu(n), z(n);
        
        // Compute intervals
        for (std::size_t i = 0; i < n-1; ++i) {
          h[i] = atmAlt_km[i+1] - atmAlt_km[i];
        }
        
        // Compute alpha values
        for (std::size_t i = 1; i < n-1; ++i) {
          alpha[i] = (3.0 / h[i]) * (atmLnNminus1[i+1] - atmLnNminus1[i]) -
                     (3.0 / h[i-1]) * (atmLnNminus1[i] - atmLnNminus1[i-1]);
        }
        
        // Solve tridiagonal system (natural spline: zero second derivative at boundaries)
        l[0] = 1.0;
        mu[0] = 0.0;
        z[0] = 0.0;
        
        for (std::size_t i = 1; i < n-1; ++i) {
          l[i] = 2.0 * (atmAlt_km[i+1] - atmAlt_km[i-1]) - h[i-1] * mu[i-1];
          mu[i] = h[i] / l[i];
          z[i] = (alpha[i] - h[i-1] * z[i-1]) / l[i];
        }
        
        l[n-1] = 1.0;
        z[n-1] = 0.0;
        atmSplineCoeffs[n-1] = 0.0;
        
        for (std::size_t j = n-1; j > 0; --j) {
          atmSplineCoeffs[j-1] = z[j-1] - mu[j-1] * atmSplineCoeffs[j];
        }
      }
    }
    CORSIKA_LOG_INFO("Loaded {} atmosphere data points from '{}'", atmAlt_km.size(), atmPath);
    if (!atmAlt_km.empty()) {
      CORSIKA_LOG_INFO("  Altitude range: {:.2f} - {:.2f} km", atmAlt_km.front(), atmAlt_km.back());
      CORSIKA_LOG_INFO("  n-1 range: {:.6e} - {:.6e}", atmNminus1.front(), atmNminus1.back());
      if (!atmSplineCoeffs.empty()) {
        CORSIKA_LOG_INFO("  Using cubic spline interpolation for smooth refractive index");
      }
    }
  } else {
    CORSIKA_LOG_WARN("Could not open atmosphere file '{}', using constant n.", atmPath);
  }

  // Create high-resolution lookup table by pre-interpolating the cubic spline
  // This provides precision better than 0.000001 for refractive index
  std::vector<double> atmAlt_km_fine;
  std::vector<double> atmNminus1_fine;
  std::vector<double> atmLnNminus1_fine;
  
  if (!atmAlt_km.empty() && !atmSplineCoeffs.empty()) {
    // Use fine step size: 0.01 km = 10 meters for high precision
    double const fine_step_km = 0.01;
    double const alt_min = atmAlt_km.front();
    double const alt_max = atmAlt_km.back();
    
    for (double alt_km = alt_min; alt_km <= alt_max; alt_km += fine_step_km) {
      // Find bracketing points using binary search
      auto it = std::upper_bound(atmAlt_km.begin(), atmAlt_km.end(), alt_km);
      if (it == atmAlt_km.begin() || it == atmAlt_km.end()) {
        // Edge cases - use boundary values
        if (it == atmAlt_km.begin()) {
          atmAlt_km_fine.push_back(alt_km);
          atmNminus1_fine.push_back(atmNminus1.front());
          atmLnNminus1_fine.push_back(atmLnNminus1.front());
        } else {
          atmAlt_km_fine.push_back(alt_km);
          atmNminus1_fine.push_back(atmNminus1.back());
          atmLnNminus1_fine.push_back(atmLnNminus1.back());
        }
        continue;
      }
      
      std::size_t i2 = std::distance(atmAlt_km.begin(), it);
      std::size_t i1 = i2 - 1;
      
      double const alt1 = atmAlt_km[i1];
      double const alt2 = atmAlt_km[i2];
      double const h = alt2 - alt1;
      
      // Cubic spline interpolation formula
      double const y1 = atmLnNminus1[i1];
      double const y2 = atmLnNminus1[i2];
      double const c1 = atmSplineCoeffs[i1];
      double const c2 = atmSplineCoeffs[i2];
      
      double const A = (alt2 - alt_km) / h;
      double const B = (alt_km - alt1) / h;
      
      double const ln_nminus1_interp = A * y1 + B * y2 + 
                                       ((A*A*A - A) * c1 + (B*B*B - B) * c2) * (h*h) / 6.0;
      double const nminus1_interp = std::exp(ln_nminus1_interp);
      
      atmAlt_km_fine.push_back(alt_km);
      atmNminus1_fine.push_back(nminus1_interp);
      atmLnNminus1_fine.push_back(ln_nminus1_interp);
    }
    
    CORSIKA_LOG_INFO("Created high-resolution lookup table with {} points (step={} km)", 
                     atmAlt_km_fine.size(), fine_step_km);
    CORSIKA_LOG_INFO("  Refractive index precision: better than 1e-6");
  }
  
  // Use the fine-resolution table if available, otherwise fall back to coarse table
  auto const& atmAlt_final = atmAlt_km_fine.empty() ? atmAlt_km : atmAlt_km_fine;
  auto const& atmNminus1_final = atmAlt_km_fine.empty() ? atmNminus1 : atmNminus1_fine;
  auto const& atmLnNminus1_final = atmAlt_km_fine.empty() ? atmLnNminus1 : atmLnNminus1_fine;

  // Lambda for refractive index: n(altitude, wavelength) -> pair<n_phase, n_group>
  // Returns both phase index (for Cherenkov angle) and group index (for timing)
  // When dispersion is enabled, applies wavelength-dependent correction factors:
  //   f_phase(λ) = 0.967 + 0.033 * (400/λ)^2.5
  //   f_group(λ) = 0.967 + 0.1155 * (400/λ)^2.5
  auto refrIndexFunc = [atmAlt = atmAlt_final,
                        atmNminus1Vals = atmNminus1_final,
                        atmLnNminus1Vals = atmLnNminus1_final,
                        useFineGrid = !atmAlt_km_fine.empty(),
                        seaLevelZ_m = (constants::EarthRadius::Mean / 1_m),
                        enableDispersion = !disable_dispersion](double abs_z_m, double wavelength_nm) 
                        -> std::pair<double, double> {
    // Compute base (n-1) at the given altitude
    double nminus1_base = 0.0003;  // fallback
    
    if (!atmAlt.empty()) {
      double alt_km = (abs_z_m - seaLevelZ_m) / 1000.0;
      
      // Edge case: below lowest altitude
      if (alt_km <= atmAlt.front()) {
        nminus1_base = atmNminus1Vals.front();
      }
      // Edge case: above highest altitude (vacuum with smooth exponential decay)
      else if (alt_km >= atmAlt.back()) {
        // Exponential extrapolation using last two points
        if (atmAlt.size() >= 2) {
          std::size_t n = atmAlt.size();
          double const alt_n = atmAlt[n-1];
          double const alt_n1 = atmAlt[n-2];
          double const ln_nminus1_n = atmLnNminus1Vals[n-1];
          double const ln_nminus1_n1 = atmLnNminus1Vals[n-2];
          
          // Scale height for exponential decay
          double const scale_height = (alt_n - alt_n1) / (ln_nminus1_n1 - ln_nminus1_n);
          double const extrapolated_ln_nminus1 = ln_nminus1_n - (alt_km - alt_n) / scale_height;
          
          // Safety check: don't go below 1e-10
          if (extrapolated_ln_nminus1 < std::log(1e-10)) {
            nminus1_base = 0.0;
          } else {
            nminus1_base = std::exp(extrapolated_ln_nminus1);
          }
        } else {
          nminus1_base = atmNminus1Vals.back();
        }
      }
      else {
        // Find bracketing points using binary search
        auto it = std::upper_bound(atmAlt.begin(), atmAlt.end(), alt_km);
        std::size_t i2 = std::distance(atmAlt.begin(), it);
        std::size_t i1 = i2 - 1;
        
        double const alt1 = atmAlt[i1];
        double const alt2 = atmAlt[i2];
        double const h = alt2 - alt1;
        double const t_local = alt_km - alt1;  // local coordinate within interval
        
        // For fine-resolution grid, simple linear interpolation is sufficient and fast
        // For coarse grid, would need cubic spline but fine grid gives us high precision
        double const t = t_local / h;
        double const ln_nminus1_interp = atmLnNminus1Vals[i1] + t * (atmLnNminus1Vals[i2] - atmLnNminus1Vals[i1]);
        nminus1_base = std::exp(ln_nminus1_interp);
      }
    }
    
    // Compute phase and group refractive indices
    double f_phase = 1.0;
    double f_group = 1.0;
    if (enableDispersion && wavelength_nm > 0.0) {
      double const ratio = 400.0 / wavelength_nm;
      double const power_term = std::pow(ratio, 2.5);
      f_phase = 0.967 + 0.033 * power_term;
      f_group = 0.967 + 0.1155 * power_term;
    }
    
    double const n_phase = 1.0 + nminus1_base * f_phase;
    double const n_group = 1.0 + nminus1_base * f_group;
    return {n_phase, n_group};
  };
  
  corsika::cherenkov::CherenkovProcess cherenkovProcess(
    /* n                 */ 1.0003,
    /* lambdaMin         */ 300.0,
    /* lambdaMax         */ 900.0,
    /* bunchSize         */ 5.0,
    /* groundZ_m         */ groundZ_m_abs,
    /* recorder          */ cherenkovRecorder,
    /* mode              */ projMode,
    /* showerAxis        */ std::make_shared<corsika::ShowerAxis>(showerAxis),
    /* showerCoreX_m     */ shower_core_x,
    /* showerCoreY_m     */ shower_core_y,
    /* showerCoreZ_m     */ shower_core_z,
    /* telescopeX_m      */ tel_x,
    /* telescopeY_m      */ tel_y,
    /* telescopeZ_m      */ tel_z,
    /* telescopeRadius_m */ tel_radius,
    /* refractiveIndexFn */ refrIndexFunc,
    /* earthRadiusM      */ constants::EarthRadius::Mean / 1_m,
    /* telescopeAzimuth  */ tel_azimuth_rad,
    /* telescopeZenith   */ tel_zenith_rad,
    /* showerAzimuth     */ phiRad,
    /* showerZenith      */ thetaRad,
    /* enableDispersion  */ !disable_dispersion,
    /* enableCurvature   */ enable_curvature,
    /* absorptionFile    */ app["--atmosphere-absorption-file"]->as<std::string>());

  
  std::stringstream args;
  for (int i = 0; i < argc; ++i) { args << argv[i] << " "; }
  // create the output manager that we then register outputs with
  OutputManager output(app["--filename"]->as<std::string>(), seed, args.str(),
                       compressOutput);

  // register energy losses as output
  EnergyLossWriter dEdX{showerAxis, dX};
  output.add("energyloss", dEdX);

  DynamicInteractionProcess<StackType> heModel;

  auto const all_elements = corsika::get_all_elements_in_universe(env);
  // have SIBYLL always for PROPOSAL photo-hadronic interactions
  auto sibyll = std::make_shared<corsika::sibyll::Interaction>(
      all_elements, corsika::setup::C7trackedParticles);

  if (auto const modelStr = app["--hadronModel"]->as<std::string>();
      modelStr == "SIBYLL-2.3d") {
    heModel = DynamicInteractionProcess<StackType>{sibyll};
  } else if (modelStr == "QGSJet-II.04") {
    heModel = DynamicInteractionProcess<StackType>{
        std::make_shared<corsika::qgsjetII::Interaction>()};
  } else if (modelStr == "QGSJet-III") {
    heModel = DynamicInteractionProcess<StackType>{
        std::make_shared<corsika::qgsjetIII::Interaction>()};
  } else if (modelStr == "EPOS-LHC") {
    heModel = DynamicInteractionProcess<StackType>{
        std::make_shared<corsika::epos::Interaction>(corsika::setup::C7trackedParticles)};
  } else if (modelStr == "EPOS-LHC-R") {
    heModel = DynamicInteractionProcess<StackType>{
        std::make_shared<corsika::EPOS_LHCR::Interaction>(
            corsika::setup::C7trackedParticles)};
  } else if (modelStr == "Pythia8") {
    heModel = DynamicInteractionProcess<StackType>{
        std::make_shared<corsika::pythia8::Interaction>(
            corsika::setup::C7trackedParticles)};
  } else {
    CORSIKA_LOG_CRITICAL("invalid choice \"{}\"; also check argument parser", modelStr);
    return EXIT_FAILURE;
  }

  InteractionCounter heCounted{heModel};

  corsika::pythia8::Decay decayPythia;
  // tau decay via TAUOLA (hard coded to left handed)
  corsika::tauola::Decay decayTauola(corsika::tauola::Helicity::LeftHanded);

  struct IsTauSwitch {
    bool operator()(const Particle& p) const {
      return (p.getPID() == Code::TauMinus || p.getPID() == Code::TauPlus);
    }
  };

  auto decaySequence = make_select(IsTauSwitch(), decayTauola, decayPythia);

  // neutrino interactions with pythia (options are: NC, CC)
  bool NC = false;
  bool CC = false;
  if (auto const nuIntStr = app["--neutrino-interaction-type"]->as<std::string>();
      nuIntStr == "neutral" || nuIntStr == "NC") {
    NC = true;
    CC = false;
  } else if (nuIntStr == "charged" || nuIntStr == "CC") {
    NC = false;
    CC = true;
  } else if (nuIntStr == "both") {
    NC = true;
    CC = true;
  }
  corsika::pythia8::NeutrinoInteraction neutrinoPrimaryPythia(
      corsika::setup::C7trackedParticles, NC, CC);

  // hadronic photon interactions in resonance region
  corsika::sophia::InteractionModel sophia;

  HEPEnergyType const emcut = 1_GeV * app["--emcut"]->as<double>();
  HEPEnergyType const hadcut = 1_GeV * app["--hadcut"]->as<double>();
  HEPEnergyType const mucut = 1_GeV * app["--mucut"]->as<double>();
  HEPEnergyType const taucut = 1_GeV * app["--taucut"]->as<double>();
  ParticleCut<SubWriter<decltype(dEdX)>> cut(emcut, emcut, hadcut, mucut, taucut,
                                             !track_neutrinos, dEdX);

  // tell proposal that we are interested in all energy losses above the particle cut
  auto const prod_threshold = std::min({emcut, hadcut, mucut, taucut});
  set_energy_production_threshold(Code::Electron, prod_threshold);
  set_energy_production_threshold(Code::Positron, prod_threshold);
  set_energy_production_threshold(Code::Photon, prod_threshold);
  set_energy_production_threshold(Code::MuMinus, prod_threshold);
  set_energy_production_threshold(Code::MuPlus, prod_threshold);
  set_energy_production_threshold(Code::TauMinus, prod_threshold);
  set_energy_production_threshold(Code::TauPlus, prod_threshold);

  // energy threshold for high energy hadronic model. Affects LE/HE switch for
  // hadron interactions and the hadronic photon model in proposal
  HEPEnergyType const heHadronModelThreshold =
      1_GeV * app["--hadronModelTransitionEnergy"]->as<double>();

  corsika::proposal::Interaction emCascade(
      env, sophia, sibyll->getHadronInteractionModel(), heHadronModelThreshold);

  // use BetheBlochPDG for hadronic continuous losses, and proposal otherwise
  corsika::proposal::ContinuousProcess<SubWriter<decltype(dEdX)>> emContinuousProposal(
      env, dEdX);
  BetheBlochPDG<SubWriter<decltype(dEdX)>> emContinuousBethe{dEdX};
  struct EMHadronSwitch {
    EMHadronSwitch() = default;
    bool operator()(const Particle& p) const { return is_hadron(p.getPID()); }
  };
  auto emContinuous =
      make_select(EMHadronSwitch(), emContinuousBethe, emContinuousProposal);

  LongitudinalWriter profile{showerAxis, dX};
  output.add("profile", profile);
  LongitudinalProfile<SubWriter<decltype(profile)>> longprof{profile};

  ProductionWriter prod_profile{showerAxis, dX};
  output.add("production_profile", prod_profile);
  ProductionProfile<SubWriter<decltype(prod_profile)>> prodprof{prod_profile};

// for ICRC2023
#ifdef WITH_FLUKA
  corsika::fluka::Interaction leIntModel{all_elements};
#else
  corsika::urqmd::UrQMD leIntModel{};
#endif
  InteractionCounter leIntCounted{leIntModel};

  // assemble all processes into an ordered process list
  struct EnergySwitch {
    HEPEnergyType cutE_;
    EnergySwitch(HEPEnergyType cutE)
        : cutE_(cutE) {}
    bool operator()(const Particle& p) const { return (p.getKineticEnergy() < cutE_); }
  };
  auto hadronSequence =
      make_select(EnergySwitch(heHadronModelThreshold), leIntCounted, heCounted);

  // observation plane; switch to vertical plane for near-horizontal showers
  double const zenithDeg = app["--zenith"]->as<double>();
  bool const isHorizontalShower = (zenithDeg > 85.0);

  DirectionVector obsPlaneNormal = DirectionVector(rootCS, {0., 0., 1.});
  DirectionVector obsPlaneXAxis = DirectionVector(rootCS, {1., 0., 0.});
  Point obsPlaneCenter = showerCore;

  if (isHorizontalShower) {
    // Use plane perpendicular to shower axis
    // Point the normal toward the incoming shower (opposite to propagation)
    obsPlaneNormal = DirectionVector(rootCS, {-nx, -ny, -nz});

    // Place the plane at the telescope location
    obsPlaneCenter = Point(rootCS, tel_x * 1_m, tel_y * 1_m, tel_z * 1_m);

    // Choose a horizontal x-axis perpendicular to the shower projection
    double const horizMag = std::sqrt(nx * nx + ny * ny);
    if (horizMag > 1e-10) {
      obsPlaneXAxis = DirectionVector(rootCS, {-ny / horizMag, nx / horizMag, 0.});
    }

    CORSIKA_LOG_INFO(
        "Using VERTICAL observation plane at telescope (zenith {:.2f} deg): center=({}, {}, {}), normal=({}, {}, {})",
        zenithDeg,
        obsPlaneCenter.getCoordinates().getX() / 1_m,
        obsPlaneCenter.getCoordinates().getY() / 1_m,
        obsPlaneCenter.getCoordinates().getZ() / 1_m,
        obsPlaneNormal.getComponents()[0],
        obsPlaneNormal.getComponents()[1],
        obsPlaneNormal.getComponents()[2]);
  } else {
    CORSIKA_LOG_INFO("Using HORIZONTAL observation plane (zenith {:.2f} deg)", zenithDeg);
  }

  Plane const obsPlane(obsPlaneCenter, obsPlaneNormal);
  ObservationPlane<TrackingType, ParticleWriterParquet> observationLevel{
      obsPlane, obsPlaneXAxis,
      true,       // plane should "absorb" particles
      1e-6 * 1_m, // ignored for absorbing planes
      false};     // do not print z-coordinate
  // register ground particle output
  output.add("particles", observationLevel);

  PrimaryWriter<TrackingType, ParticleWriterParquet> primaryWriter(observationLevel);
  output.add("primary", primaryWriter);

  int ring_number{app["--ring"]->as<int>()};
  auto const radius_{ring_number * 25_m};
  const int rr_ = static_cast<int>(radius_ / 1_m);

  // Radio observers and relevant information
  // the observer time variables
  const TimeType duration_{4e-7_s};
  const InverseTimeType sampleRate_{1e+9_Hz};

  // the observer collection for CoREAS and ZHS
  ObserverCollection<TimeDomainObserver> detectorCoREAS;
  ObserverCollection<TimeDomainObserver> detectorZHS;

  auto const showerCoreX_{showerCore.getCoordinates().getX()};
  auto const showerCoreY_{showerCore.getCoordinates().getY()};
  auto const injectionPosX_{injectionPos.getCoordinates().getX()};
  auto const injectionPosY_{injectionPos.getCoordinates().getY()};
  auto const injectionPosZ_{injectionPos.getCoordinates().getZ()};
  auto const triggerpoint_{Point(rootCS, injectionPosX_, injectionPosY_, injectionPosZ_)};

  if (ring_number != 0) {
    // setup CoREAS observers - use the for loop for star shape pattern
    for (auto phi_1 = 0; phi_1 <= 315; phi_1 += 45) {
      auto phiRad_1 = phi_1 / 180. * M_PI;
      auto const point_1{Point(rootCS, showerCoreX_ + radius_ * cos(phiRad_1),
                               showerCoreY_ + radius_ * sin(phiRad_1),
                               constants::EarthRadius::Mean)};
      std::cout << "Observer point CoREAS: " << point_1 << std::endl;
      auto triggertime_1{(triggerpoint_ - point_1).getNorm() / constants::c};
      std::string name_1 = "CoREAS_R=" + std::to_string(rr_) +
                           "_m--Phi=" + std::to_string(phi_1) + "degrees";
      TimeDomainObserver observer_1(name_1, point_1, rootCS, triggertime_1, duration_,
                                    sampleRate_, triggertime_1);
      detectorCoREAS.addObserver(observer_1);
    }

    // setup ZHS observers - use the for loop for star shape pattern
    for (auto phi_ = 0; phi_ <= 315; phi_ += 45) {
      auto phiRad_ = phi_ / 180. * M_PI;
      auto const point_{Point(rootCS, showerCoreX_ + radius_ * cos(phiRad_),
                              showerCoreY_ + radius_ * sin(phiRad_),
                              constants::EarthRadius::Mean)};
      std::cout << "Observer point ZHS: " << point_ << std::endl;
      auto triggertime_{(triggerpoint_ - point_).getNorm() / constants::c};
      std::string name_ =
          "ZHS_R=" + std::to_string(rr_) + "_m--Phi=" + std::to_string(phi_) + "degrees";
      TimeDomainObserver observer_2(name_, point_, rootCS, triggertime_, duration_,
                                    sampleRate_, triggertime_);
      detectorZHS.addObserver(observer_2);
    }
  }
  LengthType const step = 1_m;
  auto TP =
      make_tabulated_flat_atmosphere_radio_propagator(env, injectionPos, surface_, step);

  // initiate CoREAS
  RadioProcess<decltype(detectorCoREAS), CoREAS<decltype(detectorCoREAS), decltype(TP)>,
               decltype(TP)>
      coreas(detectorCoREAS, TP);

  // register CoREAS with the output manager
  output.add("CoREAS", coreas);

  // initiate ZHS
  RadioProcess<decltype(detectorZHS), ZHS<decltype(detectorZHS), decltype(TP)>,
               decltype(TP)>
      zhs(detectorZHS, TP);

  // register ZHS with the output manager
  output.add("ZHS", zhs);

  // make and register the first interaction writer
  InteractionWriter<setup::Tracking, ParticleWriterParquet> inter_writer(
      showerAxis, observationLevel);
  output.add("interactions", inter_writer);

  /* === END: SETUP PROCESS LIST === */

  // trigger the output manager to open the library for writing
  output.startOfLibrary();

  // loop over each shower
  for (int i_shower = 1; i_shower < nevent + 1; i_shower++) {

    CORSIKA_LOG_INFO("Shower {} / {} ", i_shower, nevent);

    // randomize the primary energy
    double const eSlope = app["--eslope"]->as<double>();
    PowerLawDistribution<HEPEnergyType> powerLawRng(eSlope, eMin, eMax);
    HEPEnergyType const primaryTotalEnergy =
        (eMax == eMin) ? eMin
                       : powerLawRng(RNGManager<>::getInstance().getRandomStream(
                             "primary_particle"));

    auto const eKin = primaryTotalEnergy - get_mass(beamCode);

    // set up thinning based on primary parameters
    double const emthinfrac = app["--emthin"]->as<double>();
    double const maxWeight = std::invoke([&]() {
      if (auto const wm = app["--max-weight"]->as<double>(); wm > 0)
        return wm;
      else
        return 0.5 * emthinfrac * primaryTotalEnergy / 1_GeV;
    });
    EMThinning thinning{emthinfrac * primaryTotalEnergy, maxWeight, !multithin};

    // set up the stack inspector
    StackInspector<StackType> stackInspect(10000, false, primaryTotalEnergy);

    // assemble the final process sequence
    auto sequence =
        make_sequence(stackInspect, neutrinoPrimaryPythia, hadronSequence, decaySequence,
                      emCascade,  prodprof, emContinuous, cherenkovProcess, coreas, zhs, longprof,
                      observationLevel, inter_writer, thinning, cut);

    // create the cascade object using the default stack and tracking
    // implementation
    TrackingType tracking(app["--max-deflection-angle"]->as<double>());
    StackType stack;
    Cascade EAS(env, tracking, sequence, output, stack);

    // setup particle stack, and add primary particle
    stack.clear();

    // print our primary parameters all in one place
    CORSIKA_LOG_INFO("Primary name:         {}", beamCode);
    if (app["--pdg"]->count() > 0) {
      CORSIKA_LOG_INFO("Primary PDG ID:       {}", app["--pdg"]->as<int>());
    } else {
      CORSIKA_LOG_INFO("Primary Z/A:          {}/{}", Z, A);
    }
    CORSIKA_LOG_INFO("Primary Total Energy: {}", primaryTotalEnergy);
    CORSIKA_LOG_INFO("Primary Momentum:     {}",
                     calculate_momentum(primaryTotalEnergy, get_mass(beamCode)));
    CORSIKA_LOG_INFO("Primary Direction:    {}", propDir.getNorm());
    CORSIKA_LOG_INFO("Point of Injection:   {}", injectionPos.getCoordinates());
    CORSIKA_LOG_INFO("Shower Axis Length:   {}",
                     (showerCore - injectionPos).getNorm() * 1.2);

    // add the desired particle to the stack
    auto const primaryProperties =
        std::make_tuple(beamCode, eKin, propDir.normalized(), injectionPos, 0_ns);
    stack.addParticle(primaryProperties);

    // if we want to fix the first location of the shower
    if (force_interaction) {
      CORSIKA_LOG_INFO("Fixing first interaction at injection point.");
      EAS.forceInteraction();
    }

    if (force_decay) {
      CORSIKA_LOG_INFO("Forcing the primary to decay");
      EAS.forceDecay();
    }

    primaryWriter.recordPrimary(primaryProperties);

    // run the shower
    EAS.run();

    HEPEnergyType const Efinal =
        dEdX.getEnergyLost() + observationLevel.getEnergyGround();

    CORSIKA_LOG_INFO(
        "total energy budget (GeV): {} (dEdX={} ground={}), "
        "relative difference (%): {}",
        Efinal / 1_GeV, dEdX.getEnergyLost() / 1_GeV,
        observationLevel.getEnergyGround() / 1_GeV,
        (Efinal / primaryTotalEnergy - 1) * 100);

    if (!disable_interaction_hists) {
      CORSIKA_LOG_INFO("Saving interaction histograms");
      auto const hists = heCounted.getHistogram() + leIntCounted.getHistogram();

      // directory for output of interaction histograms
      string const outdir(app["--filename"]->as<std::string>() + "/interaction_hist");
      boost::filesystem::create_directories(outdir);

      string const labHist_file = outdir + "/inthist_lab_" + to_string(i_shower) + ".npz";
      string const cMSHist_file = outdir + "/inthist_cms_" + to_string(i_shower) + ".npz";
      save_hist(hists.labHist(), labHist_file, true);
      save_hist(hists.CMSHist(), cMSHist_file, true);
    }
    CORSIKA_LOG_INFO("Total Cherenkov hits recorded: {}", cherenkovHitCount);
    if (cherenkovCsvPtr) {
      cherenkovCsvPtr->flush();
    }

  }

  // and finalize the output on disk
  output.endOfLibrary();

  // Write eventio file if that format was selected
  if (useEventio && eventioWriterPtr) {
    CORSIKA_LOG_INFO("Writing Cherenkov eventio file...");
    try {
      eventioWriterPtr->write();
      CORSIKA_LOG_INFO("Cherenkov eventio file written successfully");
    } catch (const std::exception& e) {
      CORSIKA_LOG_WARN("Failed to write eventio file: {}", e.what());
    }
  } else if (cherenkovCsvPtr) {
    cherenkovCsvPtr->close();
    CORSIKA_LOG_INFO("Cherenkov CSV file closed");
  }

  auto const run_end = std::chrono::steady_clock::now();
  auto const run_seconds =
      std::chrono::duration_cast<std::chrono::duration<double>>(run_end - run_start).count();
  CORSIKA_LOG_INFO("Run time: {:.3f} s", run_seconds);

  return EXIT_SUCCESS;
}
