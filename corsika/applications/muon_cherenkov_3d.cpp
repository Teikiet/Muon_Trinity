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
#include <corsika/framework/geometry/CoordinateSystemRegistry.hpp>
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
#include <corsika/media/magnetic/GeomagneticModel.hpp>
#include <corsika/media/refractivity/GladstoneDaleRefractiveIndex.hpp>
#include <corsika/media/density_and_composition/HomogeneousMedium.hpp>
#include <corsika/media/interfaces/IMagneticFieldModel.hpp>
#include <corsika/media/LayeredSphericalAtmosphereBuilder.hpp>
#include <corsika/media/medium/MediumPropertyModel.hpp>
#include <corsika/media/composition/NuclearComposition.hpp>
#include <corsika/media/ShowerAxis.hpp>
#include <corsika/media/magnetic/UniformMagneticField.hpp>

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
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

#include <chrono>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <limits>
#include <string>
#include <optional>
#include <set>

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

using EnvironmentInterface = media::IRefractiveIndexModel<
    media::IMediumPropertyModel<media::IMagneticFieldModel<media::IMediumModel>>>;
using EnvType = media::Environment<EnvironmentInterface>;
using StackType = setup::Stack<EnvType>;
using TrackingType = setup::Tracking;
using Particle = StackType::particle_type;

// ============================================================================
// 3-D observation-plane particle recorder (self-contained)
// ----------------------------------------------------------------------------
// A ContinuousProcess + parquet writer that records particles crossing a stack
// of evenly-spaced (in altitude) observation planes, for the 3-D time evolution
// of the shower. It never absorbs particles.
// ============================================================================
class ThreeDParticleWriter : public BaseOutput {
public:
  ParquetStreamer output_;
  unsigned int showerId_ = 0;

  // Origin of the global (X, Y, Z) output columns. Defaults to the root
  // coordinate system origin (Earth centre); set to the telescope centre so the
  // recorded global coordinates are telescope-centred.
  LengthType originX_ = 0_m;
  LengthType originY_ = 0_m;
  LengthType originZ_ = 0_m;

  void setOutputOrigin(Point const& origin) {
    auto const c = origin.getCoordinates();
    originX_ = c.getX();
    originY_ = c.getY();
    originZ_ = c.getZ();
  }

  void startOfLibrary(boost::filesystem::path const& directory) override {
    output_.initStreamer((directory / "particles_3d.parquet").string());
    output_.enableCompression();
    output_.addField("slice", parquet::Repetition::REQUIRED, parquet::Type::INT32,
                     parquet::ConvertedType::INT_32);
    output_.addField("pdg", parquet::Repetition::REQUIRED, parquet::Type::INT32,
                     parquet::ConvertedType::INT_32);
    output_.addField("kinetic_energy", parquet::Repetition::REQUIRED,
                     parquet::Type::FLOAT, parquet::ConvertedType::NONE);
    output_.addField("x", parquet::Repetition::REQUIRED, parquet::Type::FLOAT,
                     parquet::ConvertedType::NONE);
    output_.addField("y", parquet::Repetition::REQUIRED, parquet::Type::FLOAT,
                     parquet::ConvertedType::NONE);
    output_.addField("X", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,
                     parquet::ConvertedType::NONE);
    output_.addField("Y", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,
                     parquet::ConvertedType::NONE);
    output_.addField("Z", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,
                     parquet::ConvertedType::NONE);
    output_.addField("nx", parquet::Repetition::REQUIRED, parquet::Type::FLOAT,
                     parquet::ConvertedType::NONE);
    output_.addField("ny", parquet::Repetition::REQUIRED, parquet::Type::FLOAT,
                     parquet::ConvertedType::NONE);
    output_.addField("nz", parquet::Repetition::REQUIRED, parquet::Type::FLOAT,
                     parquet::ConvertedType::NONE);
    output_.addField("time", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,
                     parquet::ConvertedType::NONE);
    output_.addField("weight", parquet::Repetition::REQUIRED, parquet::Type::FLOAT,
                     parquet::ConvertedType::NONE);
    output_.addField("h", parquet::Repetition::REQUIRED, parquet::Type::FLOAT,
                     parquet::ConvertedType::NONE);
    output_.buildStreamer();
    showerId_ = 0;
  }

  void startOfShower(unsigned int const showerId) override { showerId_ = showerId; }
  void endOfShower(unsigned int const) override {}
  void endOfLibrary() override { output_.closeStreamer(); }
  YAML::Node getConfig() const override { return YAML::Node(); }

  /// @param x, y      plane-local coordinates of the crossing
  /// @param X, Y, Z   absolute (root coordinate system) position of the crossing;
  ///                  written relative to the configured output origin
  void write(Code const pid, HEPEnergyType const kineticEnergy, LengthType const x,
             LengthType const y, LengthType const X, LengthType const Y,
             LengthType const Z, double const nx, double const ny, double const nz,
             TimeType const t, double const weight, LengthType const h,
             int const slice) {
    auto const w = output_.getWriter();
    *w << showerId_ << slice << static_cast<int>(get_PDG(pid))
       << static_cast<float>(kineticEnergy / 1_GeV) << static_cast<float>(x / 1_m)
       << static_cast<float>(y / 1_m) << static_cast<double>((X - originX_) / 1_m)
       << static_cast<double>((Y - originY_) / 1_m)
       << static_cast<double>((Z - originZ_) / 1_m) << static_cast<float>(nx)
       << static_cast<float>(ny) << static_cast<float>(nz)
       << static_cast<double>(t / 1_s) << static_cast<float>(weight)
       << static_cast<float>(h / 1_m) << parquet::EndRow;
  }
};

class ThreeDPlaneRecorder : public ContinuousProcess<ThreeDPlaneRecorder>,
                            public ThreeDParticleWriter {
public:
  ThreeDPlaneRecorder(std::vector<Point> const& centers,
                      std::vector<DirectionVector> const& normals,
                      std::vector<DirectionVector> const& xaxes,
                      std::vector<DirectionVector> const& yaxes,
                      std::vector<double> const& hSlice,
                      Point const& earthCenter, LengthType const rEarth,
                      double const eminGeV, std::set<Code> const& filter,
                      Point const& outputOrigin)
      : ThreeDParticleWriter()
      , centers_(centers)
      , normals_(normals)
      , xAxes_(xaxes)
      , yAxes_(yaxes)
      , hSlice_(hSlice)
      , earthCenter_(earthCenter)
      , rEarth_(rEarth)
      , eminGeV_(eminGeV)
      , filter_(filter) {
    setOutputOrigin(outputOrigin);

    // Build the analytic description of the slice stack (see the member comments).
    if (centers_.size() >= 2) {
      sliceNormal_ = normals_[0].normalized();
      DirectionVector const n = *sliceNormal_;
      dStep_m_ = n.dot(centers_[1] - centers_[0]) / 1_m;
      analyticOk_ = (dStep_m_ != 0.0);
      for (std::size_t k = 1; analyticOk_ && k < centers_.size(); ++k) {
        double const dk = n.dot(centers_[k] - centers_[0]) / 1_m;
        double const nk = n.dot(normals_[k]); // dimensionless
        if (std::abs(dk - double(k) * dStep_m_) > 1e-6 ||
            std::abs(nk - 1.0) > 1e-9) {
          analyticOk_ = false;
        }
      }
    }
  }

  void startOfLibrary(boost::filesystem::path const& directory) override {
    ThreeDParticleWriter::startOfLibrary(directory);
    boost::filesystem::path yp = directory / "particles_3d_slices.yaml";
    std::ofstream f(yp.string());
    if (f) {
      f << "earth_radius_m: " << std::setprecision(9) << (rEarth_ / 1_m) << "\n";
      f << "output_origin_m: [" << std::setprecision(9) << (originX_ / 1_m) << ", "
        << (originY_ / 1_m) << ", " << (originZ_ / 1_m) << "]\n";
      f << "slices:\n";
      for (int k = 0; k < int(centers_.size()); ++k) {
        auto const c = centers_[k].getCoordinates();
        auto const nq = normals_[k].getComponents();
        auto const xq = xAxes_[k].getComponents();
        auto const yq = yAxes_[k].getComponents();
        // plain doubles, so the YAML holds bare numbers (a phys::units quantity
        // would print with an (empty) unit suffix, e.g. "-0.866 []")
        double const cx = c.getX() / 1_m, cy = c.getY() / 1_m, cz = c.getZ() / 1_m;
        double const nxv = nq[0], nyv = nq[1], nzv = nq[2];
        double const xxv = xq[0], xyv = xq[1], xzv = xq[2];
        double const yxv = yq[0], yyv = yq[1], yzv = yq[2];
        f << "  - {slice: " << k << ", h_m: " << std::setprecision(9) << hSlice_[k]
          << ", center_m: [" << cx << ", " << cy << ", " << cz << "]"
          << ", normal: [" << nxv << ", " << nyv << ", " << nzv << "]"
          << ", x_axis: [" << xxv << ", " << xyv << ", " << xzv << "]"
          << ", y_axis: [" << yxv << ", " << yyv << ", " << yzv << "]}\n";
      }
    }
    showerId_ = 0;
  }

  template <typename TParticle, typename TTrajectory>
  LengthType getMaxStepLength(TParticle const& particle,
                              TTrajectory const& trajectory) {
    LengthType const inf = std::numeric_limits<double>::infinity() * 1_m;
    double const duration_s = trajectory.getDuration() / 1_s;
    if (!(duration_s > 0.0)) { return inf; }

    if (analyticOk_) {
      // The slices form a uniform stack of parallel planes at signed distances
      // k*dStep from slice 0 (k = 0..N-1). With q = distance/dStep being the
      // fractional slice index, the next crossing is the next integer in the
      // direction of travel -- so we only need to test a small window of candidate
      // planes around that index instead of all N of them. O(1) instead of O(N).
      //
      // The window (rather than a single plane) makes this robust against the
      // straight-line index estimate being off by one; and since doContinuous
      // records on the sign change of the signed distance, an imperfect estimate
      // here can never lose or duplicate a record.
      DirectionVector const n = *sliceNormal_;
      double const d0_m = n.dot(particle.getPosition() - centers_[0]) / 1_m;
      double const q = d0_m / dStep_m_; // fractional slice index of the particle

      auto const velocity =
          particle.getMomentum() / particle.getEnergy() * constants::c;
      double const vn_m_s = n.dot(velocity) / (1_m / 1_s); // m/s along the normal
      if (vn_m_s == 0.0) { return inf; } // parallel to the planes -> never crosses

      double const qdot = vn_m_s / dStep_m_; // dq/dt, sign gives the direction
      double const kCenter = (qdot > 0.0) ? std::floor(q) + 1.0 : std::ceil(q) - 1.0;
      double const nSlice = double(centers_.size());
      if (kCenter < -3.0 || kCenter > nSlice + 3.0) { return inf; } // far past stack

      LengthType best = inf;
      for (int dk = -2; dk <= 2; ++dk) {
        int const kk = int(kCenter) + dk;
        if (kk < 0 || kk >= int(centers_.size())) { continue; }
        Plane const pl(centers_[kk], normals_[kk]);
        auto const inter = TrackingType::intersect(particle, pl);
        double const tEntry_s = inter.getEntry() / 1_s;
        if (!(tEntry_s > epsS_)) { continue; }
        if (!(tEntry_s <= duration_s)) { continue; }
        LengthType const L = trajectory.getLength(tEntry_s / duration_s);
        if (L < best) { best = L; }
      }
      return best;
    }

    // --- fallback: per-plane loop (only for non-uniform or single-slice stacks) ---
    LengthType best = inf;
    for (int k = 0; k < int(centers_.size()); ++k) {
      Plane const pl(centers_[k], normals_[k]);
      auto const inter = TrackingType::intersect(particle, pl);
      double const tEntry_s = inter.getEntry() / 1_s;
      if (!(tEntry_s > epsS_)) { continue; }
      if (!(tEntry_s <= duration_s)) { continue; }
      double const frac = tEntry_s / duration_s;
      LengthType const L = trajectory.getLength(frac);
      if (L < best) { best = L; }
    }
    return best;
  }

  template <typename TParticle>
  ProcessReturn doContinuous(Step<TParticle>& step, bool /*limitFlag*/) {
    Point const pPre = step.getPositionPre();
    Point const pPost = step.getPositionPost();
    TimeType const tPre = step.getTimePre();
    TimeType const tPost = step.getTimePost();
    auto const dirc = step.getDirectionPost().getComponents();
    double const nx = dirc[0], ny = dirc[1], nz = dirc[2];
    Code const pid = step.getParticlePre().getPID();
    double const ekinGeV = step.getEkinPost() / 1_GeV;
    double const weight = step.getParticlePre().getWeight();

    if (ekinGeV < eminGeV_) { return ProcessReturn::Ok; }
    if (!filter_.empty() && filter_.count(pid) == 0) { return ProcessReturn::Ok; }

    // Locate the slice range that can be crossed during this step. For a uniform
    // parallel stack the signed distances satisfy d_k(p) = d_0(p) - k*dStep, so the
    // slices whose distance changes sign form a contiguous integer range that is
    // found in O(1) -- no need to scan all N planes.
    int kLo = 0;
    int kHi = int(centers_.size()) - 1;
    if (analyticOk_) {
      DirectionVector const n = *sliceNormal_;
      double const dPre0 = n.dot(pPre - centers_[0]) / 1_m;
      double const dPost0 = n.dot(pPost - centers_[0]) / 1_m;
      double const lo = std::min(dPre0, dPost0);
      double const hi = std::max(dPre0, dPost0);
      double const s = dStep_m_;
      // k*dStep in [lo, hi]
      double const ka = (s > 0.0) ? lo / s : hi / s;
      double const kb = (s > 0.0) ? hi / s : lo / s;
      double kLoD = std::ceil(ka);
      double kHiD = std::floor(kb);
      if (!(kHiD >= 0.0) || !(kLoD <= double(centers_.size()) - 1.0)) {
        return ProcessReturn::Ok; // this step crosses no slice
      }
      kLoD = std::max(kLoD, 0.0);
      kHiD = std::min(kHiD, double(centers_.size()) - 1.0);
      kLo = int(kLoD);
      kHi = int(kHiD);
    }

    for (int k = kLo; k <= kHi; ++k) {
      DirectionVector const nk = normals_[k];
      LengthType const dPre = nk.dot(pPre - centers_[k]);
      LengthType const dPost = nk.dot(pPost - centers_[k]);
      bool const crossed =
          (dPre > 0_m && dPost <= 0_m) || (dPre < 0_m && dPost >= 0_m);
      if (!crossed) { continue; }
      double const frac = dPre / (dPre - dPost);
      Point const pos = pPre + (pPost - pPre) * frac;
      TimeType const t = tPre + (tPost - tPre) * frac;
      LengthType const xl = xAxes_[k].dot(pos - centers_[k]);
      LengthType const yl = yAxes_[k].dot(pos - centers_[k]);
      LengthType const h = (pos - earthCenter_).getNorm() - rEarth_;
      this->write(pid, step.getEkinPost(), xl, yl,
                  pos.getCoordinates().getX(), pos.getCoordinates().getY(),
                  pos.getCoordinates().getZ(), nx, ny, nz, t, weight, h, k);
    }
    return ProcessReturn::Ok;
  }

private:
  std::vector<Point> centers_;
  std::vector<DirectionVector> normals_;
  std::vector<DirectionVector> xAxes_;
  std::vector<DirectionVector> yAxes_;
  std::vector<double> hSlice_;
  Point earthCenter_;
  LengthType rEarth_;
  double eminGeV_ = 0.0;
  std::set<Code> filter_;
  double const epsS_ = 1e-12; // s ; ignore crossings closer than this

  // --- analytic slice description -------------------------------------------
  // The slices are a uniformly spaced stack of parallel planes, so the whole
  // stack is described by one normal and one signed spacing. This lets
  // getMaxStepLength() find the next crossing in O(1) instead of testing every
  // plane. analyticOk_ is false if the stack is not uniform (defensive: the
  // current slice construction always produces a uniform stack), in which case
  // the original per-plane loop is used.
  std::optional<DirectionVector> sliceNormal_;
  double dStep_m_ = 0.0;      // signed spacing along sliceNormal_
  bool analyticOk_ = false;
};

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
using MyExtraEnv = media::GladstoneDaleRefractiveIndex<
    media::MediumPropertyModel<media::UniformMagneticField<T>>>;

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
      ->check(CLI::Range(0., 100.))
      ->group("Primary");
  app.add_option("-a,--azimuth", "Primary azimuth angle (deg)")
      ->default_val(0.)
      ->check(CLI::Range(0., 360.))
      ->group("Primary");
  // Injection geometry overrides (primary only; nominal shower axis is unchanged)
  // NOTE: no CLI range limit here on purpose. The offset is validated physically
  // further down, by checking the radial distance of the resulting starting point
  // against the extent of the Earth atmosphere.
  double injection_x_m = 0.0;
  app.add_option("--injection-x", injection_x_m,
                 "Horizontal offset (m) of the primary injection point in global x. "
                 "0 = on the nominal shower axis.")
      ->default_val(0.0)
      ->group("Primary");
  double injection_y_m = 0.0;
  app.add_option("--injection-y", injection_y_m,
                 "Horizontal offset (m) of the primary injection point in global y. "
                 "0 = on the nominal shower axis.")
      ->default_val(0.0)
      ->group("Primary");
  double injection_zenith_deg = -1.0;
  app.add_option("--injection-zenith", injection_zenith_deg,
                 "Zenith angle (deg) of the primary's injection direction. "
                 "-1 = use the nominal shower zenith (--zenith).")
      ->default_val(-1.0)
      ->check(CLI::Range(-1., 100.))
      ->group("Primary");
  double injection_azimuth_deg = -1.0;
  app.add_option("--injection-azimuth", injection_azimuth_deg,
                 "Azimuth angle (deg) of the primary's injection direction. "
                 "-1 = use the nominal shower azimuth (--azimuth).")
      ->default_val(-1.0)
      ->check(CLI::Range(-1., 360.))
      ->group("Primary");
  // Which intersection of the shower axis with the injection sphere to use. The ray
  // through the observation point meets the injection sphere twice; for upgoing showers
  // (zenith > 90 deg) with the injection height BELOW the observation level both roots are
  // upstream and only the NEAR one is physical (Earth-skimming geometry). "auto" selects
  // it whenever it is upstream, otherwise the far root is used (standard downgoing case).
  std::string injectionRoot("auto");
  app.add_option("--injection-root", injectionRoot,
                 "Intersection of the shower axis with the injection sphere to use: "
                 "auto (near whenever it is upstream, i.e. upgoing Earth-skimming; else "
                 "far), near, or far.")
      ->default_val("auto")
      ->check(CLI::IsMember({"auto", "near", "far"}))
      ->group("Primary");

  //////// Config options ////////

  app.add_option("--emcut",
                 "Min. kin. energy of photons, electrons and "
                 "positrons in tracking (GeV)")
      ->default_val(0.015)
      ->check(CLI::Range(0.000001, 1.e13))
      ->group("Config");
  app.add_option("--hadcut", "Min. kin. energy of hadrons in tracking (GeV)")
      ->default_val(0.3)
      ->check(CLI::Range(0.02, 1.e13))
      ->group("Config");
  app.add_option("--mucut", "Min. kin. energy of muons in tracking (GeV)")
      ->default_val(0.015)
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
  double telescope_position_correction = 0.1;
  app.add_option("--telescope-position-correction", telescope_position_correction,
                 "Fraction of Cherenkov steps to re-center the telescope using the most energetic charged particle (0 disables).")
      ->check(CLI::Range(0.0, 1.0))
      ->group("Cherenkov");
  app.add_option("--atmosphere-absorption-file",
                 "Path to atmospheric absorption table file for Cherenkov photons. "
                 "When specified, applies Monte Carlo absorption based on the table. "
                 "File format: text file with header 'n_wavelengths n_heights', "
                 "followed by wavelength array [nm], height array [m], and ecoeff matrix. "
                 "Use convert_npz_to_dat.py to convert CHASM abstable.npz files.")
      ->default_val("")
      ->group("Cherenkov");
  bool disable_fast_timing = false;
  bool disable_curved_table = false;
  app.add_flag("--disable-cherenkov-curved-table", disable_curved_table,
               "Use the legacy curved-ray tracer instead of the CurvedRayTable. "
               "Only affects eligible curved-ground tracing; the table is the "
               "default for those rays.");
  // ---- Phase 3: exact-tracer integrator ---------------------------------------
  std::string tracer_integrator = "auto";
  app.add_option("--cherenkov-tracer", tracer_integrator,
                 "Integrator for the EXACT curved-ray tracer: dopri5 (adaptive 5(4), "
                 "default), rk4 (fixed-step reference), tangent (Phase 4, diagnostic: "
                 "regularised turning-point coordinates only, non-turning rays are "
                 "DROPPED), or auto (Phase 4 production: try tangent, fall back to "
                 "dopri5 per ray when it refuses -- a counted fallback, so no ray is "
                 "ever lost).")
      ->default_val("auto")
      ->check(CLI::IsMember({"rk4", "dopri5", "tangent", "auto"}))
      ->group("Cherenkov");
  double tracer_max_step_m = 20000.0;
  app.add_option("--cherenkov-tracer-max-step", tracer_max_step_m,
                 "Maximum adaptive step [m] for the dopri5 tracer. 0 = the RK4 "
                 "heuristic caps (H_local/20, 1000 m, plus the radial cap), which "
                 "reproduces the RK4 step pattern but not its cost.  Measured on the "
                 "frozen benchmark: 5 km -> 4.0x, 20 km -> 9.1x, 50 km -> 11.1x faster "
                 "than RK4, with max impact error 2.5 um / 0.12 mm / 0.12 mm.")
      ->default_val(20000.0)
      ->check(CLI::NonNegativeNumber)
      ->group("Cherenkov");
  double tracer_atol_pos_m = 1.0e-3;
  app.add_option("--cherenkov-tracer-atol-pos", tracer_atol_pos_m,
                 "Absolute position tolerance per component [m] for the dopri5 error "
                 "control (validated at 1 mm; the benchmark's max impact error against "
                 "a converged reference is 0.12 mm at this setting).")
      ->default_val(1.0e-3)
      ->check(CLI::PositiveNumber)
      ->group("Cherenkov");
  double tracer_atol_dir = 1.0e-9;
  app.add_option("--cherenkov-tracer-atol-dir", tracer_atol_dir,
                 "Absolute direction-cosine tolerance for the dopri5 error control.")
      ->default_val(1.0e-9)
      ->check(CLI::PositiveNumber)
      ->group("Cherenkov");
  double tracer_atol_opt_m = 1.0e-4;
  app.add_option("--cherenkov-tracer-atol-optical", tracer_atol_opt_m,
                 "Absolute group-path (optical path) tolerance [m] for the dopri5 "
                 "error control.")
      ->default_val(1.0e-4)
      ->check(CLI::PositiveNumber)
      ->group("Cherenkov");
  int trace_ab_count = 0;
  std::string cher_table_cache;
  app.add_option("--cherenkov-table-cache", cher_table_cache,
                 "Directory for the on-disk curved-table cache. The geometry table "
                 "(31 s, 229 MB) and the absorption table (18 s, 19 MB) are deterministic "
                 "functions of the atmosphere profile and the configuration, so they can be "
                 "built once and read back in ~0.3 s. A cache file is used only when its "
                 "embedded key matches this run's parameters and sampled profile; otherwise "
                 "the table is rebuilt and the file overwritten. Empty = no cache.")
      ->default_val("")
      ->group("Cherenkov");
  double cher_grazing_clearance = 2400.0;
  app.add_option("--cherenkov-grazing-clearance", cher_grazing_clearance,
                 "Grazing rejection band [m]: rays passing within this clearance of their "
                 "apparent horizon are refused by the geometry table and traced exactly. The "
                 "impact error is a cliff there: 1.05 m at mu 0.029-0.030, 0.66 m at 0.030-"
                 "0.031, 0.24 m at 0.031-0.033, then 0.011 m above 0.033 and 0.006 m above "
                 "0.040. The band is calibrated at 164 m of clearance per 0.001 of mu, so 2400 "
                 "m puts the cut at mu ~= 0.0347 and keeps every table answer millimetre-"
                 "accurate. 0 = structural guard only (pre-wave-5 behaviour).")
      ->default_val(2400.0)
      ->group("Cherenkov");
  app.add_option("--cherenkov-trace-ab", trace_ab_count,
                 "Diagnostic (2g-h2n step A): for the first N Telescope-mode photons that "
                 "the CURVED TABLE answers, also run the plane integrator from the same "
                 "emission state and write both to <cherdir>/cherenkov_ab.csv. "
                 "0 = off; the production result is unchanged.")
      ->default_val(0)
      ->group("Cherenkov");
  app.add_flag("--disable-cherenkov-fast-timing", disable_fast_timing,
               "Disable the closed-form straight-line timing path (validation reference)")
      ->group("Cherenkov");
  bool enable_cherenkov_profile = false;
  app.add_flag("--cherenkov-profile", enable_cherenkov_profile,
               "Enable the per-step Cherenkov profiler (diagnostics)")
      ->group("Cherenkov");
  int cherenkov_timing_segments = 0;
  app.add_option("--cherenkov-timing-segments", cherenkov_timing_segments,
                 "Fixed integration-segment count in the non-fast timing path "
                 "(0 = heuristic; e.g. 10000 as a high-resolution validation reference)")
      ->default_val(0)
      ->check(CLI::NonNegativeNumber)
      ->group("Cherenkov");
  double cherenkov_bunch_size = 5.0;
  app.add_option("--cherenkov-bunch-size", cherenkov_bunch_size,
                 "Reference Cherenkov bunch size (photons/bunch); with "
                 "--cherenkov-bunch-size-auto this is the value at 1 TeV")
      ->default_val(5.0)
      ->check(CLI::PositiveNumber)
      ->group("Cherenkov");
  bool no_bunch_size_auto = false;
  app.add_flag("--no-cherenkov-bunch-size-auto", no_bunch_size_auto,
               "Disable energy-scaling of the Cherenkov bunch size")
      ->group("Cherenkov");
  double cherenkov_bunch_size_max = 100.0;
  app.add_option("--cherenkov-bunch-size-max", cherenkov_bunch_size_max,
                 "Cap for the auto Cherenkov bunch size (Telescope mode is "
                 "additionally capped at 100 so one bunch cannot fill the camera)")
      ->default_val(100.0)
      ->check(CLI::PositiveNumber)
      ->group("Cherenkov");
  // 3-D observation-plane recorder options
  int observationSteps = 0;
  app.add_option("--observation-steps", observationSteps,
                 "Number of evenly-spaced 3-D observation slices between the detector "
                 "and the injection height. 0 = legacy single-plane behavior.")
      ->default_val(0)
      ->check(CLI::NonNegativeNumber)
      ->group("Config");
  double obs3d_emin = 0.0;
  app.add_option("--obs3d-emin", obs3d_emin,
                 "Min kinetic energy (GeV) of particles recorded by the 3-D recorder "
                 "(0 = no filter).")
      ->default_val(0.0)
      ->check(CLI::NonNegativeNumber)
      ->group("Config");
  std::string obs3d_particles("");
  app.add_option("--obs3d-particles", obs3d_particles,
                 "Optional comma-separated PDG codes or the preset 'muons'. Empty = all.")
      ->default_val("")
      ->group("Config");
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
  const CoordinateSystem& rootCS = env.getCoordinateSystem();
  Point const center{rootCS, 0_m, 0_m, 0_m};
  Point const surface_{rootCS, 0_m, 0_m, constants::EarthRadius::Mean};
  media::GeomagneticModel wmm(rootCS, center, corsika_data("GeoMag/WMM.COF"));

  // build an atmosphere with Keilhauer's parametrization of the
  // US standard atmosphere into `env`
    media::create_5layer_atmosphere<EnvironmentInterface, MyExtraEnv>(
      env, media::AtmosphereId::USStdBK, center, 1.000327, surface_,
      media::Medium::AirDry1Atm,
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

  // Only neutrino primaries need the neutrino interaction process in the sequence.
  bool const neutrinoPrimary = is_neutrino(beamCode);

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

  Point const showerCore{rootCS, 0_m, 0_m, observationHeight};

  // Compute distances from shower core to the injection sphere (radius = injectionHeight)
  // Sphere intersection equation: |core + s·dir|² = injH²
  //   →  s² + 2·(core·dir)·s + obsH² - injH² = 0
  // For backwardDir (opposite to propDir): core·backwardDir = obsH·cosθ
  //   →  s² + 2·obsH·cosθ·s + obsH² - injH² = 0
  // Roots: s₁₂ = -obsH·cosθ ± √(injH² - obsH²·sin²θ)
  //
  // For θ ≤ 90° (downgoing): both cosθ ≥ 0 and sinθ ≥ 0
  //   → the +root gives the far intersection (valid injection above core)
  //
  // For θ > 90° (upgoing): cosθ < 0
  //   → both roots can be upstream; see the root selection below
  auto const a_coeff = -observationHeight * cos(thetaRad);
  auto const disc = static_pow<2>(injectionHeight) -
                    static_pow<2>(sin(thetaRad) * observationHeight);
  auto const sqrt_disc = [&]() -> LengthType {
    if (disc > 0. * 1_m * 1_m) return sqrt(disc);
    return 0. * 1_m;
  }();

  // Build the backward direction (from core toward injection, opposite to propDir)
  DirectionVector const backwardDir{rootCS,
      {-sin(thetaRad) * cos(phiRad),
       -sin(thetaRad) * sin(phiRad), cos(thetaRad)}};

  // The backward ray from the core meets the injection sphere in two points,
  //   t_near = a_coeff - sqrt_disc   (the intersection CLOSER to the detector)
  //   t_far  = a_coeff + sqrt_disc
  // whose midpoint is the closest-approach (tangent) point t = a_coeff, i.e. the point of
  // minimum distance of the ray from the Earth centre.
  //
  //  * Downgoing showers (cosθ >= 0): a_coeff <= 0, hence t_near < 0 and only t_far is
  //    upstream of the core. The altitude increases monotonically along [0, t_far], so the
  //    far root is the physically meaningful one -- the standard case, unchanged.
  //  * Upgoing showers (cosθ < 0) with the injection sphere BELOW the observation level
  //    (injectionHeight < observationHeight): t_near > 0 as well. This is the
  //    Earth-skimming geometry -- the primary starts t_near in front of the detector and
  //    climbs to the observation level. Because t_near < a_coeff, the whole segment
  //    [0, t_near] lies on the monotonic part of the ray, so its lowest point is the
  //    injection point itself and the flight path never dips below the surface.
  //    Using the far root here instead puts the injection point on the OPPOSITE SIDE OF
  //    THE PLANET: the straight path to the detector then dives up to tens of km below sea
  //    level, the shower develops inside the Earth (the layered atmosphere keeps
  //    extrapolating its innermost layer downwards, so the density grows with depth) and
  //    nothing can reach the detector.
  //  * Upgoing showers with injectionHeight > observationHeight: t_near < 0, so only the
  //    far root exists -- and that path necessarily crosses the Earth (the path check
  //    after the injection validation below warns about it). For such runs put
  //    --injection-height below --observation-level.
  //
  // --injection-root selects the behaviour: "auto" (default) uses the near root whenever it
  // is upstream, "near"/"far" force one of them.
  if (!(disc > 0. * 1_m * 1_m)) {
    Point const tangentPoint = showerCore + backwardDir * a_coeff;
    CORSIKA_LOG_WARN(
        "No real intersection of the shower axis with the injection sphere (radius {} m): "
        "the closest approach of the ray to the Earth centre is r={} m, i.e. {} m above the "
        "requested injection altitude. Falling back to the tangent point, so the injected "
        "altitude will differ from --injection-height.",
        injectionHeight / 1_m, (tangentPoint - center).getNorm() / 1_m,
        ((tangentPoint - center).getNorm() - injectionHeight) / 1_m);
  }

  LengthType const tNear = a_coeff - sqrt_disc;
  LengthType const tFar = a_coeff + sqrt_disc;
  bool const useNear = (injectionRoot == "near") ||
                       (injectionRoot == "auto" && tNear > 0_m);
  LengthType const t = useNear ? tNear : tFar;

  if ((useNear && !(tNear > 0_m)) || (!useNear && !(tFar > 0_m))) {
    CORSIKA_LOG_ERROR(
        "--injection-root={} is not upstream of the detector for this geometry "
        "(t_near={} m, t_far={} m, zenith={} deg). Aborting.",
        injectionRoot, tNear / 1_m, tFar / 1_m, thetaRad * 180.0 / M_PI);
    return EXIT_FAILURE;
  }

  Point const injectionPos = showerCore + backwardDir * t;

  CORSIKA_LOG_INFO(
      "Injection sphere intersections: t_near={} m, t_far={} m -> using the {} one "
      "(t={} m, zenith={} deg); injection altitude={} m",
      tNear / 1_m, tFar / 1_m, useNear ? "NEAR (Earth-skimming)" : "FAR", t / 1_m,
      thetaRad * 180.0 / M_PI,
      ((injectionPos - center).getNorm() - constants::EarthRadius::Mean) / 1_m);
  
  // we make the axis much longer than the inj-core distance since the
  // profile will go beyond the core, depending on zenith angle
  media::ShowerAxis const showerAxis{injectionPos, (showerCore - injectionPos) * 1.2, env};
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
  // 2g-h2n Step D: the WEIGHTED yield (sum of bunch weights = estimated photons at the
  // telescope).  cherenkovHitCount counts bunches, which is not comparable across
  // absorption settings because a surviving bunch keeps its full weight.
  double cherenkovPhotonSum = 0.0;

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
        [eventioWriterPtr, &cherenkovHitCount, &cherenkovPhotonSum](
            double x_m, double y_m, double z_m,
            double lambda_nm, double z_prod_m,
            double dir_x, double dir_y, double dir_z,
            double dir_ground_x, double dir_ground_y, double dir_ground_z,
            double time_ns, double weight, double n_prod, double path_length_m) {
          ++cherenkovHitCount;
          cherenkovPhotonSum += weight;

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
        [cherenkovCsvPtr, &cherenkovHitCount, &cherenkovPhotonSum](
            double x_m, double y_m, double z_m,
            double lambda_nm, double z_prod_m,
            double dir_x, double dir_y, double dir_z,
            double dir_ground_x, double dir_ground_y, double dir_ground_z,
            double time_ns, double weight, double n_prod, double path_length_m) {
          ++cherenkovHitCount;
          cherenkovPhotonSum += weight;
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

  if (projMode != corsika::cherenkov::ProjectionMode::Telescope &&
      telescope_position_correction > 0.0) {
    CORSIKA_LOG_INFO("telescope-position-correction is enabled but projection mode is not telescope; disabling correction.");
    telescope_position_correction = 0.0;
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

  // Telescope observing-plane basis. This is deliberately identical to the basis
  // used inside corsika::cherenkov::CherenkovProcess (see
  // corsika/modules/cherenkov/Cherenkov.inl, computeTelescopeBasis()):
  //   telBasisZ = pointing direction   -> plane normal
  //   telBasisX = (-sin(az), cos(az), 0) -> plane x-axis (parallel to the ground)
  //   telBasisZ x telBasisX              -> plane y-axis (derived where needed)
  DirectionVector const telBasisZ{rootCS,
      {sin(tel_zenith_rad) * cos(tel_azimuth_rad),
       sin(tel_zenith_rad) * sin(tel_azimuth_rad),
       -cos(tel_zenith_rad)}};
  DirectionVector const telBasisX{rootCS,
      {-sin(tel_azimuth_rad), cos(tel_azimuth_rad), 0.}};
  // Telescope centre, used as the origin of all particle output coordinates.
  Point const telCenter{rootCS, tel_x * 1_m, tel_y * 1_m, tel_z * 1_m};
  
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
                        enableDispersion = !disable_dispersion](double h_asl_m, double wavelength_nm) 
                        -> std::pair<double, double> {
    // Compute base (n-1) at the given altitude
    double nminus1_base = 0.0003;  // fallback
    
    if (!atmAlt.empty()) {
      // The library passes HEIGHT ABOVE SEA LEVEL (see Cherenkov.hpp); do not subtract
      // the Earth radius here again.
      double alt_km = h_asl_m / 1000.0;
      
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

  // Vertical refractivity table for the closed-form straight-line timing path.
  // Built from the SAME fine grid the runtime refractive-index lookup uses, so
  // N(z) and A(z) match the lookup exactly (piecewise-exponential in ln N).
  std::shared_ptr<corsika::cherenkov::RefractiveIndexTable> riTable;
  if (!atmAlt_final.empty()) {
    riTable = std::make_shared<corsika::cherenkov::RefractiveIndexTable>(
        atmAlt_final, atmLnNminus1_final, constants::EarthRadius::Mean / 1_m);
    if (riTable->valid()) {
      CORSIKA_LOG_INFO("Refractive-index table: {} nodes, z=[{:.0f}, {:.0f}] m",
                       riTable->size(), riTable->zMin(), riTable->zMax());
    } else {
      CORSIKA_LOG_WARN(
          "Refractive-index table invalid (non-uniform grid?); closed-form timing disabled.");
      riTable.reset();
    }
  }

  corsika::cherenkov::CherenkovProcess cherenkovProcess(
    /* n                 */ 1.0003,
    /* lambdaMin         */ 300.0,
    /* lambdaMax         */ 900.0,
    /* bunchSize         */ 5.0,
    /* groundZ_m         */ groundZ_m_abs,
    /* recorder          */ cherenkovRecorder,
    /* mode              */ projMode,
    /* showerAxis        */ std::make_shared<corsika::media::ShowerAxis>(showerAxis),
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
    /* absorptionFile    */ app["--atmosphere-absorption-file"]->as<std::string>(),
    /* telescopePosCorr  */ telescope_position_correction,
    /* refractiveIndexTable */ riTable,
    /* enableCurvedTable   */ !disable_curved_table);

  cherenkovProcess.setProfiling(enable_cherenkov_profile);
  cherenkovProcess.setTimingSegmentsOverride(cherenkov_timing_segments);
  if (!cher_table_cache.empty()) {
    cherenkovProcess.setTableCacheDir(cher_table_cache);
    CORSIKA_LOG_INFO("Cherenkov table cache: {} (build products reused only when their key "
                     "matches this configuration)",
                     cher_table_cache);
  }
  cherenkovProcess.setGrazingClearance(cher_grazing_clearance);
  CORSIKA_LOG_INFO("Cherenkov grazing rejection band: {:.0f} m (rays nearer their apparent "
                   "horizon are traced exactly, keeping the table tau)",
                   cher_grazing_clearance);
  if (trace_ab_count > 0) {
    // 2g-h2n step A: per-photon A/B of the table route against the plane integrator.
    boost::filesystem::path const abFile = cherDir / "cherenkov_ab.csv";
    cherenkovProcess.enableTraceAB(static_cast<std::size_t>(trace_ab_count),
                                   abFile.string());
    CORSIKA_LOG_INFO("Cherenkov step-A A/B enabled: first {} table-answered photons -> {}",
                     trace_ab_count, abFile.string());
  }
  if (disable_fast_timing) {
    cherenkovProcess.setFastTiming(false);
    CORSIKA_LOG_INFO("Closed-form Cherenkov timing DISABLED (validation reference mode).");
  }
  CORSIKA_LOG_INFO("Cherenkov fast timing: {}",
                   cherenkovProcess.fastTimingEnabled() ? "on" : "off");

  // ---- Phase 3: exact-tracer integrator ---------------------------------------
  {
    using corsika::cherenkov::TracerIntegrator;
    TracerIntegrator mode = TracerIntegrator::RK4;
    if (tracer_integrator == "dopri5") mode = TracerIntegrator::DOPRI5;
    else if (tracer_integrator == "tangent") mode = TracerIntegrator::TANGENT;
    else if (tracer_integrator == "auto") mode = TracerIntegrator::AUTO;
    cherenkovProcess.setTracerIntegrator(mode);
    if (mode == TracerIntegrator::TANGENT) {
      CORSIKA_LOG_WARN(
          "--cherenkov-tracer tangent is a DIAGNOSTIC setting: rays without a turning "
          "point on their path are REFUSED and therefore dropped. Use dopri5 for "
          "production until the Phase-4 mode selector lands.");
    }
    cherenkovProcess.setDopriStepLimits(1.0e-3, tracer_max_step_m);
    cherenkovProcess.setDopriTolerances(tracer_atol_pos_m, tracer_atol_dir,
                                        tracer_atol_opt_m);
    CORSIKA_LOG_INFO(
        "Cherenkov exact tracer: {} (max step {} m, atol pos {} m / dir {} / optical "
        "{} m; dopri5 requires a separable index and a refractive-index table)",
        (mode == TracerIntegrator::DOPRI5)    ? "dopri5 adaptive 5(4)"
        : (mode == TracerIntegrator::TANGENT) ? "tangent only (turning-point coordinates)"
        : (mode == TracerIntegrator::AUTO)    ? "auto (tangent, dopri5 fallback)"
                                              : "rk4 fixed step",
        (tracer_max_step_m > 0.0) ? tracer_max_step_m : -1.0, tracer_atol_pos_m,
        tracer_atol_dir, tracer_atol_opt_m);
  }

  
  std::stringstream args;
  for (int i = 0; i < argc; ++i) { args << argv[i] << " "; }
  // create the output manager that we then register outputs with
  OutputManager output(app["--filename"]->as<std::string>(), seed, args.str(),
                       compressOutput);

  // register energy losses as output
  EnergyLossWriter dEdX{showerAxis, dX};
  output.add("energyloss", dEdX);

  DynamicInteractionProcess<StackType> heModel;

  auto const all_elements = corsika::media::get_all_elements_in_universe(env);
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
  // NOTE: neutrino interactions are handled per shower by
  // corsika::pythia8::SingleUseNeutrinoInteraction (constructed inside the shower
  // loop below). The old corsika::pythia8::NeutrinoInteraction re-initialises
  // Pythia on every call and must not be reused across events.

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

  DirectionVector obsPlaneNormal = DirectionVector(rootCS, {0., 0., 1.});
  DirectionVector obsPlaneXAxis = DirectionVector(rootCS, {1., 0., 0.});
  Point obsPlaneCenter = showerCore;

  // The absorbing "stopping" surface follows the Cherenkov projection mode, so the
  // particle record and the photon record refer to the same detector surface.
  //
  // NOTE ON THE NORMAL SIGN: ObservationPlane treats the positive-normal side as
  // the "catching" side and, in its leap-frog fallback, writes+removes any particle
  // found on the negative side. The normal therefore has to point back towards the
  // incoming shower (i.e. opposite to the shower propagation direction), which is
  // the convention the original ground and near-horizontal branches already used.
  if (projMode == corsika::cherenkov::ProjectionMode::Telescope) {
    // Telescope observing plane: perpendicular to the pointing direction, through
    // the telescope centre -- the same surface used for photon projection.
    // telBasisZ points downstream, so the normal is its negative.
    obsPlaneNormal = DirectionVector(rootCS,
        {-sin(tel_zenith_rad) * cos(tel_azimuth_rad),
         -sin(tel_zenith_rad) * sin(tel_azimuth_rad),
         cos(tel_zenith_rad)});
    obsPlaneXAxis = telBasisX;
    obsPlaneCenter = telCenter;
    CORSIKA_LOG_INFO(
        "Stopping plane = TELESCOPE observing plane (zenith {:.2f} deg, azimuth {:.2f} deg): center=({}, {}, {}), normal=({}, {}, {})",
        tel_zenith_rad * 180.0 / M_PI, tel_azimuth_rad * 180.0 / M_PI,
        obsPlaneCenter.getCoordinates().getX() / 1_m,
        obsPlaneCenter.getCoordinates().getY() / 1_m,
        obsPlaneCenter.getCoordinates().getZ() / 1_m,
        obsPlaneNormal.getComponents()[0], obsPlaneNormal.getComponents()[1],
        obsPlaneNormal.getComponents()[2]);
  } else if (projMode == corsika::cherenkov::ProjectionMode::ShowerAxis) {
    // Shower-coordinate surface: perpendicular to the shower axis, through the core.
    // Normal points back along the shower, towards the incoming shower.
    obsPlaneNormal = DirectionVector(rootCS, {-nx, -ny, -nz});
    // In-plane x-axis: the same construction as Cherenkov.inl (ShowerAxis mode),
    // built from the shower direction (nx, ny, nz).
    double b1x, b1y, b1z;
    if (std::abs(nx) < 0.9) {
      b1x = 0.0;
      b1y = nz;
      b1z = -ny;
    } else {
      b1x = -nz;
      b1y = 0.0;
      b1z = nx;
    }
    double const n1 = std::sqrt(b1x * b1x + b1y * b1y + b1z * b1z);
    obsPlaneXAxis = DirectionVector(rootCS, {b1x / n1, b1y / n1, b1z / n1});
    obsPlaneCenter = showerCore;
    CORSIKA_LOG_INFO("Stopping plane = SHOWER-axis surface (zenith {:.2f} deg)",
                     zenithDeg);
  } else {
    // Ground surface: horizontal plane at the observation level.
    obsPlaneNormal = DirectionVector(rootCS, {0., 0., 1.});
    obsPlaneXAxis = DirectionVector(rootCS, {1., 0., 0.});
    obsPlaneCenter = showerCore;
    CORSIKA_LOG_INFO("Stopping plane = GROUND surface (zenith {:.2f} deg)", zenithDeg);
  }

  Plane const obsPlane(obsPlaneCenter, obsPlaneNormal);
  ObservationPlane<TrackingType, ParticleWriterParquet> observationLevel{
      obsPlane, obsPlaneXAxis,
      true,       // plane should "absorb" particles
      1e-6 * 1_m, // ignored for absorbing planes
      true};     // do not print z-coordinate
  // Record positions/directions in global coordinates with the telescope centre as
  // the origin (instead of plane-local coordinates).
  observationLevel.setGlobalOutput(telCenter);
  // register ground particle output
  output.add("particles", observationLevel);

 CORSIKA_LOG_INFO(
        "Observation Plane Works");
  
  PrimaryWriter<TrackingType, ParticleWriterParquet> primaryWriter(observationLevel);
  output.add("primary", primaryWriter);

CORSIKA_LOG_INFO(
        "Primary Writer Works");

  /* === START: PRIMARY INJECTION OVERRIDE (position + direction) === */
  // Nominal direction/point define the detector geometry (unchanged).
  // The primary may be injected elsewhere and/or pointing elsewhere.
  DirectionVector injDir = propDir;
  bool const directionOverridden =
      (injection_zenith_deg >= 0.0) || (injection_azimuth_deg >= 0.0);
  if (directionOverridden) {
    double const injTheta =
        (injection_zenith_deg >= 0.0 ? injection_zenith_deg
                                     : app["--zenith"]->as<double>()) /
        180. * M_PI;
    double const injPhi =
        (injection_azimuth_deg >= 0.0 ? injection_azimuth_deg
                                      : app["--azimuth"]->as<double>()) /
        180. * M_PI;
    injDir = DirectionVector(rootCS,
                             {sin(injTheta) * cos(injPhi),
                              sin(injTheta) * sin(injPhi), -cos(injTheta)});
    injDir = injDir.normalized();
  }

  // Injection position offsets are applied in GLOBAL horizontal directions, i.e.
  // the root coordinate system's x and y axes. When x/y offsets are given, the
  // --injection-height is interpreted as the RADIAL altitude above sea level
  // (polar coordinate), so the injection z is recomputed to keep the injection
  // point on the sphere of radius R + injection_height:
  //     z = sqrt((R + injection_height)^2 - x^2 - y^2)
  //
  // NOTE: with --zenith 0 the nominal injection point is (0, 0, injectionHeight),
  // so the x/y offsets are the absolute horizontal coordinates of the production
  // point, and its radial distance is exactly R + injection_height.
  DirectionVector const injE1{rootCS, {1., 0., 0.}};   // still used below (diag/YAML)
  DirectionVector const injE2{rootCS, {0., 1., 0.}};
  Point injectionPosPrimary = injectionPos;
  bool const positionOverridden = (injection_x_m != 0.0) || (injection_y_m != 0.0);
  if (positionOverridden) {
    double const xNew = injectionPos.getCoordinates().getX() / 1_m + injection_x_m;
    double const yNew = injectionPos.getCoordinates().getY() / 1_m + injection_y_m;
    double const rNew = injectionHeight / 1_m;            // R + injection_height
    double const rho2 = xNew * xNew + yNew * yNew;
    if (rho2 >= rNew * rNew) {
      CORSIKA_LOG_ERROR(
          "injection-x/y offset too large: horizontal distance {} m is beyond "
          "the injection sphere of radius {} m.", std::sqrt(rho2), rNew);
      return EXIT_FAILURE;
    }
    double const zNew = std::sqrt(rNew * rNew - rho2);    // keep |P| = R + h
    injectionPosPrimary = Point(rootCS, xNew * 1_m, yNew * 1_m, zNew * 1_m);
  }

  /* === validate the primary's starting point against the atmosphere === */
  // Radial distance of the starting point from the Earth centre, expressed as
  // height above sea level. The CORSIKA 7 5-layer atmospheres (USStdBK, ...) are
  // valid for heights up to 112.8 km above sea level; above that the outer layer is
  // of constant density -- see corsika/media/CORSIKA7Atmospheres.hpp.
  // This check is unconditional: it also validates the nominal --injection-height.
  constexpr double kAtmosphereTopAboveSeaLevel_m = 112.8e3;
  double const injRadius_m = (injectionPosPrimary - center).getNorm() / 1_m;
  double const hInj_m = injRadius_m - constants::EarthRadius::Mean / 1_m;

  if (hInj_m < 0.0) {
    CORSIKA_LOG_ERROR(
        "Injection point is BELOW sea level: r={} m, h={} m, "
        "injectionPosPrimary=({}, {}, {}) m. Aborting.",
        injRadius_m, hInj_m,
        injectionPosPrimary.getCoordinates().getX() / 1_m,
        injectionPosPrimary.getCoordinates().getY() / 1_m,
        injectionPosPrimary.getCoordinates().getZ() / 1_m);
    return EXIT_FAILURE;
  }
  if (hInj_m > kAtmosphereTopAboveSeaLevel_m) {
    CORSIKA_LOG_ERROR(
        "Injection point is ABOVE the atmosphere: r={} m, h={} m > {} m "
        "(CORSIKA 7 5-layer atmosphere validity top). Aborting.",
        injRadius_m, hInj_m, kAtmosphereTopAboveSeaLevel_m);
    return EXIT_FAILURE;
  }

  /* === check that the straight flight path to the detector stays above sea level === */
  // The path is a straight line in the root frame from the detector (showerCore) to the
  // injection point, so its minimum distance from the Earth centre is reached at the
  // projection of the centre onto the segment, clamped into the segment. If that point is
  // below sea level the primary has to cross the Earth: the shower then develops inside
  // the planet (the CORSIKA 7 layered atmosphere keeps extrapolating its innermost layer
  // downwards, i.e. the density GROWS with depth), its electromagnetic component is
  // absorbed within metres, only muons leak away, and the detector sees nothing.
  // This is pure geometry and also covers the --injection-x/-y overrides.
  auto const pathVec = injectionPosPrimary - showerCore;              // vector [m]
  double const pathLength_m = pathVec.getNorm() / 1_m;
  double const pathFrac =
      (pathLength_m > 0.0)
          ? std::clamp(pathVec.normalized().dot(center - showerCore) / (1_m * pathLength_m),
                       0.0, 1.0)
          : 0.0;
  Point const pathDeepest = showerCore + pathVec * pathFrac;
  double const pathDeepestAltitude_m =
      (pathDeepest - center).getNorm() / 1_m - constants::EarthRadius::Mean / 1_m;
  if (pathDeepestAltitude_m < 0.0) {
    CORSIKA_LOG_WARN(
        "The straight path from the injection point to the detector dips {} m BELOW sea "
        "level (deepest point {} m from the detector): the primary would have to cross the "
        "Earth and the detector would see nothing. For upgoing showers use "
        "--injection-height below --observation-level, or a larger zenith angle.",
        pathDeepestAltitude_m, pathFrac * pathLength_m);
  }

  if (positionOverridden || directionOverridden) {
    auto const pc = propDir.getComponents();
    auto const ic = injDir.getComponents();
    // NOTE: getComponents()[i] returns a phys::units quantity (implicitly
    // convertible to double). Convert explicitly here, because std::clamp and
    // std::acos require plain doubles (template deduction / overload resolution
    // would otherwise fail).
    double const pcx = pc[0], pcy = pc[1], pcz = pc[2];
    double const icx = ic[0], icy = ic[1], icz = ic[2];
    double const cosang = std::clamp(pcx * icx + pcy * icy + pcz * icz, -1.0, 1.0);
    double const angDeg = std::acos(cosang) * 180.0 / M_PI;

    CORSIKA_LOG_INFO(
        "Primary injection override: position=({}, {}, {}) m, r={} m, h={} m, "
        "direction=({}, {}, {}), angle to shower axis={:.3f} deg",
        injectionPosPrimary.getCoordinates().getX() / 1_m,
        injectionPosPrimary.getCoordinates().getY() / 1_m,
        injectionPosPrimary.getCoordinates().getZ() / 1_m, injRadius_m, hInj_m,
        icx, icy, icz, angDeg);

    auto const denom = obsPlaneNormal.dot(injDir);
    double hit_x_m = 0.0;
    double hit_y_m = 0.0;
    double hit_r_m = 0.0;
    bool crossesPlane = false;
    // NOTE: quantity::magnitude() returns the SIGNED value, so take the absolute
    // value here -- otherwise every downgoing shower (denom = -1) would be
    // misclassified as "parallel to the plane".
    if (std::fabs(denom.magnitude()) < 1e-12) {
      CORSIKA_LOG_WARN(
          "Injection direction is parallel to the observation plane: the primary "
          "will never cross it (in straight-line geometry).");
    } else {
      auto const num = obsPlaneNormal.dot(obsPlaneCenter - injectionPosPrimary);
      double const s_m = (num / denom) / 1_m;
      if (s_m <= 0.0) {
        CORSIKA_LOG_WARN(
            "Primary does not approach the observation plane along its injection "
            "direction (it is injected already past the plane, or moving away from "
            "it): it will never cross it.");
      } else {
        crossesPlane = true;
        Point const hit = injectionPosPrimary + injDir * (s_m * 1_m);
        auto const d = hit - obsPlaneCenter;
        hit_x_m = injE1.dot(d) / 1_m;
        hit_y_m = injE2.dot(d) / 1_m;
        hit_r_m = std::hypot(hit_x_m, hit_y_m);
        CORSIKA_LOG_INFO(
            "Primary ray crosses the observation plane at global horizontal "
            "(x,y)=({}, {}), r={} m, after {} m of flight",
            hit_x_m, hit_y_m, hit_r_m, s_m);
        if (projMode == corsika::cherenkov::ProjectionMode::Telescope) {
          if (hit_r_m > tel_radius) {
            CORSIKA_LOG_WARN(
                "Primary crosses {} m from the telescope centre, OUTSIDE the "
                "{} m aperture: it MISSES the telescope.",
                hit_r_m, tel_radius);
          } else {
            CORSIKA_LOG_INFO(
                "Primary crosses {} m from the telescope centre, inside the "
                "{} m aperture.",
                hit_r_m, tel_radius);
          }
        } else {
          CORSIKA_LOG_INFO(
              "Note: the particle-recording plane is infinite; a non-parallel ray "
              "is always recorded, just at non-zero global horizontal (x,y).");
        }
      }
    }

    boost::filesystem::path injYaml = cherDir / "injection_geometry.yaml";
    std::ofstream injOut(injYaml.string());
    if (injOut) {
      injOut << "earth_radius_m: " << std::setprecision(9)
             << (constants::EarthRadius::Mean / 1_m) << "\n";
      injOut << "nominal:\n";
      injOut << "  direction: [" << pcx << ", " << pcy << ", " << pcz << "]\n";
      injOut << "  injection_point_m: ["
             << injectionPos.getCoordinates().getX() / 1_m << ", "
             << injectionPos.getCoordinates().getY() / 1_m << ", "
             << injectionPos.getCoordinates().getZ() / 1_m << "]\n";
      injOut << "primary:\n";
      injOut << "  direction: [" << icx << ", " << icy << ", " << icz << "]\n";
      injOut << "  injection_point_m: ["
             << injectionPosPrimary.getCoordinates().getX() / 1_m << ", "
             << injectionPosPrimary.getCoordinates().getY() / 1_m << ", "
             << injectionPosPrimary.getCoordinates().getZ() / 1_m << "]\n";
      injOut << "  radius_m: " << injRadius_m << "\n";
      injOut << "  altitude_m: " << hInj_m << "\n";
      injOut << "offset_global_m: {x: " << injection_x_m << ", y: " << injection_y_m
             << "}\n";
      injOut << "angle_to_axis_deg: " << angDeg << "\n";
      injOut << "observation_plane:\n";
      injOut << "  center_m: ["
             << obsPlaneCenter.getCoordinates().getX() / 1_m << ", "
             << obsPlaneCenter.getCoordinates().getY() / 1_m << ", "
             << obsPlaneCenter.getCoordinates().getZ() / 1_m << "]\n";
      // plain doubles, so the YAML holds bare numbers (a phys::units quantity
      // would print with an (empty) unit suffix, e.g. "0 []")
      auto const onc = obsPlaneNormal.getComponents();
      auto const e1c = injE1.getComponents();
      auto const e2c = injE2.getComponents();
      double const onx = onc[0], ony = onc[1], onz = onc[2];
      double const e1x = e1c[0], e1y = e1c[1], e1z = e1c[2];
      double const e2x = e2c[0], e2y = e2c[1], e2z = e2c[2];
      injOut << "  normal: [" << onx << ", " << ony << ", " << onz << "]\n";
      injOut << "  x_axis: [" << e1x << ", " << e1y << ", " << e1z << "]\n";
      injOut << "  y_axis: [" << e2x << ", " << e2y << ", " << e2z << "]\n";
      injOut << "  crosses_plane: " << (crossesPlane ? "true" : "false") << "\n";
      if (crossesPlane) {
        injOut << "  crossing_local_m: {x: " << hit_x_m << ", y: " << hit_y_m
               << ", r: " << hit_r_m << "}\n";
      } else {
        // no crossing: do not emit zeros, which could be mistaken for a crossing
        // exactly at the plane centre
        injOut << "  crossing_local_m: null\n";
      }
      injOut << "  telescope_radius_m: " << tel_radius << "\n";
    } else {
      CORSIKA_LOG_WARN("Failed to write injection geometry YAML: {}", injYaml.string());
    }
  }

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
  // LengthType const step = 1_m;
  CORSIKA_LOG_INFO(
        "Trigger points Works");
  //auto TP =
  //    make_tabulated_flat_atmosphere_radio_propagator(env, injectionPos, surface_, step);
  
  // initiate CoREAS
  //RadioProcess<decltype(detectorCoREAS), CoREAS<decltype(detectorCoREAS), decltype(TP)>,
  //             decltype(TP)>
  //    coreas(detectorCoREAS, TP);

  // register CoREAS with the output manager
  //output.add("CoREAS", coreas);
  // initiate ZHS
  //RadioProcess<decltype(detectorZHS), ZHS<decltype(detectorZHS), decltype(TP)>,
  //             decltype(TP)>
  //    zhs(detectorZHS, TP);
  
  // register ZHS with the output manager
  //output.add("ZHS", zhs);

  // make and register the first interaction writer
  InteractionWriter<setup::Tracking, ParticleWriterParquet> inter_writer(
      showerAxis, observationLevel);
  output.add("interactions", inter_writer);

  /* === 3-D observation recorder (optional) === */
  std::optional<ThreeDPlaneRecorder> obs3d;

  if (observationSteps > 0) {
    // parse --obs3d-particles (comma-separated PDG codes or the preset 'muons')
    std::set<Code> obs3dFilter;
    if (!obs3d_particles.empty()) {
      std::stringstream ss(obs3d_particles);
      std::string tok;
      while (std::getline(ss, tok, ',')) {
        auto const lpos = tok.find_first_not_of(" \t");
        if (lpos == std::string::npos) { continue; }
        auto const rpos = tok.find_last_not_of(" \t");
        tok = tok.substr(lpos, rpos - lpos + 1);
        if (tok.empty()) { continue; }
        if (tok == "muons") {
          obs3dFilter.insert(Code::MuMinus);
          obs3dFilter.insert(Code::MuPlus);
        } else {
          std::istringstream itok(tok);
          int pdgVal = 0;
          itok >> pdgVal;
          if (itok.fail()) {
            CORSIKA_LOG_ERROR("--obs3d-particles: cannot parse '{}' as a PDG code", tok);
            return EXIT_FAILURE;
          }
          obs3dFilter.insert(convert_from_PDG(PDGCode(pdgVal)));
        }
      }
    }

    const int N = observationSteps;

    std::vector<Point> cen;
    std::vector<DirectionVector> nor;
    std::vector<DirectionVector> xa;
    std::vector<DirectionVector> ya;
    std::vector<double> hSlice;

    // All slices are planes parallel to the stopping surface, evenly spaced along
    // its normal, with slice 0 sitting exactly ON the stopping surface. Slice N-1
    // is one step short of the injection point, so no slice coincides with the
    // primary's starting point.
    //
    // The injection point lies *upstream*, i.e. on the negative-normal side for the
    // telescope/shower surfaces, so the stack marches in the direction that reaches
    // the injection point rather than blindly along +normal.
    double const sInj_signed_m =
        (injectionPosPrimary - obsPlaneCenter).dot(obsPlaneNormal) / 1_m;
    double const sInj_m = std::abs(sInj_signed_m);
    if (!(sInj_m > 0.0)) {
      CORSIKA_LOG_ERROR(
          "obs3d: the primary injection point lies IN the stopping plane, so there "
          "is no span for a stack of planes parallel to it (s={} m). This happens "
          "when the nominal shower axis is perpendicular to the telescope pointing "
          "-- e.g. --telescope-zenith 90 with --zenith unset. Give a non-zero "
          "--injection-x/--injection-y, or set --zenith to the shower direction.",
          sInj_signed_m);
      return EXIT_FAILURE;
    }
    double const Delta_m = sInj_m / double(N);
    double const stackSign = (sInj_signed_m >= 0.0) ? 1.0 : -1.0;

    DirectionVector const sliceNormal = obsPlaneNormal.normalized();
    DirectionVector const sliceX = obsPlaneXAxis.normalized();
    DirectionVector const sliceY = sliceNormal.cross(sliceX);

    for (int k = 0; k < N; ++k) {
      double const s_m = stackSign * double(k) * Delta_m;
      Point const ck = obsPlaneCenter + sliceNormal * (s_m * 1_m);
      cen.push_back(ck);
      nor.push_back(sliceNormal);
      xa.push_back(sliceX);
      ya.push_back(sliceY);
      hSlice.push_back(((ck - center).getNorm() - constants::EarthRadius::Mean) / 1_m);
    }

    obs3d.emplace(cen, nor, xa, ya, hSlice, center, constants::EarthRadius::Mean,
                  obs3d_emin, obs3dFilter, telCenter);
    output.add("particles_3d", *obs3d);
    CORSIKA_LOG_INFO(
        "3-D observation recorder active: {} slices parallel to the stopping plane "
        "(spacing={:g} m, altitude span={:g} m)",
        N, Delta_m, hSlice.back() - hSlice.front());
  }

  /* === END: SETUP PROCESS LIST === */

  // trigger the output manager to open the library for writing
  output.startOfLibrary();

  double cascadeSeconds = 0.0; // accumulated EAS.run() wall time (P0 diagnostic)

  // loop over each shower
  for (int i_shower = 1; i_shower < nevent + 1; i_shower++) {

    CORSIKA_LOG_INFO("Shower {} / {} ", i_shower, nevent);

    cherenkovProcess.resetTelescopePositionCorrectionForEvent();

    // randomize the primary energy
    double const eSlope = app["--eslope"]->as<double>();
    PowerLawDistribution<HEPEnergyType> powerLawRng(eSlope, eMin, eMax);
    HEPEnergyType const primaryTotalEnergy =
        (eMax == eMin) ? eMin
                       : powerLawRng(RNGManager<>::getInstance().getRandomStream(
                             "primary_particle"));

    auto const eKin = primaryTotalEnergy - get_mass(beamCode);

    // Energy-dependent Cherenkov bunch size (R1): keep the bunch count roughly
    // constant with primary energy (CORSIKA 7 CERSIZ-like).
    {
      double bunchSize = cherenkov_bunch_size;
      if (!no_bunch_size_auto) {
        double const cap =
            (projMode == corsika::cherenkov::ProjectionMode::Telescope)
                ? std::min(cherenkov_bunch_size_max, 100.0)
                : cherenkov_bunch_size_max;
        double const eTeV = static_cast<double>(primaryTotalEnergy / 1_TeV);
        bunchSize = std::clamp(cherenkov_bunch_size * eTeV, 1.0, cap);
      }
      cherenkovProcess.setBunchSize(bunchSize);
      CORSIKA_LOG_INFO("Cherenkov bunch size: {:.1f} photons/bunch", bunchSize);
    }

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
  // generic per-shower runner: the process sequence is provided by the caller
  auto run_shower = [&](auto& sequence) {
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
      CORSIKA_LOG_INFO("Primary Direction:    {}", injDir.getNorm());
      CORSIKA_LOG_INFO("Point of Injection:   {}", injectionPosPrimary.getCoordinates());
      // NOTE: intentionally the NOMINAL axis length (injectionPos, not
      // injectionPosPrimary): the ShowerAxis in use is anchored on the nominal axis.
      CORSIKA_LOG_INFO("Shower Axis Length:   {}",
                       (showerCore - injectionPos).getNorm() * 1.2);

      // add the desired particle to the stack
      auto const primaryProperties =
          std::make_tuple(beamCode, eKin, injDir.normalized(), injectionPosPrimary, 0_ns);
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
      auto const tRun0 = std::chrono::steady_clock::now();
      EAS.run();
      cascadeSeconds +=
          std::chrono::duration<double>(std::chrono::steady_clock::now() - tRun0).count();

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

      if (telescope_position_correction > 0.0 &&
          projMode == corsika::cherenkov::ProjectionMode::Telescope) {
        auto const correction = cherenkovProcess.getTelescopePositionCorrectionResult();
        boost::filesystem::path corrFile = cherDir / "telescope_position_correction.yaml";
        std::ofstream corrOut(corrFile.string(), std::ios::app);
        if (corrOut) {
          corrOut << "---\n";
          corrOut << "event: " << i_shower << "\n";
          corrOut << "offset_local_m:\n";
          corrOut << "  x: " << correction.offset_local_x_m << "\n";
          corrOut << "  y: " << correction.offset_local_y_m << "\n";
          corrOut << "  z: " << correction.offset_local_z_m << "\n";
          corrOut << "has_update: " << (correction.has_update ? "true" : "false") << "\n";
          corrOut << "last_update_pid: " << correction.last_update_pid << "\n";
        } else {
          CORSIKA_LOG_WARN("Failed to write telescope position correction YAML: {}", corrFile.string());
        }
      }
      if (cherenkovCsvPtr) {
        cherenkovCsvPtr->flush();
      }

  };
  if (neutrinoPrimary) {
    corsika::pythia8::SingleUseNeutrinoInteraction neutrinoPrimaryPythia(
        corsika::setup::C7trackedParticles, NC, CC);
    neutrinoPrimaryPythia.reset();
    if (obs3d) {
      auto sequence = make_sequence(stackInspect, neutrinoPrimaryPythia, hadronSequence,
                                    decaySequence, emCascade, prodprof, emContinuous,
                                    cherenkovProcess, longprof, *obs3d, observationLevel,
                                    inter_writer, thinning, cut);
      run_shower(sequence);
    } else {
      auto sequence = make_sequence(stackInspect, neutrinoPrimaryPythia, hadronSequence,
                                    decaySequence, emCascade, prodprof, emContinuous,
                                    cherenkovProcess, longprof, observationLevel,
                                    inter_writer, thinning, cut);
      run_shower(sequence);
    }
  } else {
    if (obs3d) {
      auto sequence = make_sequence(stackInspect, hadronSequence, decaySequence,
                                    emCascade, prodprof, emContinuous, cherenkovProcess,
                                    longprof, *obs3d, observationLevel, inter_writer,
                                    thinning, cut);
      run_shower(sequence);
    } else {
      auto sequence = make_sequence(stackInspect, hadronSequence, decaySequence,
                                    emCascade, prodprof, emContinuous, cherenkovProcess,
                                    longprof, observationLevel, inter_writer, thinning,
                                    cut);
      run_shower(sequence);
    }
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
  CORSIKA_LOG_INFO("Cascade (EAS.run) time: {:.3f} s ({:.1f}% of total)", cascadeSeconds,
                   100.0 * cascadeSeconds / std::max(1e-9, run_seconds));
  if (enable_cherenkov_profile) {
    auto const& p = cherenkovProcess.getProfile();
    CORSIKA_LOG_INFO(
        "Cherenkov profile: total={:.3f} s, production={:.3f} s, bunchloop={:.3f} s, "
        "timing={:.3f} s, recorder={:.3f} s, bunches={}, recorded={}",
        p.total_s, p.production_s, p.bunchloop_s, p.timing_s, p.recorder_s, p.bunches,
        p.recorded);
  }

  // ---- 2g-h2n Step D: curved-route provenance -----------------------------------
  // Pull-based by design: CherenkovProcess has no finalize hook, so the experiment is
  // the only place these can be reported.  logCurvedTableStats() is the INFO block with
  // the gate/table/rejection split; the summary below is the one-line form that is
  // diffed across configurations.
  if (enable_curvature) {
    cherenkovProcess.logCurvedTableStats();
  }
  {
    auto const s = cherenkovProcess.curvedTableStats();
    std::size_t rejectedTotal = 0;
    for (int k = 0; k < 9; ++k) rejectedTotal += s.rejected[k];
    CORSIKA_LOG_INFO(
        "Cherenkov curved-route summary: tableAvailable={} attempted={} bypassed={} hits={} "
        "(above={} below={}) rejected={} legacyOK={} legacyFail={} "
        "blindLegacyOK={} blindLegacyFail={} dropped={}",
        s.tableAvailable, s.tableAttempted, s.bypassed, s.hits, s.tableHitAbove, s.tableHitBelow,
        rejectedTotal, s.legacyFallbackSucceeded, s.legacyFallbackFailed,
        s.bypassedLegacySucceeded, s.bypassedLegacyFailed, s.finalDrop);
    if (!s.tableAvailable && !s.unavailableReason.empty()) {
        CORSIKA_LOG_WARN("Cherenkov curved table UNAVAILABLE: {}", s.unavailableReason);
    }
    CORSIKA_LOG_INFO(
        "Cherenkov tau provenance: absorption table={} exact chord={} missing={} | "
        "absorption draws={}",
        s.absTableTau, s.absFallbackTau, s.tauMissing,
        cherenkovProcess.absorptionDrawCount());
    if (cherenkovProcess.traceABRows() > 0) {
        CORSIKA_LOG_INFO("Cherenkov step-A A/B rows written: {}",
                         cherenkovProcess.traceABRows());
    }
    auto const pb = cherenkovProcess.planeBackStats();
    if (pb.calls > 0) {
        CORSIKA_LOG_INFO(
            "Cherenkov plane-crossing test (step C): calls={} backward={} "
            "(|tPlane|<=1 m: {}) rejected={} min tPlane={:.2f} m",
            pb.calls, pb.backward, pb.nearZero, pb.rejected, pb.minTPlane);
    }
  }
  // ---- Phase 3: exact-tracer counters ------------------------------------------
  {
    auto const t = cherenkovProcess.tracerStats();
    CORSIKA_LOG_INFO(
        "Cherenkov exact tracer: {} dopri5 / {} tangent / {} rk4 rays | tangent "
        "refused {} ({} of them re-traced by dopri5) | steps {} accepted / {} rejected, "
        "{} derivative evals, h[min/max]={:.4g}/{:.4g} m, capped {}, crossing iters {}, "
        "min-step failures {}",
        t.dopriRays, t.tangentRays, t.rk4Rays, t.tangentUncertified, t.tangentFallbacks,
        t.accepted, t.rejected, t.fEvaluations, t.hMin, t.hMax, t.cappedSteps,
        t.crossingIterations, t.minStepFailures);
  }
  CORSIKA_LOG_INFO("Cherenkov yield: bunches={} photons(sum of weights)={:.1f}",
                   cherenkovHitCount, cherenkovPhotonSum);
  {
    // Deterministic absorption: the recorded yield is the weighted sum, so the
    // keep-floor truncation must be visible as a fraction of the incoming weight.
    auto const& aw = cherenkovProcess.absorptionWeightStats();
    double const cutFrac =
        (aw.emittedWeight > 0.0L)
            ? static_cast<double>(aw.cutWeight / aw.emittedWeight)
            : 0.0;
    CORSIKA_LOG_INFO(
        "Cherenkov absorption (deterministic, weight attenuation): bunches={} cut={} | "
        "weight in={:.1f} out={:.1f} cut={:.1f} (cut fraction={:.3e})",
        aw.evaluatedBunches, aw.cutBunches, static_cast<double>(aw.emittedWeight),
        static_cast<double>(aw.recordedWeight), static_cast<double>(aw.cutWeight), cutFrac);

    // WAVE 5: the tau/survival distribution that actually reached the weight.  This is
    // the number that settles whether the run-level attenuation is consistent with the
    // per-bunch optical depths, or whether the bookkeeping itself is wrong.
    double const tauMean =
        aw.tauCount ? aw.tauSum / static_cast<double>(aw.tauCount) : 0.0;
    double const tauVar = aw.tauCount
                              ? std::max(0.0, aw.tauSumSq / static_cast<double>(aw.tauCount) -
                                                  tauMean * tauMean)
                              : 0.0;
    double const survMean =
        aw.evaluatedBunches ? aw.survSum / static_cast<double>(aw.evaluatedBunches) : 0.0;
    CORSIKA_LOG_INFO(
        "Cherenkov tau/survival: N_tau={} tau_mean={:.4f} tau_rms={:.4f} tau_max={:.2f} | "
        "surv_mean={:.4f} (implied absorbed={:.2f}%) | surv buckets <0.1:{} 0.1-0.5:{} "
        "0.5-0.9:{} >=0.9:{}",
        aw.tauCount, tauMean, std::sqrt(tauVar), aw.tauMax, survMean,
        100.0 * (1.0 - survMean), aw.bSurvLt01, aw.bSurvLt05, aw.bSurvLt09, aw.bSurvHigh);

    // WAVE 5: do the two routes agree on xi and the perigee flag for the same photon?
    CORSIKA_LOG_INFO(
        "Cherenkov xi comparison: N={} mean|dxi|={:.6f} max|dxi|={:.6f} "
        "perigee-flag mismatches={} ({:.1f}%)",
        aw.xiCmpCount, aw.xiCmpCount ? aw.xiDiffSum / static_cast<double>(aw.xiCmpCount) : 0.0,
        aw.xiDiffMax, aw.xiFlagMismatch,
        aw.xiCmpCount ? 100.0 * static_cast<double>(aw.xiFlagMismatch) /
                            static_cast<double>(aw.xiCmpCount)
                      : 0.0);

    // WAVE 5: the gated substitution's tau against the tracer's, on the legacy route.
    CORSIKA_LOG_INFO(
        "Cherenkov tau substitution: N={} tracer_mean={:.4f} table_mean={:.4f} "
        "mean(diff)={:+.4f} max|diff|={:.3f} rejected={}",
        aw.substCount, aw.substCount ? aw.substTracerSum / static_cast<double>(aw.substCount) : 0.0,
        aw.substCount ? aw.substTableSum / static_cast<double>(aw.substCount) : 0.0,
        aw.substCount ? (aw.substTableSum - aw.substTracerSum) / static_cast<double>(aw.substCount)
                      : 0.0,
        aw.substDiffMaxAbs, aw.substRejected);
  }

  return EXIT_SUCCESS;
}
